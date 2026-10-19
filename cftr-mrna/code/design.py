"""
2단계 — 후보 생성기
참조 CFTR 서열에서 CDS를 추출하고, 코돈 최적화 방식과 GC 상한을 바꿔 가며
CDS 변이체를 만든 뒤, UTR 라이브러리와 조합해 완결된 mRNA 후보를 생성한다.
"""
import csv, itertools, json, random, re, time
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from dnachisel import (DnaOptimizationProblem, EnforceTranslation, EnforceGCContent,
                       AvoidPattern, AvoidHairpins, CodonOptimize)

ROOT = Path(__file__).resolve().parent.parent
DATA, RESULTS = ROOT / "data", ROOT / "results"

CDS_START, CDS_END = 70, 4513          # 0-based slice = 1-based 71..4513
POLY_A = 120

# match_codon_usage는 사용 빈도 분포에서 코돈을 뽑으므로 실행할 때마다 다른 서열이
# 나온다. 같은 후보의 개시부 MFE가 재실행 사이에 최대 12 이상 벌어져 순위가 바뀌었다.
# 씨앗을 고정해 재현되게 한다. use_best_codon과 wildtype은 원래 결정적이다.
SEED = 20251019

# 설계 축 — 여기를 늘리면 후보가 늘어난다
CODON_METHODS = ["use_best_codon", "match_codon_usage"]
GC_CAPS = [0.55, 0.58, 0.60]

# UTR 라이브러리 확장 이전의 조합 집합. 4(5UTR) x 2(KOZAK) x 3(3UTR) x 7(CDS)
# = 168개 후보로, 기존 순위표·델파이 집계가 모두 이 집합 위에서 만들어졌다.
# prior_3utr이 들어 있으므로 include_controls=True와 함께 써야 재현된다.
LEGACY_168 = {
    "hba1_5utr", "hbb_5utr", "mrna1273_5utr", "minimal_5utr",
    "kozak_strong", "kozak_used",
    "hba1_3utr", "hbb_3utr", "prior_3utr",
}


def load_reference():
    rec = SeqIO.read(DATA / "reference" / "cftr_ref_mrna.fasta", "fasta")
    full = str(rec.seq).upper()
    cds = full[CDS_START:CDS_END]
    protein = str(SeqIO.read(DATA / "reference" / "cftr_protein.fasta", "fasta").seq)

    assert len(cds) % 3 == 0, "CDS가 3의 배수가 아님"
    assert cds.startswith("ATG"), "ATG로 시작하지 않음"
    assert str(Seq(cds).translate())[:-1] == protein, "번역 결과가 참조 단백질과 다름"
    assert protein[507] == "F", "508번 잔기가 F가 아님"
    print(f"참조 {len(full)} nt / CDS {len(cds)} nt / 단백질 {len(protein)} aa — 검증 통과")
    return cds, protein


def has_polya_signal(seq):
    """3'UTR에 폴리아데닐화 신호(AATAAA 또는 변이형 ATTAAA)가 있는지.

    위치는 보지 않고 존재 여부만 본다. 생체 내에서 이 신호는 보통 3'UTR
    말단 10~30 nt 이내에 오지만, BNT162b2 유래 서열은 끝에서 76 nt 앞에
    있어 위치 조건을 걸면 실제 임상에서 쓰이는 서열이 탈락한다.
    """
    return re.search(r"A[AT]TAAA", seq) is not None


def load_utr_library(require_polya=True):
    """조립에 쓰는 세 칸(5UTR/KOZAK/3UTR)만 담는다.

    3UTR_CONTROL(음성 대조군)은 전수 조합에서 제외하고 따로 보관한다.
    build_candidates(include_controls=True)로 검증 실행 때만 투입한다.

    서열이 비어 있는 3'UTR 항목(아직 확보하지 못한 자리표시자)은 제외한다.
    빈 3'UTR을 그대로 조합하면 3'UTR이 없는 후보가 오류 없이 만들어져
    결과만 조용히 틀어진다. 5'UTR의 빈 서열(minimal_5utr)은 Kozak만 쓰는
    의도된 음성 대조군이므로 그대로 둔다.

    require_polya=True면 폴리아데닐화 신호가 없는 3'UTR을 후보에서 뺀다.
    구간별 접근성은 분해 효소의 접근 난이도를 재는 지표라 서열이 3'UTR로
    기능할 수 있는지는 재지 못한다. 실제로 제한효소 자리만 모인 대조군
    prior_3utr은 회문·고GC 서열이 단단히 접혀 접근성이 가장 낮게 나오고,
    방향 부호가 음이므로 최고점을 받아 1위에 올랐다. 구조적 안정성과
    기능적 적격성은 축이 다르므로 한 점수에 섞지 않고 제약으로 가른다.
    False로 두면 제약 도입 전 상태를 재현할 수 있다.
    """
    lib = {"5UTR": [], "KOZAK": [], "3UTR": []}
    controls, skipped, no_pas = [], [], []
    with open(DATA / "elements" / "utr_library.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            item = {"id": row["id"], "seq": row["sequence"].strip().upper(),
                    "source": row["source"], "region": row["region"]}
            if row["region"].startswith("3UTR") and not item["seq"]:
                skipped.append(item)
                continue
            if (row["region"].startswith("3UTR") and require_polya
                    and not has_polya_signal(item["seq"])):
                no_pas.append(item)
                continue
            if row["region"] in lib:
                lib[row["region"]].append(item)
            else:
                controls.append(item)
    if skipped:
        print(f"  미확보 제외(빈 3'UTR): {', '.join(s['id'] for s in skipped)}")
    if no_pas:
        print(f"  제약 탈락(폴리아데닐화 신호 없음): "
              f"{', '.join(s['id'] for s in no_pas)}")
    if controls:
        print(f"  기본 조합 제외(대조군): {', '.join(c['id'] for c in controls)}")
    lib["_controls"] = controls
    lib["_no_pas"] = no_pas
    return lib


def split_embedded_kozak(seq, kozak_motifs):
    """5'UTR 끝에 Kozak 서열이 이미 붙어 있으면 떼어내 돌려준다.

    일부 공개 5'UTR(예: tpl_5utr_mod)은 설계 단계에서 Kozak까지 포함해
    발표돼 있다. 그대로 Kozak 칸을 덧붙이면 GCCACC가 두 번 들어가
    AUG 앞 간격이 어긋나므로, 떼어낸 뒤 Kozak 칸이 공급하도록 한다.
    가장 긴 모티프부터 맞춰 본다.
    """
    for motif in sorted(kozak_motifs, key=len, reverse=True):
        if motif and seq.endswith(motif):
            return seq[: -len(motif)], motif
    return seq, None


def optimize_cds(cds, protein, method, gc_cap):
    """코돈 최적화 — 아미노산 서열 보존을 제약으로 강제한다."""
    # 조합마다 같은 씨앗에서 출발시켜, 앞선 조합의 난수 소비가 뒤 결과를
    # 바꾸지 않도록 한다(설계 축을 늘리거나 줄여도 각 서열이 그대로 나온다).
    random.seed(SEED)
    np.random.seed(SEED)
    problem = DnaOptimizationProblem(
        sequence=cds,
        constraints=[
            EnforceTranslation(),
            EnforceGCContent(mini=0.40, maxi=gc_cap),
            EnforceGCContent(mini=0.25, maxi=0.75, window=50),
            AvoidPattern("AAAAAAAAA"), AvoidPattern("TTTTTTTTT"),
            AvoidPattern("GGGGGGGGG"), AvoidPattern("CCCCCCCCC"),
            AvoidHairpins(stem_size=10, hairpin_window=60),
        ],
        objectives=[CodonOptimize(species="h_sapiens", method=method)],
        logger=None,
    )
    problem.resolve_constraints()
    problem.optimize()
    opt = problem.sequence
    assert str(Seq(opt).translate())[:-1] == protein, f"{method}/{gc_cap}: 아미노산이 바뀜"
    return opt


def build_candidates(methods=None, caps=None, keep=None, include_controls=False,
                     require_polya=True):
    """methods/caps를 좁히거나 keep(허용 id 집합)을 주면 후보 수를 줄일 수 있다.

    include_controls=True면 3UTR_CONTROL(음성 대조군)을 3'UTR 칸에 함께 넣어
    조합한다. 점수 체계가 기능 없는 3'UTR을 걸러 내는지 검증할 때 쓴다.
    단 require_polya=True면 대조군은 그 제약에서 먼저 탈락하므로, 제약 도입
    전 상태를 재현하려면 require_polya=False와 함께 써야 한다.
    """
    methods = methods or CODON_METHODS
    caps = caps or GC_CAPS
    cds, protein = load_reference()
    lib = load_utr_library(require_polya=require_polya)
    if include_controls and lib["_controls"]:
        lib["3UTR"] = lib["3UTR"] + lib["_controls"]
        print(f"  대조군 투입: {', '.join(c['id'] for c in lib['_controls'])}")
    if keep:
        for region in ("5UTR", "KOZAK", "3UTR"):
            lib[region] = [x for x in lib[region] if x["id"] in keep] or lib[region]

    # 0) 5'UTR에 이미 박혀 있는 Kozak 떼어내기 — 중복 삽입 방지
    kozak_motifs = [k["seq"] for k in lib["KOZAK"]]
    for u5 in lib["5UTR"]:
        body, embedded = split_embedded_kozak(u5["seq"], kozak_motifs)
        u5["seq"], u5["embedded_kozak"] = body, embedded
        if embedded:
            print(f"  {u5['id']}: 끝에 Kozak({embedded}) 내장 — 떼어내고 {len(body)} nt로 사용")

    # 1) CDS 변이체
    variants = {}
    for method, cap in itertools.product(methods, caps):
        key = f"{method}_gc{int(cap*100)}"
        t = time.time()
        variants[key] = optimize_cds(cds, protein, method, cap)
        gc = (variants[key].count("G") + variants[key].count("C")) / len(variants[key])
        changed = sum(a != b for a, b in zip(cds, variants[key])) / len(cds)
        print(f"  {key:32s} GC {gc:.3f}  변경 {changed:5.1%}  {time.time()-t:4.1f}s")
    variants["wildtype"] = cds          # 대조군: 최적화하지 않은 원본

    # 2) UTR 조합
    candidates = []
    for cds_id, cds_seq in variants.items():
        for u5, koz, u3 in itertools.product(lib["5UTR"], lib["KOZAK"], lib["3UTR"]):
            seq = u5["seq"] + koz["seq"] + cds_seq + u3["seq"] + "A" * POLY_A
            # AUG 바로 앞에 Kozak이 한 번만 오는지 — 5'UTR 본체가 또 Kozak으로
            # 끝나면 연달아 두 번 들어간다. 서열 내부에 우연히 나타나는
            # GCCACC는 AUG와 떨어져 있어 문제가 아니므로 접합부만 본다.
            assert not any(u5["seq"].endswith(m) for m in kozak_motifs), \
                f"{u5['id']}|{koz['id']}: AUG 앞에 Kozak이 연달아 두 번 삽입됨"
            candidates.append({
                "id": f"{cds_id}|{u5['id']}|{koz['id']}|{u3['id']}",
                "cds_variant": cds_id, "utr5": u5["id"],
                "kozak": koz["id"], "utr3": u3["id"],
                "utr5_embedded_kozak": u5.get("embedded_kozak"),
                "cds_start": len(u5["seq"]) + len(koz["seq"]),   # AUG의 0-based 위치
                "length": len(seq), "sequence": seq,
            })
    return candidates


if __name__ == "__main__":
    RESULTS.mkdir(exist_ok=True)
    cands = build_candidates()
    print(f"\n후보 {len(cands)}개 생성 (길이 {min(c['length'] for c in cands)}"
          f"~{max(c['length'] for c in cands)} nt)")
    with open(RESULTS / "candidates.json", "w") as f:
        json.dump(cands, f)
    print(f"저장: {RESULTS / 'candidates.json'}")
