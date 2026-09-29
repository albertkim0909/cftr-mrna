"""
2단계 — 후보 생성기
참조 CFTR 서열에서 CDS를 추출하고, 코돈 최적화 방식과 GC 상한을 바꿔 가며
CDS 변이체를 만든 뒤, UTR 라이브러리와 조합해 완결된 mRNA 후보를 생성한다.
"""
import csv, itertools, json, time
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from dnachisel import (DnaOptimizationProblem, EnforceTranslation, EnforceGCContent,
                       AvoidPattern, AvoidHairpins, CodonOptimize)

ROOT = Path(__file__).resolve().parent.parent
DATA, RESULTS = ROOT / "data", ROOT / "results"

CDS_START, CDS_END = 70, 4513          # 0-based slice = 1-based 71..4513
POLY_A = 120

# 설계 축 — 여기를 늘리면 후보가 늘어난다
CODON_METHODS = ["use_best_codon", "match_codon_usage"]
GC_CAPS = [0.55, 0.58, 0.60]


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


def load_utr_library():
    lib = {"5UTR": [], "KOZAK": [], "3UTR": []}
    with open(DATA / "elements" / "utr_library.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            lib[row["region"]].append({"id": row["id"],
                                       "seq": row["sequence"].strip().upper(),
                                       "source": row["source"]})
    return lib


def optimize_cds(cds, protein, method, gc_cap):
    """코돈 최적화 — 아미노산 서열 보존을 제약으로 강제한다."""
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


def build_candidates(methods=None, caps=None, keep=None):
    """methods/caps를 좁히거나 keep(허용 id 집합)을 주면 후보 수를 줄일 수 있다."""
    methods = methods or CODON_METHODS
    caps = caps or GC_CAPS
    cds, protein = load_reference()
    lib = load_utr_library()
    if keep:
        for region in lib:
            lib[region] = [x for x in lib[region] if x["id"] in keep] or lib[region]

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
            candidates.append({
                "id": f"{cds_id}|{u5['id']}|{koz['id']}|{u3['id']}",
                "cds_variant": cds_id, "utr5": u5["id"],
                "kozak": koz["id"], "utr3": u3["id"],
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
