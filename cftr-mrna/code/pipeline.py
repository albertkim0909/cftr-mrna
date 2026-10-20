"""
5단계 — 임의 단백질에 대한 역설계 파이프라인

run_all.py는 참조 CFTR 서열을 전제로 동작한다. 웹에서는 사용자가 넣은 임의의
단백질 서열을 받아야 하므로, CDS 출발점을 만드는 부분만 갈아 끼우고 이후
단계(코돈 최적화 → UTR 조합 → 지표 → 점수)는 기존 모듈을 그대로 부른다.

설계 판단은 바꾸지 않는다. 지표 7개(구간별 접근성 분할), 가중치는 델파이 1차
제출값의 평균, ⑤ 승계 비율 1:1:2, 폴리아데닐화 신호 제약, 1·2위 판정 규칙 1번의
동일 서열 병합까지 results/delphi/final_1120/README.md에 기록된 것을 따른다.
"""
import itertools
import random
import re

import numpy as np
from Bio.Seq import Seq

import design
import metrics
import score

# 점수 체계 — results/delphi/final_1120/README.md 참조
METRIC_SET = "split3utr"
SPLIT_RATIO = (1, 1, 2)          # ⑤ acc_mean 가중치를 5'UTR : CDS : 3'UTR로 승계
MERGE_PATTERN = re.compile(r"match_codon_usage_gc(57|60)\|")


def mean_weights_split():
    """델파이 1차 제출값의 평균을 7지표로 승계한 가중치."""
    out = []
    for w in design_delphi().values():
        nw = {k: w[k] for k in ("acc_kozak", "dG_init", "cai", "acc_start30")}
        t = sum(SPLIT_RATIO)
        nw["acc_utr5"] = w["acc_mean"] * SPLIT_RATIO[0] / t
        nw["acc_cds"] = w["acc_mean"] * SPLIT_RATIO[1] / t
        nw["acc_utr3"] = w["acc_mean"] * SPLIT_RATIO[2] / t
        out.append(nw)
    return score.mean_weights(out)


def design_delphi():
    return score.DELPHI_ROUND1


# ---------------------------------------------------------------------------
# CDS 출발점
# ---------------------------------------------------------------------------
def backtranslate(protein, seed=design.SEED):
    """단백질에서 CDS 출발점을 만든다 — 인간 코돈 사용 빈도에 따른 가중 추출.

    최빈 코돈으로 역번역하면 출발점이 use_best_codon 결과와 같아져 비교 축이
    무너진다. 빈도 분포에서 뽑아 '자연 서열에 준하는' 출발점을 만들고, 씨앗을
    고정해 재현되게 한다. CFTR처럼 실제 참조 서열이 있으면 그쪽을 쓴다.
    """
    rng = random.Random(seed)
    table = metrics._table
    pools = {}
    for aa, codons in table.items():
        items = sorted(codons.items())
        pools[aa] = ([c for c, _ in items], [f for _, f in items])
    out = []
    for aa in protein:
        if aa not in pools:
            raise ValueError(f"표준 아미노산이 아닙니다: {aa}")
        codons, freqs = pools[aa]
        out.append(rng.choices(codons, weights=freqs, k=1)[0])
    stop_codons, stop_freqs = pools["*"]
    out.append(rng.choices(stop_codons, weights=stop_freqs, k=1)[0])
    return "".join(out)


def clean_protein(text):
    """FASTA 머리글과 공백을 걷어내고 대문자 아미노산 서열만 남긴다."""
    lines = [l.strip() for l in text.strip().splitlines()]
    body = "".join(l for l in lines if not l.startswith(">"))
    seq = re.sub(r"[^A-Za-z]", "", body).upper()
    if seq.endswith("*"):
        seq = seq[:-1]
    bad = set(seq) - set("ACDEFGHIKLMNPQRSTVWY")
    if bad:
        raise ValueError(f"표준 아미노산이 아닌 문자: {', '.join(sorted(bad))}")
    if len(seq) < 20:
        raise ValueError("단백질 서열이 너무 짧습니다(20 aa 이상).")
    return seq


# ---------------------------------------------------------------------------
# 후보 생성
# ---------------------------------------------------------------------------
def build_candidates(protein, methods=None, caps=None, require_polya=True,
                     progress=None):
    """임의 단백질에 대해 후보 전수를 조립한다. design.build_candidates의 변형."""
    methods = methods or design.CODON_METHODS
    caps = caps or design.GC_CAPS
    base_cds = backtranslate(protein)
    lib = design.load_utr_library(require_polya=require_polya)

    kozak_motifs = [k["seq"] for k in lib["KOZAK"]]
    for u5 in lib["5UTR"]:
        body, embedded = design.split_embedded_kozak(u5["seq"], kozak_motifs)
        u5["seq"], u5["embedded_kozak"] = body, embedded

    variants = {}
    combos = list(itertools.product(methods, caps))
    for i, (method, cap) in enumerate(combos, 1):
        key = f"{method}_gc{int(cap*100)}"
        if progress:
            progress("optimize", i - 1, len(combos), f"코돈 최적화 {key}")
        variants[key] = design.optimize_cds(base_cds, protein, method, cap)
    variants["backtranslated"] = base_cds   # 최적화하지 않은 출발점(대조)
    if progress:
        progress("optimize", len(combos), len(combos), "코돈 최적화 완료")

    out = []
    for cds_id, cds_seq in variants.items():
        for u5, koz, u3 in itertools.product(lib["5UTR"], lib["KOZAK"], lib["3UTR"]):
            seq = u5["seq"] + koz["seq"] + cds_seq + u3["seq"] + "A" * design.POLY_A
            out.append({
                "id": f"{cds_id}|{u5['id']}|{koz['id']}|{u3['id']}",
                "cds_variant": cds_id, "utr5": u5["id"],
                "kozak": koz["id"], "utr3": u3["id"],
                "cds_start": len(u5["seq"]) + len(koz["seq"]),
                "cds_len": len(cds_seq),
                "length": len(seq), "sequence": seq,
            })
    return out


# ---------------------------------------------------------------------------
# 전체 실행
# ---------------------------------------------------------------------------
def run(protein, require_polya=True, merge_duplicates=True, progress=None,
        cancelled=None):
    """단백질 서열 하나를 받아 후보 전수의 순위표를 돌려준다."""
    cands = build_candidates(protein, require_polya=require_polya,
                             progress=progress)
    rows = []
    n = len(cands)
    for i, c in enumerate(cands, 1):
        if cancelled and cancelled():
            raise RuntimeError("실행이 취소되었습니다.")
        m = metrics.compute(c)
        rows.append({"id": c["id"], "cds_variant": c["cds_variant"],
                     "utr5": c["utr5"], "kozak": c["kozak"], "utr3": c["utr3"],
                     "length": c["length"], "metrics": m})
        if progress and (i % 5 == 0 or i == n):
            progress("metrics", i, n, f"지표 계산 {i}/{n}")

    scored = list(rows)
    if merge_duplicates:
        # 1·2위 판정 규칙 1번 — 동일 서열 병합. match_codon_usage는 GC 상한
        # 0.55/0.58/0.60에서 같은 서열을 내므로 최저 상한판만 남긴다.
        scored = [r for r in scored if not MERGE_PATTERN.match(r["id"])]

    w = mean_weights_split()
    scored = score.score_all(scored, weights=w, metric_set=METRIC_SET)
    return {"protein_len": len(protein),
            "n_generated": n,
            "n_ranked": len(scored),
            "weights": w,
            "rows": scored}
