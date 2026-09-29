"""
1단계 — 관통 골격
후보 3개로 생성 → 지표 → 점수 → 표 출력까지 끝에서 끝까지 한 번 통과시킨다.
층 사이의 데이터 형식을 확정하는 것이 목적이며, 결과의 의미는 아직 따지지 않는다.
"""
import csv, time
from pathlib import Path
import design, metrics, score

RESULTS = Path(__file__).resolve().parent.parent / "results"


def main():
    t0 = time.time()

    print("[1/4] 후보 생성")
    cands = design.build_candidates(
        methods=["use_best_codon"],
        caps=[0.60],
        keep={"mrna1273_5utr", "hba1_5utr", "kozak_used", "prior_3utr"},
    )
    # wildtype + 최적화본 각각에서 2조합 -> 4개 중 3개만 사용
    cands = cands[:3]
    for c in cands:
        c["cds_len"] = 4443
        print(f"   {c['id']}  ({c['length']} nt)")

    print("[2/4] 지표 계산")
    rows = []
    for c in cands:
        t = time.time()
        m = metrics.compute(c)
        rows.append({"id": c["id"], "metrics": m})
        print(f"   {c['id'][:46]:46s} {time.time()-t:4.1f}s")

    print("[3/4] 점수화")
    rows = score.score_all(rows)

    print("[4/4] 표 출력\n")
    cols = ["rank", "id", "score", "acc_kozak", "dG_init", "cai", "acc_mean", "gc_cds"]
    print(f"{'순위':<4}{'후보':<46}{'점수':>7}{'개시부접근':>10}"
          f"{'ΔG개시':>9}{'CAI':>7}{'평균접근':>9}{'GC':>7}")
    print("-" * 99)
    for r in rows:
        m = r["metrics"]
        print(f"{r['rank']:<6}{r['id'][:44]:<46}{r['score']:6.3f}"
              f"{m['acc_kozak']:10.3f}{m['dG_init']:9.1f}{m['cai']:7.3f}"
              f"{m['acc_mean']:9.3f}{m['gc_cds']:7.3f}")

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "skeleton_table.csv"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r["rank"], r["id"], round(r["score"], 4)] +
                       [round(r["metrics"][c], 4) for c in cols[3:]])
    print(f"\n저장: {out}\n총 {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
