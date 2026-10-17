"""
2~4단계 전체 실행
후보 전수를 생성하고 지표를 계산한 뒤 점수를 매겨 결과를 저장한다.

사용법:
    python run_all.py            # 전체 (약 13분)
    python run_all.py --quick    # 축소판 (약 1분, 형식 점검용)
"""
import csv, json, sys, time
from pathlib import Path
import design, metrics, score

RESULTS = Path(__file__).resolve().parent.parent / "results"
CDS_LEN = design.CDS_END - design.CDS_START

METRIC_COLS = ["acc_kozak", "dG_init", "efe_init", "cai",
               "acc_start30", "acc_mean", "gc_cds", "gc_total", "length"]


def main(quick=False):
    t0 = time.time()
    RESULTS.mkdir(exist_ok=True)

    print("=" * 70)
    print("[1/3] 후보 생성")
    cands = design.build_candidates(
        methods=["use_best_codon"] if quick else None,
        caps=[0.60] if quick else None,
    )
    if quick:
        cands = cands[::8]
    for c in cands:
        c["cds_len"] = CDS_LEN
    print(f"  후보 {len(cands)}개")

    print(f"\n[2/3] 지표 계산 (후보당 약 4.5초)")
    rows = []
    for i, c in enumerate(cands, 1):
        m = metrics.compute(c)
        rows.append({"id": c["id"], "cds_variant": c["cds_variant"],
                     "utr5": c["utr5"], "kozak": c["kozak"], "utr3": c["utr3"],
                     "is_control": c.get("is_control", False),
                     "metrics": m})
        if i % 10 == 0 or i == len(cands):
            el = time.time() - t0
            print(f"  {i:4d}/{len(cands)}  경과 {el/60:4.1f}분  "
                  f"남은 예상 {(el/i)*(len(cands)-i)/60:4.1f}분")

    print("\n[3/3] 점수화")
    rows = score.score_all(rows)

    # 저장
    csv_path = RESULTS / ("candidates_quick.csv" if quick else "candidates_scored.csv")
    cols = (["rank", "id", "score", "is_control", "would_rank",
             "cds_variant", "utr5", "kozak", "utr3"] + METRIC_COLS)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r["rank"], r["id"], round(r["score"], 4),
                        r["is_control"], r.get("would_rank", ""),
                        r["cds_variant"], r["utr5"], r["kozak"], r["utr3"]] +
                       [round(r["metrics"][k], 4) for k in METRIC_COLS])
    with open(RESULTS / "candidates_scored.json", "w") as f:
        json.dump(rows, f)

    # 상위 10개 출력 (정상 후보만)
    main = [r for r in rows if not r["is_control"]]
    ctrl = [r for r in rows if r["is_control"]]
    print(f"\n{'순위':<5}{'후보':<58}{'점수':>7}{'개시부':>8}{'CAI':>7}{'GC':>7}")
    print("-" * 92)
    for r in main[:10]:
        m = r["metrics"]
        print(f"{r['rank']:<6}{r['id'][:56]:<58}{r['score']:6.3f}"
              f"{m['acc_kozak']:8.3f}{m['cai']:7.3f}{m['gc_cds']:7.3f}")

    # 음성 대조군 — 순위에는 넣지 않고, 끼워 넣었다면 몇 위였을지만 본다
    if ctrl:
        best = ctrl[0]
        print(f"\n음성 대조군 {len(ctrl)}개 (순위 제외)")
        print(f"  최고 점수 {best['score']:.3f} — 정상 후보 {len(main)}개 중 "
              f"{best['would_rank']}위에 해당 ({best['id'][:56]})")
        if best["would_rank"] <= max(5, len(main) // 10):
            print("  경고: 대조군이 상위권에 올라와 있다. "
                  "점수 체계가 나쁜 3'UTR을 걸러내지 못하고 있다는 뜻이다.")

    print(f"\n저장: {csv_path}")
    print(f"총 {(time.time()-t0)/60:.1f}분")
    return rows


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
