"""
2~4단계 전체 실행
후보 전수를 생성하고 지표를 계산한 뒤 점수를 매겨 결과를 저장한다.

사용법:
    python run_all.py               # 확장 라이브러리 전체
    python run_all.py --quick       # 축소판 (약 1분, 형식 점검용)
    python run_all.py --legacy168   # 확장 이전 168개 집합 재현 (대조군 자동 포함)
    python run_all.py --controls    # 음성 대조군(prior_3utr)을 조합에 포함
"""
import csv, json, sys, time
from pathlib import Path
import design, metrics, score

RESULTS = Path(__file__).resolve().parent.parent / "results"
CDS_LEN = design.CDS_END - design.CDS_START

METRIC_COLS = ["acc_kozak", "dG_init", "efe_init", "cai",
               "acc_start30", "acc_mean", "gc_cds", "gc_total", "length"]


def main(quick=False, legacy168=False, controls=False):
    t0 = time.time()
    RESULTS.mkdir(exist_ok=True)

    print("=" * 70)
    print("[1/3] 후보 생성")
    # 168개 집합에는 prior_3utr이 포함되어 있으므로 대조군을 함께 투입해야 재현된다
    cands = design.build_candidates(
        methods=["use_best_codon"] if quick else None,
        caps=[0.60] if quick else None,
        keep=design.LEGACY_168 if legacy168 else None,
        include_controls=controls or legacy168,
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
                     "metrics": m})
        if i % 10 == 0 or i == len(cands):
            el = time.time() - t0
            print(f"  {i:4d}/{len(cands)}  경과 {el/60:4.1f}분  "
                  f"남은 예상 {(el/i)*(len(cands)-i)/60:4.1f}분")

    print("\n[3/3] 점수화")
    rows = score.score_all(rows)

    # 저장
    # 실행 조건마다 파일을 따로 둔다 — 기존 순위표를 덮어쓰면 비교가 불가능해진다
    if quick:
        stem = "candidates_quick"
    elif legacy168:
        stem = "candidates_scored_legacy168"
    elif controls:
        stem = "candidates_scored_full_ctrl"
    else:
        stem = "candidates_scored_full"
    csv_path = RESULTS / f"{stem}.csv"
    cols = ["rank", "id", "score", "cds_variant", "utr5", "kozak", "utr3"] + METRIC_COLS
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r["rank"], r["id"], round(r["score"], 4),
                        r["cds_variant"], r["utr5"], r["kozak"], r["utr3"]] +
                       [round(r["metrics"][k], 4) for k in METRIC_COLS])
    with open(RESULTS / f"{stem}.json", "w") as f:
        json.dump(rows, f)

    # 상위 10개 출력
    print(f"\n{'순위':<5}{'후보':<58}{'점수':>7}{'개시부':>8}{'CAI':>7}{'GC':>7}")
    print("-" * 92)
    for r in rows[:10]:
        m = r["metrics"]
        print(f"{r['rank']:<6}{r['id'][:56]:<58}{r['score']:6.3f}"
              f"{m['acc_kozak']:8.3f}{m['cai']:7.3f}{m['gc_cds']:7.3f}")

    print(f"\n저장: {csv_path}")
    print(f"총 {(time.time()-t0)/60:.1f}분")
    return rows


if __name__ == "__main__":
    main(quick="--quick" in sys.argv,
         legacy168="--legacy168" in sys.argv,
         controls="--controls" in sys.argv)
