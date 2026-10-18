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
        include_controls=True,      # prior_3utr(음성 대조군)을 전수 조합에 포함
    )
    if quick:
        cands = cands[::8]
    for c in cands:
        c["cds_len"] = CDS_LEN
    print(f"  후보 {len(cands)}개")

    # 이어 달리기 — 계산이 끝난 후보는 cache에 한 줄씩 적어 두고 재실행 때 건너뛴다.
    # 난수 시드가 고정돼 있어 같은 id는 같은 서열을 뜻한다.
    cache_path = RESULTS / ("metrics_cache_quick.jsonl" if quick else "metrics_cache.jsonl")
    done = {}
    if cache_path.exists():
        with open(cache_path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    done[r["id"]] = r
        print(f"  이어 달리기: 기존 {len(done)}건 재사용")

    print(f"\n[2/3] 지표 계산 (후보당 약 2.6초)")
    rows = []
    cache_f = open(cache_path, "a", encoding="utf-8")
    for i, c in enumerate(cands, 1):
        if c["id"] in done:
            rows.append(done[c["id"]])
            continue
        m = metrics.compute(c)
        row = {"id": c["id"], "cds_variant": c["cds_variant"],
               "utr5": c["utr5"], "kozak": c["kozak"], "utr3": c["utr3"],
               "metrics": m}
        rows.append(row)
        cache_f.write(json.dumps(row) + "\n")
        cache_f.flush()
        if i % 10 == 0 or i == len(cands):
            el = time.time() - t0
            print(f"  {i:4d}/{len(cands)}  경과 {el/60:4.1f}분  "
                  f"남은 예상 {(el/i)*(len(cands)-i)/60:4.1f}분")
    cache_f.close()

    print("\n[3/3] 점수화")
    rows = score.score_all(rows)

    # 저장
    csv_path = RESULTS / ("candidates_quick.csv" if quick else "candidates_scored.csv")
    cols = ["rank", "id", "score", "cds_variant", "utr5", "kozak", "utr3"] + METRIC_COLS
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r["rank"], r["id"], round(r["score"], 4),
                        r["cds_variant"], r["utr5"], r["kozak"], r["utr3"]] +
                       [round(r["metrics"][k], 4) for k in METRIC_COLS])
    with open(RESULTS / "candidates_scored.json", "w") as f:
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
    main(quick="--quick" in sys.argv)
