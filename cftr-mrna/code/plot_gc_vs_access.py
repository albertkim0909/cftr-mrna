"""
1주차 산출물 — 코돈 최적화 강도(GC)와 번역 개시부 접근성의 관계
가설: 코돈을 강하게 최적화할수록 GC가 올라가 구조가 단단해지고,
      그 결과 리보솜이 진입해야 할 개시부가 닫힌다.
"""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

# 검증된 3색 (all-pairs 통과)
COLORS = {"use_best_codon": "#2a78d6", "match_codon_usage": "#eb6834", "wildtype": "#1baf7a"}
LABELS = {"use_best_codon": "최빈 코돈 사용", "match_codon_usage": "코돈 분포 일치", "wildtype": "야생형 (비최적화)"}
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e5e4e0", "#fcfcfb"


def group(cds_variant):
    if cds_variant == "wildtype":
        return "wildtype"
    return "match_codon_usage" if "match_codon_usage" in cds_variant else "use_best_codon"


def main():
    rows = list(csv.DictReader(open(RESULTS / "candidates_scored.csv", encoding="utf-8")))
    for r in rows:
        r["g"] = group(r["cds_variant"])

    plt.rcParams["font.family"] = ["NanumGothic","Noto Sans CJK KR","DejaVu Sans"]
    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=200)
    fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)

    for g in ["use_best_codon", "match_codon_usage", "wildtype"]:
        pts = [r for r in rows if r["g"] == g]
        ax.scatter([float(r["gc_cds"]) * 100 for r in pts],
                   [float(r["acc_kozak"]) for r in pts],
                   s=46, c=COLORS[g], label=LABELS[g],
                   edgecolors=SURFACE, linewidths=1.6, zorder=3, alpha=0.92)

    ax.set_xlabel("CDS의 GC 함량 (%)", color=INK2, fontsize=10)
    ax.set_ylabel("개시부 비결합 확률", color=INK2, fontsize=10)
    ax.set_title("코돈 최적화가 올린 GC는 번역 개시부를 닫는다",
                 color=INK, fontsize=12.5, pad=12, loc="left")
    ax.grid(True, color=GRID, linewidth=0.9, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)

    leg = ax.legend(frameon=False, fontsize=9.5, loc="upper right")
    for t in leg.get_texts():
        t.set_color(INK)

    fig.tight_layout()
    out = RESULTS / "fig1_gc_vs_accessibility.png"
    fig.savefig(out, facecolor=SURFACE)
    print("저장:", out)

    # 표 형태로도 남긴다 (색만으로 정보를 전달하지 않기 위함)
    with open(RESULTS / "fig1_data.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["group", "gc_cds_percent", "acc_kozak", "cai", "id"])
        for r in rows:
            w.writerow([LABELS[r["g"]], round(float(r["gc_cds"]) * 100, 2),
                        r["acc_kozak"], r["cai"], r["id"]])


if __name__ == "__main__":
    main()
