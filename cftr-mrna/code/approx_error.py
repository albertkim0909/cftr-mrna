"""
전체 길이 RNAfold(전역 폴딩) vs RNAplfold(창 150) 근사 오차 측정.

3단계 metrics.py가 전체 길이 구조를 RNAplfold 창 150으로 근사한 것이
타당한지 검증한다. 같은 서열에 두 방법을 모두 돌려 위치별 unpaired
확률을 비교하고, 오차의 크기와 방향, 후보 간 순위 보존 여부를 낸다.

모델 설정은 metrics.py와 동일하다.
  - 기본 에너지 모델 (noLP 미사용, dangling end 기본값)
  - 37도
  - RNAplfold: W=150, L=120, u=1

입력:  /home/claude/cands5.json
       [{"rank":1, "id":"...", "cds_start":56, "sequence":"AUG..."}, ...]
       design.py로 후보를 재생성해 만든다. cds_start는 5'UTR+Kozak 길이.

출력:  approx_results.json   후보별 요약 지표
       pos_<rank>.csv        위치별 unpaired 확률과 차이
       struct_<rank>.txt     서열, MFE 구조, centroid 구조

소요:  후보 1개당 약 6분 (4,700 nt 기준). RNAplfold는 3초.

의존:  pip install ViennaRNA   (2.7.2에서 검증)
"""
import json
import math
import time
import RNA

# ─── 설정 (metrics.py와 일치해야 함) ───────────────────────────────
INIT_UP, INIT_DOWN = 30, 90      # 개시창: AUG 기준 -30 ~ +90
PLFOLD_W, PLFOLD_L = 150, 120    # RNAplfold 창 크기, 최대 염기쌍 거리
KT = 0.6163207755                # kcal/mol, 37도

IN_PATH = "/home/claude/cands5.json"
OUT_DIR = "/home/claude"

cands = json.load(open(IN_PATH))
out = []

for c in cands:
    seq = c["sequence"]
    n = len(seq)
    s = c["cds_start"]               # 개시 코돈 A의 0-based 위치
    t0 = time.time()
    print(f"[{c['rank']}위] {n} nt — 전체 길이 RNAfold 시작", flush=True)

    # ─── 전체 길이 RNAfold (거리 제한 없음) ────────────────────────
    fc = RNA.fold_compound(seq)
    mfe_struct, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)       # 긴 서열에서 수치 오버플로 방지
    _, efe = fc.pf()                 # 앙상블 자유에너지
    t_fold = time.time() - t0
    print(f"    MFE {mfe:.2f} / EFE {efe:.2f} / {t_fold:.1f}s", flush=True)

    # 위치별 결합 확률 누적.
    # 전체 확률 행렬(n x n)은 메모리를 많이 쓰므로 희소 리스트로 받는다.
    # 위치 i는 (i,j)와 (j,i) 양쪽에 나타나므로 둘 다 더해야 한다.
    paired = [0.0] * (n + 1)
    for p in fc.plist_from_probs(1e-6):
        if p.type == 0:              # 0 = 염기쌍 확률 (그 외는 다른 종류)
            paired[p.i] += p.p
            paired[p.j] += p.p
    up_fold = [min(1.0, max(0.0, 1.0 - paired[i])) for i in range(1, n + 1)]

    diversity = fc.mean_bp_distance()
    mfe_freq = math.exp((efe - mfe) / KT)
    cent_struct, cent_dist = fc.centroid()

    # ─── RNAplfold (창 150) ───────────────────────────────────────
    t1 = time.time()
    up = RNA.pfl_fold_up(seq, 1, PLFOLD_W, PLFOLD_L)
    up_pl = [up[i][1] for i in range(1, n + 1)]
    t_pl = time.time() - t1
    print(f"    RNAplfold {t_pl:.1f}s (배속 {t_fold/t_pl:.0f}x)", flush=True)

    def win(v, a, b):
        """구간 [a, b)의 평균. metrics.py의 접근성 정의와 동일."""
        seg = v[max(0, a):b]
        return sum(seg) / len(seg)

    rec = {
        "rank": c["rank"], "id": c["id"], "length": n, "cds_start": s,
        "mfe": mfe, "efe": efe, "mfe_freq": mfe_freq,
        "diversity": diversity, "centroid_dist": cent_dist,
        "t_fold": t_fold, "t_plfold": t_pl,
        "acc_kozak_fold": win(up_fold, s - 12, s + 3),
        "acc_kozak_pl":   win(up_pl,   s - 12, s + 3),
        "acc_start30_fold": win(up_fold, s, s + 30),
        "acc_start30_pl":   win(up_pl,   s, s + 30),
        "acc_mean_fold": sum(up_fold) / n,
        "acc_mean_pl":   sum(up_pl) / n,
        "init_window_fold": win(up_fold, s - INIT_UP, s + INIT_DOWN),
        "init_window_pl":   win(up_pl,   s - INIT_UP, s + INIT_DOWN),
    }

    # ─── 오차 지표 ────────────────────────────────────────────────
    # diff > 0 이면 RNAplfold가 더 열려 있다고 판단한 것
    diff = [up_pl[i] - up_fold[i] for i in range(n)]
    rec["mae"] = sum(abs(d) for d in diff) / n          # 어긋남의 크기
    rec["bias"] = sum(diff) / n                          # 어긋남의 방향
    mi = max(range(n), key=lambda i: abs(diff[i]))
    rec["max_err"] = diff[mi]
    rec["max_err_pos"] = mi + 1
    iw = list(range(max(0, s - INIT_UP), s + INIT_DOWN))
    rec["mae_init"] = sum(abs(diff[i]) for i in iw) / len(iw)
    rec["bias_init"] = sum(diff[i] for i in iw) / len(iw)

    mx = sum(up_fold) / n
    my = sum(up_pl) / n
    sxy = sum((up_fold[i] - mx) * (up_pl[i] - my) for i in range(n))
    sxx = sum((up_fold[i] - mx) ** 2 for i in range(n))
    syy = sum((up_pl[i] - my) ** 2 for i in range(n))
    rec["pearson"] = sxy / math.sqrt(sxx * syy)

    out.append(rec)

    # ─── 저장 (후보 하나 끝날 때마다 즉시 기록) ────────────────────
    with open(f"{OUT_DIR}/pos_{c['rank']:03d}.csv", "w") as f:
        f.write("pos,unpaired_rnafold,unpaired_rnaplfold,diff\n")
        for i in range(n):
            f.write(f"{i+1},{up_fold[i]:.6f},{up_pl[i]:.6f},{diff[i]:.6f}\n")
    with open(f"{OUT_DIR}/struct_{c['rank']:03d}.txt", "w") as f:
        f.write(f">{c['id']} rank{c['rank']}\n{seq}\n")
        f.write(f"{mfe_struct} ({mfe:.2f})\n")
        f.write(f"{cent_struct} centroid d={cent_dist:.2f}\n")
    json.dump(out, open(f"{OUT_DIR}/approx_results.json", "w"), indent=1)

    print(f"    MAE {rec['mae']:.4f} / bias {rec['bias']:+.4f} / "
          f"r {rec['pearson']:.4f}\n", flush=True)

print("완료")
