"""
3단계 — 지표 계산
후보 전수에는 두 가지만 돌린다.
  (1) 개시부 창(-30..+90)의 RNAfold  : 번역 개시 장벽
  (2) 전체 길이 RNAplfold            : 위치별 비결합 확률
전체 길이 RNAfold는 최종 후보에만 쓴다(6분/건).
"""
import math
import RNA
import python_codon_tables as pct

INIT_UP, INIT_DOWN = 30, 90      # AUG 기준 상류/하류 창
PLFOLD_W, PLFOLD_L = 150, 120    # RNAplfold 창 / 최대 염기쌍 거리

# ---- CAI 가중치 테이블 --------------------------------------------------
_table = pct.get_codons_table("h_sapiens")
_W = {}          # 코돈 -> 상대 적응도 w
for aa, codons in _table.items():
    if aa == "*" or len(codons) < 2:   # 종결코돈과 단일코돈 아미노산(M, W) 제외
        continue
    best = max(codons.values())
    for c, f in codons.items():
        _W[c] = (f / best) if best else 0.0


def cai(cds):
    """Codon Adaptation Index — w값들의 기하평균."""
    vals = [math.log(max(_W[cds[i:i+3]], 1e-6))
            for i in range(0, len(cds) - 3, 3)     # 종결코돈 제외
            if cds[i:i+3] in _W]
    return math.exp(sum(vals) / len(vals)) if vals else float("nan")


def gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)


def compute(candidate):
    """후보 하나의 지표를 계산해 dict로 돌려준다."""
    seq = candidate["sequence"]
    s = candidate["cds_start"]                   # AUG의 0-based 위치
    out = {}

    # (1) 개시부 창 — MFE와 앙상블 자유에너지
    w = seq[max(0, s - INIT_UP): s + INIT_DOWN]
    fc = RNA.fold_compound(w)
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    _, efe = fc.pf()
    out["dG_init"] = mfe
    out["efe_init"] = efe

    # (2) 전체 길이 국소 폴딩 — 위치별 비결합 확률
    up = RNA.pfl_fold_up(seq, 1, PLFOLD_W, PLFOLD_L)
    p = [up[i][1] for i in range(1, len(seq) + 1)]     # 1-based -> 0-based
    out["acc_kozak"] = sum(p[max(0, s - 12): s + 3]) / len(p[max(0, s - 12): s + 3])
    out["acc_start30"] = sum(p[s: s + 30]) / 30
    out["acc_mean"] = sum(p) / len(p)

    # (2-1) 구간별 평균 접근성
    # 전체 평균 하나만 쓰면 CDS 4,443 nt에 3'UTR 수백 nt의 차이가 묻혀,
    # 3'UTR만 다른 후보들의 점수가 사실상 같아진다. 구간을 나누어 희석을 막는다.
    e = s + candidate.get("cds_len", 4443)          # CDS 끝(3'UTR 시작)
    def _mean(seg):
        return (sum(seg) / len(seg)) if seg else float("nan")
    out["acc_utr5"] = _mean(p[:s])
    out["acc_cds"] = _mean(p[s:e])
    out["acc_utr3"] = _mean(p[e:])

    # (3) 서열 지표
    cds = seq[s: s + candidate.get("cds_len", 4443)]
    out["gc_cds"] = gc(cds)
    out["gc_total"] = gc(seq)
    out["cai"] = cai(cds)
    out["length"] = len(seq)
    return out
