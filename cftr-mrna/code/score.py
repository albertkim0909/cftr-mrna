"""
4단계 — 점수화·순위·결론 안정성 검증

지표를 후보 집단 안에서 0~1로 정규화한 뒤 가중합해 순위를 낸다.

1-2의 목적은 "정답 가중치를 찾는 것"이 아니라 "가중치를 누구 기준으로 잡든
상위 후보가 유지되는가"를 확인하는 것이다. 따라서 이 모듈은 가중치 한 벌로
순위를 내는 데 그치지 않고, 여러 벌을 넣어 상위권이 유지되는지까지 본다.

방향 부호는 가중치와 달리 사람마다 달라서는 안 되는 모형 구조 항목이므로
METRIC_SETS에 고정해 둔다. 투표 대상이 아니라 합의 대상이다.

사용 예:
    import score
    score.score_all(rows)                       # 기본 지표·기본 가중치
    score.score_all(rows, metric_set="split3utr", weights=w)
    rep = score.stability_report(rows, WEIGHT_VECTORS)
"""
import itertools
import math
import random

# ---------------------------------------------------------------------------
# 지표 집합 — (지표명, 방향). +1은 클수록 좋음, -1은 작을수록 좋음
# 지표 개수가 바뀌어도 코드를 고치지 않도록 이름을 바깥에서 받는다.
# ---------------------------------------------------------------------------
METRIC_SETS = {
    # 현행 5지표. 델파이 1차가 이 집합을 전제로 응답을 받았다.
    "base": [
        ("acc_kozak",   +1),   # 개시부가 열려 있을수록 리보솜 진입이 쉽다
        ("dG_init",     +1),   # 개시부 MFE가 높을수록(덜 음수) 장벽이 낮다
        ("cai",         +1),   # 코돈이 인간에 적합할수록 신장이 고르다
        ("acc_start30", +1),   # 개시 직후 구간의 접근성
        ("acc_mean",    -1),   # 전체가 너무 열려 있으면 분해에 취약하다
    ],
    # acc_mean을 구간별로 쪼개 3'UTR 기여가 CDS 길이에 묻히지 않게 한 집합.
    "split3utr": [
        ("acc_kozak",   +1),
        ("dG_init",     +1),
        ("cai",         +1),
        ("acc_start30", +1),
        ("acc_utr5",    -1),
        ("acc_cds",     -1),
        ("acc_utr3",    -1),
    ],
    # ④ 개시 직후 접근성을 뺀 집합(1-3에서 삭제 확정 시 사용).
    "base_no_start30": [
        ("acc_kozak",   +1),
        ("dG_init",     +1),
        ("cai",         +1),
        ("acc_mean",    -1),
    ],
}

DEFAULT_METRIC_SET = "base"

# 김도윤이 정한 잠정 가중치 — 중립 기준이 아니므로 비교 기준으로 쓰지 않는다.
PROVISIONAL_WEIGHTS = {
    "acc_kozak": 0.30, "dG_init": 0.20, "cai": 0.25,
    "acc_start30": 0.10, "acc_mean": 0.15,
}

# 델파이 1차 응답(각자 제출값). 방향 부호는 위 METRIC_SETS가 고정한다.
DELPHI_ROUND1 = {
    "양우진": {"acc_kozak": 0.30, "dG_init": 0.25, "cai": 0.15,
               "acc_start30": 0.10, "acc_mean": 0.20},
    "윤준우": {"acc_kozak": 0.30, "dG_init": 0.15, "cai": 0.35,
               "acc_start30": 0.00, "acc_mean": 0.20},
    "이지유": {"acc_kozak": 0.35, "dG_init": 0.15, "cai": 0.25,
               "acc_start30": 0.15, "acc_mean": 0.10},
    "김도윤": {"acc_kozak": 0.30, "dG_init": 0.20, "cai": 0.25,
               "acc_start30": 0.10, "acc_mean": 0.15},
}

TOP_K = 5             # 판정은 상위 5개 기준으로만 한다
STABILITY_CUT = 0.80  # 무작위 가중치에서 상위 5위 진입 비율이 이 값 이상이면 안정
N_RANDOM = 2000       # 무작위 가중치 추출 벌 수


# ---------------------------------------------------------------------------
# 가중치 처리
# ---------------------------------------------------------------------------
def mean_weights(vectors):
    """여러 응답의 지표별 평균. 합이 1에서 벗어나면 재정규화한다.

    방향 부호는 METRIC_SETS가 고정하므로 여기서는 크기만 평균한다.
    부호를 섞어 평균하면 합이 1에서 벗어나 재정규화가 필요해지고,
    그 보정을 '반올림 보정'으로 적는 잘못이 생긴다.
    """
    keys = list(next(iter(vectors)).keys()) if isinstance(vectors, list) else None
    if keys is None:
        vectors = list(vectors)
        keys = list(vectors[0].keys())
    out = {k: sum(abs(v[k]) for v in vectors) / len(vectors) for k in keys}
    return normalize_weights(out)


def normalize_weights(w):
    s = sum(w.values())
    if s <= 0:
        raise ValueError("가중치 합이 0 이하입니다")
    return {k: v / s for k, v in w.items()}


def random_weights(keys, n=N_RANDOM, seed=0):
    """합이 1인 공간에서 치우침 없이 n벌 뽑는다(Dirichlet(1,...,1)).

    지수분포 난수를 지표 수만큼 뽑아 합으로 나누면 단체(simplex) 위의
    균등분포가 된다. 네 사람 응답만으로는 "가중치에 좌우되지 않는다"를
    주장하기에 표본이 너무 작아, 가능한 가중치를 전부 훑어보는 장치다.
    """
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        e = [rng.expovariate(1.0) for _ in keys]
        t = sum(e)
        out.append({k: v / t for k, v in zip(keys, e)})
    return out


def perturb(w, pct=0.20, seed=0):
    """가중치를 ±pct 흔든 벌을 돌려준다(민감도 분석용)."""
    rng = random.Random(seed)
    out = {k: v * (1 + rng.uniform(-pct, pct)) for k, v in w.items()}
    return normalize_weights(out)


# ---------------------------------------------------------------------------
# 정규화와 점수
# ---------------------------------------------------------------------------
def normalize(values, how="minmax"):
    if how == "minmax":
        lo, hi = min(values), max(values)
        if hi - lo < 1e-12:
            return [0.5] * len(values)
        return [(v - lo) / (hi - lo) for v in values]
    if how == "zscore":
        n = len(values)
        mu = sum(values) / n
        sd = math.sqrt(sum((v - mu) ** 2 for v in values) / (n - 1)) if n > 1 else 0.0
        if sd < 1e-12:
            return [0.0] * n
        return [(v - mu) / sd for v in values]
    raise ValueError(f"알 수 없는 정규화 방식: {how}")


def build_norm(rows, metric_set=DEFAULT_METRIC_SET, how="minmax"):
    """후보별 정규화 지표값 행렬을 만든다. 가중치와 무관하므로 한 번만 계산한다."""
    spec = METRIC_SETS[metric_set]
    norm = {}
    for key, direction in spec:
        vals = [r["metrics"][key] for r in rows]
        nv = normalize(vals, how)
        if direction < 0:
            nv = [(1 - v) if how == "minmax" else -v for v in nv]
        for r, v in zip(rows, nv):
            norm.setdefault(r["id"], {})[key] = v
    return norm


def rank_with(rows, weights, norm=None, metric_set=DEFAULT_METRIC_SET, how="minmax"):
    """가중치 한 벌로 순위표를 낸다. [(후보 id, 점수)] 내림차순."""
    spec = METRIC_SETS[metric_set]
    norm = norm if norm is not None else build_norm(rows, metric_set, how)
    w = normalize_weights({k: weights.get(k, 0.0) for k, _ in spec})
    scored = [(r["id"], sum(w[k] * norm[r["id"]][k] for k, _ in spec)) for r in rows]
    scored.sort(key=lambda t: -t[1])
    return scored


def score_all(rows, weights=None, metric_set=DEFAULT_METRIC_SET, how="minmax"):
    """기존 호출 호환 — rows에 norm·score·rank를 채우고 정렬해 돌려준다."""
    spec = METRIC_SETS[metric_set]
    w = normalize_weights({k: (weights or PROVISIONAL_WEIGHTS).get(k, 0.0)
                           for k, _ in spec})
    norm = build_norm(rows, metric_set, how)
    for r in rows:
        r["norm"] = norm[r["id"]]
        r["score"] = sum(w[k] * r["norm"][k] for k, _ in spec)
    rows.sort(key=lambda r: -r["score"])
    for rank, r in enumerate(rows, 1):
        r["rank"] = rank
    return rows


# ---------------------------------------------------------------------------
# 일치도 통계
# ---------------------------------------------------------------------------
def _ranks(scored):
    return {cid: i + 1 for i, (cid, _) in enumerate(scored)}


def spearman(a, b):
    """두 순위표의 스피어만 ρ. 동점이 없다고 보고 단순식을 쓴다."""
    keys = [k for k in a if k in b]
    n = len(keys)
    d2 = sum((a[k] - b[k]) ** 2 for k in keys)
    return 1 - 6 * d2 / (n * (n * n - 1))


def kendall_tau(a, b):
    keys = [k for k in a if k in b]
    con = dis = 0
    for x, y in itertools.combinations(keys, 2):
        s = (a[x] - a[y]) * (b[x] - b[y])
        if s > 0:
            con += 1
        elif s < 0:
            dis += 1
    return (con - dis) / (con + dis) if (con + dis) else float("nan")


def kendall_w(rank_tables):
    """Kendall's W — 응답자 m명이 후보 n개를 매긴 순위의 일치도(0~1)."""
    m = len(rank_tables)
    keys = sorted(set.intersection(*[set(t) for t in rank_tables]))
    n = len(keys)
    if m < 2 or n < 2:
        return float("nan")
    tot = {k: sum(t[k] for t in rank_tables) for k in keys}
    mean = sum(tot.values()) / n
    s = sum((v - mean) ** 2 for v in tot.values())
    return 12 * s / (m * m * (n ** 3 - n))


# ---------------------------------------------------------------------------
# 결론 안정성 보고
# ---------------------------------------------------------------------------
def verdict(overlap, mean_rank, top_k=TOP_K):
    """판정 기준 — 집계 전에 확정해 둔 것. 결과를 보고 바꾸지 않는다."""
    four = [c for c, n in overlap.items() if n >= 4]
    three = [c for c, n in overlap.items() if n == 3]
    for c in four:
        if mean_rank.get(c, 99) <= 3:
            return ("가중치 무관", f"{c}가 전원 상위 {top_k}이며 평균 가중치 순위 "
                                   f"{mean_rank[c]}위 — 결론이 가중치 설정에 좌우되지 않음")
    if len(three) >= 2:
        return ("중간", f"3명 공통 후보 {len(three)}개({', '.join(three)}) — "
                        "소수 후보군으로 좁혀지나 단일 최적은 확정되지 않음")
    return ("조건부", "단일 최적 서열을 제시하지 않고, 무엇을 중시하느냐에 따라 "
                      "추천이 갈린다는 조건부 결론")


def stability_report(rows, weight_vectors, metric_set=DEFAULT_METRIC_SET,
                     how="minmax", top_k=TOP_K, n_random=N_RANDOM, seed=0):
    """응답자별 순위표, 겹침, 상관, Kendall's W, 무작위 가중치 안정도를 한 번에 낸다."""
    spec = METRIC_SETS[metric_set]
    keys = [k for k, _ in spec]
    norm = build_norm(rows, metric_set, how)

    tables = {who: rank_with(rows, w, norm, metric_set, how)
              for who, w in weight_vectors.items()}
    tables["평균"] = rank_with(rows, mean_weights(list(weight_vectors.values())),
                              norm, metric_set, how)
    ranks = {who: _ranks(t) for who, t in tables.items()}

    tops = {who: [c for c, _ in t[:top_k]] for who, t in tables.items()}
    people = [w for w in weight_vectors]

    # 겹침은 응답자 본인 표로만 센다. 평균 표는 네 사람 가중치에서 산술적으로
    # 파생된 값이라 독립된 한 표가 아니며, 이를 세면 다수 쪽 견해가 두 번
    # 계산되어 "4명 전원"이라는 판정 기준이 사실과 어긋나게 된다.
    overlap = {}
    for who in people:
        for c in tops[who]:
            overlap[c] = overlap.get(c, 0) + 1
    overlap = dict(sorted(overlap.items(), key=lambda kv: -kv[1]))

    pairs = {f"{a}-{b}": {"spearman": spearman(ranks[a], ranks[b]),
                          "kendall": kendall_tau(ranks[a], ranks[b]),
                          "top_overlap": len(set(tops[a]) & set(tops[b]))}
             for a, b in itertools.combinations(people, 2)}
    w_stat = kendall_w([ranks[p] for p in people])

    hits = {r["id"]: 0 for r in rows}
    for w in random_weights(keys, n_random, seed):
        for c, _ in rank_with(rows, w, norm, metric_set, how)[:top_k]:
            hits[c] += 1
    stability = {c: n / n_random for c, n in sorted(hits.items(), key=lambda kv: -kv[1])
                 if n}
    stable = [c for c, p in stability.items() if p >= STABILITY_CUT]

    kind, why = verdict(overlap, _ranks(tables["평균"]), top_k)
    return {"metric_set": metric_set, "normalization": how, "top_k": top_k,
            "tables": tables, "ranks": ranks, "top": tops, "overlap": overlap,
            "pairwise": pairs, "kendall_w": w_stat,
            "stability": stability, "stable": stable,
            "verdict": kind, "verdict_reason": why}


def sensitivity(rows, weights, pct=0.20, n=200, metric_set=DEFAULT_METRIC_SET,
                how="minmax", top_k=TOP_K, seed=0):
    """가중치를 ±pct 흔들었을 때 상위 top_k 구성이 얼마나 유지되는지."""
    norm = build_norm(rows, metric_set, how)
    base = [c for c, _ in rank_with(rows, weights, norm, metric_set, how)[:top_k]]
    keep = []
    for i in range(n):
        t = [c for c, _ in rank_with(rows, perturb(weights, pct, seed + i),
                                     norm, metric_set, how)[:top_k]]
        keep.append(len(set(base) & set(t)) / top_k)
    return {"base_top": base, "mean_retained": sum(keep) / len(keep),
            "min_retained": min(keep), "pct": pct, "n": n}
