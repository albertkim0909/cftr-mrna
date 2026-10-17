"""
4단계 — 점수화
지표를 후보 집단 안에서 0~1로 정규화한 뒤 가중합한다.
가중치는 잠정값이며, 민감도 분석으로 결론이 이 값에 좌우되지 않음을 확인해야 한다.
"""

# (지표, 방향, 가중치) — 방향 +1은 클수록 좋음, -1은 작을수록 좋음
WEIGHTS = [
    ("acc_kozak",   +1, 0.30),   # 개시부가 열려 있을수록 리보솜 진입이 쉽다
    ("dG_init",     +1, 0.20),   # 개시부 MFE가 높을수록(덜 음수) 장벽이 낮다
    ("cai",         +1, 0.25),   # 코돈이 인간에 적합할수록 신장이 고르다
    ("acc_start30", +1, 0.10),   # 개시 직후 구간의 접근성
    ("acc_mean",    -1, 0.15),   # 전체가 너무 열려 있으면 분해에 취약하다
]


def normalize(values):
    lo, hi = min(values), max(values)
    if hi - lo < 1e-12:
        return [0.5] * len(values)
    return [(v - lo) / (hi - lo) for v in values]


def score_all(rows):
    """rows: [{'id':..., 'metrics': {...}}, ...] — 각 row에 score와 norm을 채워 넣는다.

    is_control=True인 행은 음성 대조군이다. 정규화 기준(최소·최대)은 정상 후보만으로
    잡고 대조군에는 같은 기준을 적용한다. 대조군이 기준 범위를 벗어나면 점수가
    0~1 밖으로 나올 수 있는데, 그 자체가 판별 결과이므로 자르지 않는다.
    순위는 정상 후보에만 매기고, 대조군에는 would_rank(끼워 넣었다면 몇 위인지)를 적는다.
    """
    main = [r for r in rows if not r.get("is_control")]
    ctrl = [r for r in rows if r.get("is_control")]
    ref = main or rows

    for key, direction, _ in WEIGHTS:
        vals = [r["metrics"][key] for r in ref]
        lo, hi = min(vals), max(vals)
        span = hi - lo
        for r in rows:
            v = 0.5 if span < 1e-12 else (r["metrics"][key] - lo) / span
            r.setdefault("norm", {})[key] = 1 - v if direction < 0 else v

    for r in rows:
        r["score"] = sum(w * r["norm"][k] for k, _, w in WEIGHTS)

    main.sort(key=lambda r: -r["score"])
    for rank, r in enumerate(main, 1):
        r["rank"] = rank
    ctrl.sort(key=lambda r: -r["score"])
    for r in ctrl:
        r["rank"] = None
        r["would_rank"] = sum(1 for m in main if m["score"] > r["score"]) + 1
    return main + ctrl
