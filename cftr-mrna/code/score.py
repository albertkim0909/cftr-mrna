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
    """rows: [{'id':..., 'metrics': {...}}, ...] — 각 row에 score와 norm을 채워 넣는다."""
    for key, direction, _ in WEIGHTS:
        vals = [r["metrics"][key] for r in rows]
        norm = normalize(vals)
        if direction < 0:
            norm = [1 - v for v in norm]
        for r, v in zip(rows, norm):
            r.setdefault("norm", {})[key] = v

    for r in rows:
        r["score"] = sum(w * r["norm"][k] for k, _, w in WEIGHTS)
    rows.sort(key=lambda r: -r["score"])
    for rank, r in enumerate(rows, 1):
        r["rank"] = rank
    return rows
