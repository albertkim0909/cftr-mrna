"""
5-1·5-2 — 웹 최소판 백엔드

입력 칸 하나, 실행, 결과 표와 그래프까지만 둔다. 설계 판단은 전부 기존 모듈과
results/delphi/final_1120/README.md의 결정을 따르며, 이 파일은 그것을 호출해
JSON으로 내보내는 얇은 층이다.

실행:
    cd code && uvicorn api:app --host 0.0.0.0 --port 8000
"""
import csv
import threading
import traceback
import uuid
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

import pipeline
import score

ROOT = Path(__file__).resolve().parent.parent
WEBAPP = Path(__file__).resolve().parent / "webapp"
PRECOMPUTED = ROOT / "results" / "3utr_pas" / "ranking_pas1120_112.csv"

METRIC_KEYS = ["acc_kozak", "dG_init", "cai", "acc_start30",
               "acc_utr5", "acc_cds", "acc_utr3"]

NOTICES = [
    "점수 체계는 CFTR을 기준으로 지표와 가중치를 정한 것입니다. 다른 단백질에 적용하려면 "
    "치료 역치와 가중치를 다시 설정해야 합니다.",
    "구간 평균 접근성 지표는 3'UTR의 우열을 변별하지 못하며, 기능 요소(폴리아데닐화 신호) "
    "제약으로 최소 요건만 거를 수 있습니다. 상위 후보 선정은 나머지 지표가 주도합니다.",
]

app = FastAPI(title="CFTR mRNA 역설계")
_jobs = {}
_lock = threading.Lock()


@app.exception_handler(Exception)
def _any_error(request, exc):
    """예상 못 한 오류도 JSON으로 돌려준다.

    기본 처리는 평문 "Internal Server Error"를 보내는데, 프런트가 JSON으로
    읽다가 SyntaxError만 띄우게 되어 진짜 원인이 터미널에만 남는다. 학교
    발표 환경에서는 터미널을 못 볼 수도 있으므로 화면에 원인을 드러낸다.
    """
    # format_exc()는 핸들러 안에서 활성 예외가 없어 비어 나온다. 예외 객체에서
    # 직접 추적을 뽑아야 한다.
    trace = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
    return JSONResponse(
        status_code=500,
        content={"detail": f"{type(exc).__name__}: {exc}", "trace": trace[-1500:]},
    )


class RunRequest(BaseModel):
    protein: str
    require_polya: bool = True
    merge_duplicates: bool = True


# ---------------------------------------------------------------------------
# 결과 정리
# ---------------------------------------------------------------------------
def _contributions(row, weights):
    """1위 후보의 지표별 기여도 — 가중치 × 정규화값을 점수 대비 비율로."""
    out = []
    total = sum(weights[k] * row["norm"][k] for k in METRIC_KEYS)
    for k in METRIC_KEYS:
        part = weights[k] * row["norm"][k]
        out.append({"metric": k, "weight": round(weights[k], 4),
                    "normalized": round(row["norm"][k], 4),
                    "contribution": round(part, 4),
                    "share": round(part / total, 4) if total else 0.0})
    return out


def _pack(result, top_n=20):
    rows = result["rows"]
    weights = result["weights"]
    top = []
    for r in rows[:top_n]:
        top.append({"rank": r["rank"], "id": r["id"], "score": round(r["score"], 4),
                    "cds_variant": r["cds_variant"], "utr5": r["utr5"],
                    "kozak": r["kozak"], "utr3": r["utr3"],
                    "length": r.get("length"),
                    "metrics": {k: round(r["metrics"][k], 4) for k in METRIC_KEYS}})
    return {
        "protein_len": result.get("protein_len"),
        "n_generated": result.get("n_generated"),
        "n_ranked": result.get("n_ranked"),
        "weights": {k: round(v, 4) for k, v in weights.items()},
        "top": top,
        "scores": [round(r["score"], 4) for r in rows],
        "contributions": _contributions(rows[0], weights) if rows else [],
        "notices": NOTICES,
    }


# ---------------------------------------------------------------------------
# 미리 계산해 둔 CFTR 결과
# ---------------------------------------------------------------------------
def _load_precomputed():
    if not PRECOMPUTED.exists():
        raise HTTPException(503, "미리 계산해 둔 결과 파일이 없습니다.")
    rows = []
    with open(PRECOMPUTED, encoding="utf-8") as f:
        for d in csv.DictReader(f):
            rows.append({
                "id": d["id"], "cds_variant": d["cds_variant"], "utr5": d["utr5"],
                "kozak": d["kozak"], "utr3": d["utr3"],
                "length": int(d["length"]), "score": float(d["score"]),
                "metrics": {k: float(d[k]) for k in METRIC_KEYS},
            })
    merged = [r for r in rows if not pipeline.MERGE_PATTERN.match(r["id"])]
    w = pipeline.mean_weights_split()
    scored = score.score_all(merged, weights=w, metric_set=pipeline.METRIC_SET)
    packed = _pack({"protein_len": 1480, "n_generated": len(rows),
                    "n_ranked": len(scored), "weights": w, "rows": scored})
    packed["source"] = ("results/3utr_pas/ranking_pas1120_112.csv — 전수 1,120개를 "
                        "미리 계산해 둔 결과. 동일 서열 병합 후 재점수화")
    return packed


# ---------------------------------------------------------------------------
# 작업 관리
# ---------------------------------------------------------------------------
def _worker(job_id, protein, require_polya, merge_duplicates):
    def progress(phase, done, total, message):
        with _lock:
            j = _jobs[job_id]
            j.update(phase=phase, done=done, total=total, message=message)

    def cancelled():
        with _lock:
            return _jobs[job_id].get("cancel", False)

    try:
        result = pipeline.run(protein, require_polya=require_polya,
                              merge_duplicates=merge_duplicates,
                              progress=progress, cancelled=cancelled)
        with _lock:
            _jobs[job_id].update(state="done", message="완료", result=_pack(result))
    except Exception as e:
        with _lock:
            _jobs[job_id].update(state="error", message=str(e),
                                 detail=traceback.format_exc()[-2000:])


@app.post("/api/run")
def start(req: RunRequest):
    try:
        protein = pipeline.clean_protein(req.protein)
    except ValueError as e:
        raise HTTPException(400, str(e))
    job_id = uuid.uuid4().hex[:12]
    with _lock:
        _jobs[job_id] = {"state": "running", "phase": "optimize", "done": 0,
                         "total": 0, "message": "준비 중", "protein_len": len(protein)}
    t = threading.Thread(target=_worker, daemon=True,
                         args=(job_id, protein, req.require_polya, req.merge_duplicates))
    t.start()
    return {"job_id": job_id, "protein_len": len(protein)}


@app.get("/api/job/{job_id}")
def status(job_id: str):
    with _lock:
        j = _jobs.get(job_id)
        if not j:
            raise HTTPException(404, "없는 작업입니다.")
        return JSONResponse(dict(j))


@app.post("/api/job/{job_id}/cancel")
def cancel(job_id: str):
    with _lock:
        if job_id not in _jobs:
            raise HTTPException(404, "없는 작업입니다.")
        _jobs[job_id]["cancel"] = True
    return {"ok": True}


@app.get("/api/example")
def example():
    return _load_precomputed()


@app.get("/")
def index():
    return FileResponse(WEBAPP / "index.html")


app.mount("/static", StaticFiles(directory=WEBAPP), name="static")
