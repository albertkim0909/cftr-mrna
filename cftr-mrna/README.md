# CFTR 맞춤형 mRNA 치료제 자동 설계 프로그램

2026 교과심화탐구 발표회 / 팀 "가랑이 사이로 빠져버린 청포도"

정상 CFTR 단백질을 암호화하는 치료용 mRNA 후보를 자동으로 설계하고,
구조 지표와 번역 동역학으로 평가해 순위를 매기는 프로그램이다.

---

## 설치

```bash
pip install ViennaRNA dnachisel biopython matplotlib
```

ViennaRNA 2.7.2, dnaChisel 3.2.16에서 동작을 확인하였다.

## 실행

```bash
cd code
python run_skeleton.py        # 후보 3개, 약 15초 — 배관 점검용
python run_all.py --quick     # 축소판, 약 2분
python run_all.py             # 전체 168개, 약 14분
python plot_gc_vs_access.py   # 그림 생성
```

### 웹 (5-1·5-2)

```bash
pip install fastapi uvicorn
cd code
uvicorn api:app --port 8000     # 브라우저에서 http://localhost:8000
```

임의의 목표 단백질 서열을 받아 후보 전수를 생성하고 순위를 낸다. 실행은 백그라운드로 돌아가며
진행률이 표시된다. CFTR은 전수 1,120개를 미리 계산해 둔 결과가 있어 「CFTR 결과 보기」로 즉시
확인할 수 있다.

소요 시간은 단백질 길이에 비례한다. 80 aa 기준 1,120개에 약 7.5분이며, 300~400 aa면 15~25분을
잡아야 한다. 심사 자리에서 돌려 보일 때는 이 점을 감안할 것.

| 파일 | 역할 |
| --- | --- |
| `code/pipeline.py` | 임의 단백질 → 후보 전수 → 지표 → 점수. CDS 출발점만 역번역으로 만들고 이후는 기존 모듈 호출 |
| `code/api.py` | FastAPI 엔드포인트. 실행·진행률·취소·예시 |
| `code/webapp/index.html` | 단일 페이지. 입력 칸, 실행, 상위 후보 표, 그래프 두 개 |

CFTR처럼 참조 mRNA가 있는 경우와 달리, 사용자가 단백질만 넣으면 CDS 출발점이 없다. 인간 코돈
사용 빈도에 따른 가중 추출로 출발점을 만들고 `backtranslated`로 표기한다(최빈 코돈 역번역을 쓰면
출발점이 `use_best_codon` 결과와 같아져 비교 축이 무너진다). 씨앗은 `design.SEED`로 고정한다.

## 폴더 구조

```
data/
  reference/    참조 서열 (수정 금지)
  prior/        선행 탐구 산출물 — 비교 기준선
  elements/     UTR·Kozak 라이브러리
  validation/   검증용 (백신 서열, CFTR2 변이 목록)
code/
  design.py                후보 생성
  metrics.py               지표 계산
  score.py                 점수화
  run_skeleton.py          관통 골격 (1단계)
  run_all.py               전체 실행 (2~4단계)
  plot_gc_vs_access.py     그림
results/
  candidates_scored.csv    후보별 지표와 순위
  fig1_*.png / *.csv       그림과 그 원자료
```

---

## 설계 축

후보는 세 축의 조합으로 생성된다. 현재 168개.

| 축 | 값 | 개수 |
|---|---|---|
| 코돈 최적화 방식 | `use_best_codon`, `match_codon_usage` | 2 |
| GC 상한 | 0.55, 0.58, 0.60 | 3 |
| 대조군 | 야생형(비최적화) | +1 |
| 5′UTR | HBA1, HBB, mRNA-1273, 없음 | 4 |
| Kozak | `GCCGCCACC`, `GCCACC` | 2 |
| 3′UTR | HBA1, HBB, 선행 탐구본 | 3 |

(2×3 + 1) × 4 × 2 × 3 = 168

## 지표

후보 전수에는 두 가지 계산만 돌린다. 후보당 약 4.5초.

- **개시부 창 RNAfold** — AUG 기준 −30 ~ +90 구간의 MFE와 앙상블 자유에너지
- **전체 길이 RNAplfold** — 위치별 비결합 확률 (창 150, 최대 염기쌍 거리 120)

전체 길이 RNAfold는 5 kb에서 약 6분이 걸리므로 최종 후보에만 사용한다.
RNAplfold로 바꾸면 약 3초로, 128배 빠르다.

여기에 CAI와 GC 함량을 더해 총 9개 지표를 산출한다.

## 점수

지표를 후보 집단 안에서 0~1로 정규화한 뒤 가중합한다.
가중치는 **잠정값**이며, 민감도 분석으로 결론이 이 값에 좌우되지 않음을
확인해야 한다. `score.py`의 `WEIGHTS`에서 수정한다.

---

## 검증

- 참조 CDS와 최적화 CDS의 번역 결과가 UniProt P13569과 일치하는지 매 실행 확인
- 508번 잔기가 F인지 확인
- 최적화 후에도 아미노산 서열이 보존되는지 `assert`로 강제

## 자료 출처

`data/*/SOURCE.md` 참조.

- 참조 서열: NCBI NM_000492.4 (6,070 nt, CDS 71..4513)
- 단백질: UniProt P13569 (1,480 aa)
- HBA1: NM_000558.5 / HBB: NM_000518.5
- mRNA-1273 5′UTR: GenBank OR134578.1
- BNT162b2: GenBank OR134577.1

## 한계

- 점수 가중치가 실측 데이터로 학습된 값이 아니다 (MPRA 회귀는 후속 과제)
- ViennaRNA의 2차 구조 예측은 세포 내 환경을 완전히 반영하지 못한다
- 단일 세포를 가정하며 폐 상피로의 전달 효율을 반영하지 않는다
- 실제 발현량 검증에는 세포 실험이 필요하다. 본 시스템은 스크리닝 도구이다.
