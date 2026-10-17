# 조절 요소 라이브러리

`utr_library.csv` — 5'UTR, Kozak, 3'UTR 부품 목록. 조립 모듈(3-2)이 이 파일을 읽어 후보를 전수 생성한다.

## 열 설명

| 열 | 내용 |
|---|---|
| `id` | 부품 식별자 |
| `region` | `5UTR` / `KOZAK` / `3UTR` / `3UTR_CONTROL`(음성 대조군, 조립 시 제외) |
| `length` | 염기 수 |
| `sequence` | 서열 (DNA 표기, T 사용) |
| `gc_percent` | GC 함량 (%) |
| `source` | 유래 |
| `accession` | 공개 데이터베이스 등록번호 |
| `trimmed` | 말단 제한효소 인식부위 제거 여부 (`yes`/`no`) |
| `note` | 특이사항 |

## 말단 제한효소 인식부위 처리 규칙 (2-7)

공개 서열에서 떼어낸 UTR에는 클로닝용 제한효소 인식부위가 섞여 들어오는 경우가 있다. 이는 발현 벡터 조립에 쓰인 흔적이며 UTR 본래 기능과 무관하므로, 다음 규칙으로 처리한다.

1. 서열 **양 끝 10염기 이내**에 나타나는 인식부위는 벡터 유래로 보고, 해당 부위와 그 바깥 염기를 함께 제거한 판을 별도 항목으로 등록한다 (`trimmed = yes`).
2. 서열 **내부**에 나타나는 인식부위는 원서열의 일부로 보고 제거하지 않는다. 여섯 염기 조합은 긴 서열에서 우연히 등장할 수 있다.
3. 원본과 제거판을 **모두 후보로 유지**한다. 어느 쪽을 쓸지는 사람이 판단하지 않고, 4-1 단일 변인 비교 실험에서 점수·순위 차이를 측정해 결정한다.

검사 대상 인식부위: XhoI(CTCGAG), NheI(GCTAGC), SalI(GTCGAC), SpeI(ACTAGT), BamHI(GGATCC), EcoRI(GAATTC), HindIII(AAGCTT), NotI(GCGGCCGC)

## 출처 유형 표기

라이브러리 항목의 근거는 두 종류이며, 신뢰도가 다르므로 보고서에서 구분해 기술한다.

- **논문 검증** — `cyba_*`(Ferizi et al. 2016, Sci Rep), `mtrnr1_emcv_3utr`(Kirshina et al. 2023, Biomolecules)
- **특허 기재** — `rabb_*`·`mtrnr1_emcv_3utr`(RU2831165C1), `tpl_5utr_mod`·`aes_mtrnr1_3utr_mod*`(RU2831163C1, Sirius University 2024). 동료 심사를 거치지 않은 출원 문서이므로 성능 수치는 출원인 주장으로 취급한다.

## 현재 적용 결과

- `bnt162b2_5utr` — SpeI(8-13) 포함. 1-13 제거한 `bnt162b2_5utr_trim`(34 nt) 추가.
- `bnt162b2_3utr` — XhoI(4-9), NheI(293-298) 포함. 양 끝 제거한 `bnt162b2_3utr_trim`(283 nt) 추가. 내부 NheI(31)은 규칙 2에 따라 보존.
- `prior_3utr` — 41 nt 중 인식부위 4곳(XhoI 5, SalI 10, SpeI 14, XhoI 36)이 앞쪽에 연속으로 몰린 다중 클로닝 자리 형태. 치료용 3'UTR로서의 기능이 없다고 판단해 **정상 후보군에서 제외하고 음성 대조군으로 전환**(`region = 3UTR_CONTROL`). 조립 모듈(3-2)은 이 항목을 조합 대상에서 빼고, 4-1 비교 실험에서만 따로 투입해 점수 체계가 양·불량을 구별하는지 확인한다. 당시 어느 벡터에서 가져온 서열인지는 기록이 없어 미확인.
- `aes_mtrnr1_3utr_mod` — 양 끝에 XhoI(1-6), NheI(289-294) 포함. 양쪽을 제거한 `aes_mtrnr1_3utr_mod_trim`(282 nt) 추가. 내부 NheI(28)은 규칙 2에 따라 보존.
- `mtrnr1_emcv_3utr` — NheI(147)이 있으나 서열 한가운데이므로 규칙 2에 따라 보존. 제거판 없음.
- `tpl_5utr_mod` — XhoI(175), HindIII(210)이 있으나 모두 내부이므로 규칙 2에 따라 보존. 제거판 없음.
- `cyba_5utr` / `cyba_3utr` / `rabb_5utr` / `rabb_5utr_mod` — 검사 대상 인식부위 없음.
- `mrna1273_3utr` — 미확보. 공개 서열에 UTR 구간 좌표 주석이 없어 분리 실패. 종결코돈 기준 직접 분리를 시도 중이며, 확보 시 이 줄의 `sequence`·`length`·`gc_percent`를 채운다.
