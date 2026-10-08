# Mito Analyzer 알고리즘 설명 (v1.0)

예시 이미지가 들어간 단계별 그림 설명은 [`Mito_Analyzer_algorithm.pdf`](Mito_Analyzer_algorithm.pdf), 결과 표의 지표 설명은 [`metrics.md`](metrics.md)에 있습니다.
이 문서는 각 단계가 **무엇을, 어떤 값으로, 왜** 하는지를 코드 기준으로 적었습니다. 특히 세포 ROI를 정하는 1단계를 가능한 한 자세히 설명합니다.

## 목차
0. 입력과 픽셀 크기
1. 세포 ROI 설정 (자세히)
2. green+ / green− 세포 분류
3. 미토콘드리아 분할과 MiNA
4. 미토콘드리아 위의 green
5. 출력, 그룹별 출력, 여러 그룹 분석
6. 임계값 수동 조정과 미리보기
7. 기호와 상수 요약

---

## 0. 입력과 픽셀 크기

### 입력
같은 시야의 3채널 TIFF (ImageJ RGB 내보내기 또는 단일 채널 8/16-bit)입니다. 파일 이름의 `…Ch1_Red / …Ch2_Green / …Ch3_Blue`로 세트를 묶고, 폴더는 하위 폴더까지 찾습니다.
- Red = 미토콘드리아, Green = 관심 단백질(POI), Blue = 핵
- 흰색 주석(스케일 바, R·G·B 모두 > 200)은 지우고 계산합니다 (2 px 팽창).

### 픽셀 크기 (픽셀 → µm)
모든 길이·넓이(µm, µm²)와 크기 상수가 픽셀 크기에 의존합니다. `cell_roi.read_pixel_size`가 다음 순서로 읽습니다.
1. OME-XML의 `PhysicalSizeX` (+ `PhysicalSizeXUnit`)
2. ImageJ가 저장한 `unit=` (micron, µm, nm, mm, cm, m, inch) + `XResolution`
3. TIFF `ResolutionUnit` (2 = inch, 3 = cm) + `XResolution`
4. 결과가 0.005–10 µm/px 밖이면(예: 72 dpi 화면 해상도) **없음**으로 봅니다. 단위가 없는 `XResolution`도 없음입니다.

없을 때
- GUI: 이미지 세트 표의 **µm/px** 칸에 빨간 `?`가 표시됩니다. Run all을 누르면 그 세트들의 µm/px를 한 번 묻고(기본값 0.099), 칸을 더블클릭해 세트마다 직접 입력할 수도 있습니다. Options의 *Pixel size*에 값을 넣으면 모든 세트에 그 값을 씁니다.
- 명령줄: `--pixel-size-um`. 없으면 0.099 µm/px로 계산하고 경고합니다.
- 사용한 값과 출처는 각 샘플 `_results.xlsx`의 settings 시트 `pixel_size_um`, `pixel_size_source` (TIFF / entered / assumed)에 남습니다.

### 크기 상수의 배율 k
상수는 기준 데이터(픽셀 0.0990 µm)에서 정했습니다. 다른 배율에서도 같은 물리 크기로 동작하도록 **k = 0.099 / 픽셀 크기**(0.2–5로 제한)를 곱합니다. 아래에서 "px"는 기준 픽셀이며, 실제로는 × k (넓이는 × k²)입니다.

---

## 1. 세포 ROI 설정 (자세히)

### 1-0. 문제와 전략
세포막 염색이 없습니다. 그래서 세포 경계를 직접 보지 못하고, 다음 세 가지 단서로 추정합니다.

| 단서 | 채널 | 쓰임 |
|---|---|---|
| 핵 | blue | 세포마다 하나. 세포의 중심이자 **씨앗(seed)** |
| 미토콘드리아 밀도 | red | 세포질을 채움 → **세포 영역**. 세포와 세포 사이에 생기는 가는 "미토가 없는 선" → **경계** |
| 약한 자가형광 | green | 세포 사이의 어두운 골짜기 → **보조 경계** |

전체 흐름
```
blue ──► 1-1 핵 찾기 (초점 핵 + 흐린 핵) ─────────────┐
red  ──► 1-2 세포 영역(전경) + 미토가 끊어지는 선 ────┤
green ─► 1-2 자가형광 골짜기 (보조 경계) ─────────────┤
                                                      ▼
                    1-3 모든 핵에서 동시에 compact watershed
                                                      ▼
                    1-4 ROI 다듬기 (열림, 구멍, 가장 큰 조각)
                                                      ▼
                    1-5 가장자리 끝 잘라내기 (선택, 기본 켜짐)
                                                      ▼
                    1-6 QC: 가장자리 세포 제외 (edge), 약한 경계 표시
```
핵심 가정: **세포 하나 = 핵 하나**. ROI 개수는 찾은 핵(초점 + 흐린 핵) 개수와 같습니다.

### 1-1. 핵 찾기 (`mito_app/nuclei.py`, 기본 *Texture + merge*)
핵 염색에는 세포질에 퍼진 흐린 blue, 샷 노이즈, 핵 안의 어두운 구멍(염색질·핵소체)이 섞여 있습니다. 진짜 핵은 **염색질 무늬(texture)** 가 있고 세포질 blue는 매끈하다는 차이를 씁니다.

1. **후보**: blue를 σ = 2 px로 부드럽게 → 25 µm 크기의 배경을 빼기(4배 축소 영상에서 grey opening 후 확대) → **0.5 × Otsu** 이상인 곳.
   열림(반지름 4 px), 구멍 채우기, 0.4 × 3000 px² 미만 조각 제거.
   *(수동 조정: Thresholds → Nuclei가 이 0.5 × Otsu 값을 대신합니다.)*
2. **조각 나누기**: 후보 마스크의 거리 지도(σ 2 px)에서 h-maxima(높이 ≈ 0.3 µm) 씨앗으로 watershed → 붙은 핵 후보가 조각으로 나뉩니다.
3. **무늬 검사 (조각마다)**: σ 1.5 px 영상에서, 조각을 10 px 깎은 안쪽의
   `무늬 = √(분산 − 잡음분산) / 중앙값` (잡음은 후보 밖 고역 통과 영상의 MAD로 추정).
   - 기준 밝기 ref = 무늬 ≥ 0.25인 조각들의 중앙값 밝기
   - **남김**: 무늬 ≥ 0.22, 또는 무늬 ≥ 0.17이면서 밝기 ≥ 0.7 × ref
   - 나머지(매끈한 세포질·배경 덩어리)는 버림
4. **조각 합치기**: 남은 이웃 조각을 합칩니다. 단, 접촉선이 **어둡고**(경계 중앙값 < 0.8 × 작은 조각 중앙값) **좁으면**(접촉 길이 < 0.6 × 작은 조각 단축) 따로 둡니다.
   → 휘거나 콩팥 모양인 핵은 하나로 남고, 정말 붙어 있는 두 핵은 둘로 남습니다.
5. **윤곽 다듬기 (핵마다)**: 핵 중심부와 주변 고리의 중간값을 그 핵만의 임계값으로 다시 자르고(넓이가 0.6–1.5배일 때만 채택), 볼록 껍질 안에서 닫힘(10 px)으로 가장자리 홈을 메웁니다. solidity < 0.8이면 버립니다.
6. **초점 / 흐린 핵**: 무늬 검사를 다시 하고, 밝기 ≥ 0.5 × ref이고 넓이 ≥ 3000 px²이면 **초점 맞은 핵**, 아니면 **흐린 핵**(넓이 ≥ max(4000 px², 초점 핵 중앙 넓이의 절반)일 때만).
   흐린 핵도 자기 세포를 가집니다(옵션 *Dim nuclei also get their own cell*, 기본 켜짐). 그 세포는 표의 `seed` = `dim nucleus`입니다.

다른 방식: *Texture, no merge*(4번 생략), *Global Otsu (old)*(부드럽게 한 blue의 Otsu + 흐린 핵은 Otsu 40 % 이상이고 둥근 덩어리).

### 1-2. 세포 영역과 경계 단서 (`cell_roi.segment_morph`, `cell_roi.segment`)
**(a) 세포 영역(전경)**
- 미토 밀도 `md3` = red를 σ = 3 µm로 부드럽게 ÷ 그 99번째 백분위수
- **전경 = md3 > 0.12** (또는 핵 위). 닫힘(2 µm)과 200 µm² 미만 구멍 채우기
- 의미: 미토콘드리아가 어느 정도 있는 곳이 세포입니다. 0.12는 "가장 빽빽한 곳의 12 %".
- *(수동 조정: Thresholds → Cell area가 0.12를 대신합니다. 세포가 너무 작게 잡히면 낮추고, 배경이 세포에 붙으면 높입니다.)*

**(b) 미토 밀도 지도와 끊어지는 선**
- `md` = red를 σ = 1.5 µm로 부드럽게 ÷ 99번째 백분위수 (0–1로 자름)
- `ridge` = md에 Sato 어두운 능선 필터(σ = 1.5, 2.5 µm) ÷ 99.5번째 백분위수 (0–2로 자름)
  → 미토콘드리아 사이를 지나는 가는 어두운 선에서 큽니다. 세포와 세포 사이 경계가 보통 여기에 있습니다.

**(c) green 경계 지도** (`segment`의 landscape)
- green을 자동 보정(배경 0, 세포질 기준 밝기 맞춤) → log1p(최대 20) → σ 4 px
- Sato 어두운 능선(σ 10, 16 px)을 99.5번째 백분위수로 나눈 `vn` + 0.5 × (1 − 정규화한 green 밀도)
  → green 자가형광이 끊어지는 골짜기와 green이 거의 없는 곳에서 큽니다.

### 1-3. 핵에서 동시에 퍼지는 compact watershed
- **비용** = (1 − md) + ridge + green 경계 지도
  - 미토가 빽빽한 곳: 비용이 낮아 쉽게 넘어감
  - 끊어지는 선, green 골짜기: 비용이 높아 경계가 거기에 놓임
- **씨앗** = 모든 핵(초점 + 흐린 핵)의 픽셀, **마스크** = 1-2(a)의 전경
- `skimage.segmentation.watershed(cost, nuclei, mask=전경, compactness = 0.003 / k)`
- **compactness가 하는 일**: 핵에서 멀어질수록 비용에 거리 × 0.003을 더합니다.
  - 경계 단서가 강한 곳(끊어지는 선이 뚜렷함): 선이 이깁니다.
  - 단서가 약한 곳(미토가 고르게 이어짐): 두 핵의 대략 중간에서 만납니다.
  - 그래서 한 세포가 비용이 낮은 길을 따라 길게 뻗어 이웃을 삼키지 않고, 각 세포가 핵을 중심으로 고르게 자랍니다.
- 결과적으로 **모든 ROI에 핵이 정확히 하나** 들어갑니다. 핵을 찾지 못한 세포의 세포질은 이웃 ROI에 붙습니다(핵 검출이 중요한 이유).

### 1-4. ROI 다듬기 (`cell_roi.smooth_labels`)
각 ROI마다
- 열림(반지름 9 px): 가는 돌기와 1–2 px 다리 제거
- 핵 픽셀은 절대 잃지 않음
- 구멍 채우기, 가장 큰 연결 조각만 남김

### 1-5. 가장자리 끝 잘라내기 (`cell_roi.trim_edge_cells`, 옵션 *Keep cells whose tip touches the border*, 기본 켜짐)
핵은 안쪽에 있는데 ROI의 끝부분만 이미지 가장자리(3 px 이내)에 닿는 세포를 살립니다.
1. 대상: ROI가 가장자리에 닿고, 핵이 가장자리에서 **3 µm 이상** 떨어진 세포
2. ROI 안에서 씨앗 두 개로 watershed (같은 비용 지도):
   - 씨앗 1 = 핵 + 핵 둘레 2 µm (가장자리 근처 제외)
   - 씨앗 2 = 가장자리에 닿은 픽셀
3. 핵 쪽 조각을 남기고, 아래 조건을 **모두** 만족하면 가장자리 쪽을 잘라냅니다.
   - 남는 부분 ≥ 원래 ROI의 **55 %**, 세포질 ≥ 0.5 × 핵 넓이
   - 자른 선이 **어두움**: 선 위 미토 밀도 중앙값 ≤ **0.8** × 남는 세포질의 중앙값 (또는 잘린 부분이 ROI의 10 % 미만인 작은 조각)
   - 핵이 자른 선에서 3 µm 이상 떨어짐
   - 남는 부분이 더는 가장자리를 따라 붙어 있지 않음 (접촉 ≤ 3 µm, 선의 양 끝 제외)
   - 남는 부분이 둥근 편: solidity ≥ **0.70**
4. 잘린 넓이는 표의 `edge_trimmed_um2`, summary 이미지에서는 빗금으로 표시됩니다.

실제 데이터에서 잘려 나간 부분은 대부분 핵을 찾지 못한 옆 세포의 미토였습니다. 그래서 이 단계는 너무 크게 잡힌 ROI를 바로잡는 효과도 있습니다.

### 1-6. QC
- **edge**: ROI가 가장자리 3 px 이내에 닿으면 분석에서 제외 (옵션 *Also analyse cells touching the image border*로 포함 가능)
- **weak_border_with**: 이웃 ROI와의 경계가 green 밀도 기준으로 어둡지 않으면(상대 대비 < 0.0) 목록에 남깁니다. *Exclude cells that share a weak border*를 켜면 이 세포들을 binucleate로 제외합니다.
- ROI는 Fiji ROI Manager용 `_RoiSet.zip`(전부)과 `_RoiSet_filtered.zip`(QC 통과), 라벨 TIFF, `_rois.csv`, `_summary.png`로 저장됩니다.

### 1-6b. 수동 검토 (옵션 *Review and edit the cell ROIs before the analysis*, 기본 켜짐)
1-1~1-5로 찾은 ROI를 분석 전에 사용자가 확인·수정합니다 (`mito_app/roi_review.py`, `mito_app/roi_edit.py`).
- 도구: 선택, 영역 더하기, 빼기, 새 세포 그리기, 삭제, 닿아 있는 두 세포 합치기, 되돌리기.
- 세포는 항상 구멍 없는 한 덩어리로 유지됩니다. 편집으로 갈라지면 원래 세포(또는 핵)와 겹치는 조각만 남깁니다.
- 확정하면 라벨을 1..n으로 다시 매기고, 핵 라벨도 같은 번호로 맞춘 뒤, 1-6 QC부터 이후 단계가 이 ROI로 진행됩니다. 손으로 그린 세포는 `seed = manual` (핵 없음).
- 확정한 ROI는 `<sample>_roi_review.npz`로 저장되어 다음 실행에서 다시 불러올 수 있습니다.

### 1-7. ROI가 이상할 때 확인할 것
| 증상 | 가능한 원인 | 조정 |
|---|---|---|
| 세포 두 개가 한 ROI | 한쪽 핵을 못 찾음 | Thresholds → Nuclei를 낮춰 미리보기, 또는 *Dim nuclei* 켜기 |
| 배경이 ROI에 붙음 | 세포 영역 기준이 낮음 | Thresholds → Cell area를 높임 (예: 0.15–0.2) |
| ROI가 세포보다 작음 | 세포 영역 기준이 높음, 미토가 희미함 | Cell area를 낮춤 (예: 0.08) |
| 가짜 세포 (세포질 덩어리에 ROI) | 매끈한 blue 덩어리가 핵으로 잡힘 | Nuclei를 높이거나 *Texture + merge* 방식 확인 |
| 가장자리 세포가 많이 빠짐 | 끝 잘라내기 조건 불만족 | summary의 빗금 확인, 필요하면 edge 세포 포함 옵션 |

---

## 2. green+ / green− 세포 분류 (`mito_app/green_cells.py`)
- 밝은 green 마스크: 보정한 green을 σ 1.5 px로 부드럽게 → log → Otsu 두 번(세포질 vs 신호, 약한 vs 강한 신호). 단, 세포질 기준 밝기의 3배 미만은 밝다고 보지 않습니다.
  *(수동 조정: Thresholds → Green+)*
- 세포 점수 = 세포질(세포 − 핵) 중 마스크가 덮는 비율. **≥ 2 %**(옵션)이면 green+.
- ROI 모양은 바꾸지 않습니다. 이후 모든 표·상관은 green+ / green−를 따로도 계산합니다.

## 3. 미토콘드리아 분할과 MiNA (`mito_app/mito_objects.py`, `mito_app/mina.py`)
기본 *Split objects*
1. 전처리: rolling ball(반지름 1.5 µm) 배경 제거 → σ 0.7 px
2. 전경 (세포마다): 국소 평균(1.1 µm 창)보다 밝고, **세포 Otsu × 0.5** 이상. 6 px² 이하 구멍 채움, 8 px² 미만 조각 제거.
   *(수동 조정: Thresholds → Mitochondria. 켜면 이 두 조건 대신 전처리한 red > 값 하나로 자릅니다. 아래 조각 나누기는 그대로입니다.)*
3. 밝기 봉우리(h-maxima, h는 이미지 밝기 범위에 비례)에서 watershed로 잘게 나눔
4. 다시 합치기: 두 조각은 접촉점이 어둡고(접촉 최대 밝기 / 어두운 봉우리 < 0.75) **동시에** 좁을 때(접촉 길이 < 얇은 쪽 폭)만 따로 둠
5. 서로 닿는 조각 사이에 1 px 틈을 둬 MiNA가 따로 셉니다.

*MiNA classic*: 세포마다 red Otsu 하나 (수동이면 전처리한 red > 값).
그다음 골격화(Lee) → 가지·분기점·끝점·네트워크·도넛 수와 길이 (Fiji MiNA와 같은 정의, 모집단 SD). 조각화 지표는 [`metrics.md`](metrics.md) 참고.

## 4. 미토콘드리아 위의 green (`mito_app/green.py`)
- green에서 배경(세포 사이 어두운 틈 밝기)을 뺌
- 세포 단위: 미토 마스크 위/밖 평균, 농축비, 미토 위 비율, Pearson, Manders
- green puncta: σ 1 px → white top-hat(반지름 8 px) → 이미지 하나에 임계값 하나 = max(Otsu, 중앙값 + 6 MAD) × sensitivity (분석 세포 안에서 계산)
  *(수동 조정: Thresholds → Green puncta, top-hat 단위)*
- 미토 조각 단위: 형태 + 그 위의 green, 겹치는 puncta

## 5. 출력, 그룹별 출력, 여러 그룹 분석
### 샘플마다
`<출력>/<dataset>/<preset>-<dataset>/<sample>/`에 `_results.xlsx`(green_pos / green_neg / all, settings), `_per_cell / _per_mito / _green_puncta / _correlations` CSV(각 green+/− 버전), ROI zip, 오버레이 이미지가 저장됩니다.

### 그룹
그룹은 이미지 세트 표의 **Group** 열(기본 = dataset)입니다. `<출력>/groups.csv`에 저장되고, Analysis 탭의 *Groups…* 에서 바꿀 수 있습니다. 명령줄: `--group`.

### 그룹별 출력 (배치가 끝나면 자동, 또는 Analysis → *Export by group…*)
```
<출력>/groups/
  groups_overview.csv                  그룹별 이미지·세포·green+ 세포·미토 조각 수
  group_comparison.xlsx                지표 × 그룹 비교표 (아래)
  group_comparison_cells.csv           (세포 단위, 전체 세포)
  group_comparison_cells_per_image.csv (이미지 단위)
  <그룹>/
    <그룹>_cells.csv, _mito.csv                    그 그룹의 모든 세포 / 미토 조각
    <그룹>_cells_green_pos.csv, _cells_green_neg.csv (미토도 같음)
    <그룹>_per_image.csv                           이미지당 평균 (이미지 하나 = 한 줄)
    <그룹>_correlations.csv                        모든 지표 쌍의 회귀 / Spearman
    <그룹>_results.xlsx                            위 표 전부
```
`group_comparison.xlsx`는 **긴 형식**(한 줄 = 지표 × 그룹)입니다. 그룹이 많아도 열이 늘지 않습니다. 열 구성:
- n, mean, SD, SEM, 95 % CI, median, Q1, Q3
- 전체 그룹 검정: 2그룹이면 Welch t, Mann-Whitney U / 3그룹 이상이면 ANOVA, Kruskal-Wallis
- 대조군을 고르면(`--control` 또는 내보내기 창): 각 그룹 vs 대조군 Welch, Mann-Whitney와 **Holm 보정** p
- 시트: 전체 / green+ / green−, 각각 세포 / 이미지 / 미토 조각 단위

### 여러 그룹 분석 (Analysis 탭)
- **Groups: mean / median**
  - *Control*: 각 그룹을 대조군과 비교 (Mann-Whitney, Holm 보정; 그림에 *, **, ***)
  - *Groups shown*: 보여줄 그룹 선택
  - 그룹이 6개를 넘으면 x축 이름을 기울이고 점을 작게 그림
  - *Table: this metric only*로 지금 지표만 표에 표시
  - *Unit = images*: 이미지 하나를 한 번만 셈
- **Groups: correlation**
  - *A vs B*: 두 그룹 히트맵 + 차이(Fisher z 검정). ◀ ▶로 B를 차례로 바꿈
  - *All groups*: 그룹마다 작은 히트맵 + 모든 그룹에서 상관이 같은지 보는 **Cochran's Q 검정** 지도
  - 그룹별 상관 행렬은 한 번 계산하면 다시 써서 빠르게 바뀝니다.
- 많은 칸을 동시에 검정하면 약 5 %는 우연히 p < 0.05가 됩니다. 화면에 "몇 칸 중 몇 칸, 우연 기대 수"를 함께 표시합니다.

## 6. 임계값 수동 조정과 미리보기
기본은 이미지마다 자동입니다. Options → **Thresholds**에서 *Manual*을 켜면 그 값 하나를 배치의 **모든 이미지에 일괄 적용**합니다.

| 임계값 | 적용되는 영상 | 자동 값 | 수동 값의 단위 |
|---|---|---|---|
| Nuclei | 배경 제거한 blue (Otsu 방식: 부드럽게 한 blue) | 0.5 × Otsu (Otsu 방식: Otsu) | blue 밝기 |
| Cell area | 미토 밀도 (red σ 3 µm ÷ p99) | 0.12 | 0–1 |
| Mitochondria | 전처리한 red (rolling ball 1.5 µm, σ 0.7 px) | 세포마다 국소 평균 AND 0.5 × Otsu | red 밝기 |
| Green+ (bright green) | 부드럽게 한 보정 green | Otsu 두 번, ≥ 세포질 3배 | 보정 green 단위 |
| Green puncta | green top-hat | max(Otsu, 중앙값 + 6 MAD) × sensitivity | top-hat 단위 |

**미리보기** (*Preview / adjust thresholds…*)
- 이미지 세트와 임계값을 고르고 슬라이더(또는 숫자)를 움직이면 마스크가 바로 다시 그려집니다. 그 이미지의 자동 값이 옆에 표시됩니다. 확대·이동 도구가 있고, 마스크 / 윤곽선 / 원본 보기를 고를 수 있습니다.
- **Apply to all images**: Options의 수동 값으로 설정 (모든 이미지)
- **This image only**: 그 이미지 세트에만 적용 (표의 *Thresholds* 칸에 표시). 개별 값이 일괄 값보다 우선합니다.
- 미리보기는 분석과 같은 함수(`mito_app/thresholds.py`)로 계산하므로, 미리보기에서 맞춘 값이 분석에서도 같은 마스크를 만듭니다. 단, Mitochondria 자동 값은 세포마다 다르므로 미리보기의 "auto"는 참고값이고, Green puncta 자동 값은 분석 세포 안에서 계산되므로 조금 다를 수 있습니다.
- 사용한 값은 settings 시트에 자동 값과 함께 남습니다 (`nuclei_threshold`, `nuclei_threshold_auto`, `cell_fg_level`, `mito_threshold`, `puncta_threshold`, `green_bright_threshold`, `param_thr_*`).
- 명령줄: `--thr-nuclei`, `--cell-fg-level`, `--thr-mito`, `--thr-green-bright`, `--thr-puncta` (0 = 자동).

## 7. 기호와 상수 요약
| 기호 / 상수 | 값 | 뜻 |
|---|---|---|
| k | 0.099 / 픽셀 크기 | 크기 상수 배율 |
| 세포 영역 기준 | 0.12 × p99 | 미토 밀도(σ 3 µm) |
| compactness | 0.003 / k | 핵에서의 거리 벌점 |
| 끝 잘라내기 | 핵 ≥ 3 µm, 남는 ≥ 55 %, 선 ≤ 0.8 × 세포질, solidity ≥ 0.70 | 가장자리 세포 살리기 조건 |
| 핵 무늬 | ≥ 0.22, 또는 ≥ 0.17 & 밝기 ≥ 0.7 × ref | 진짜 핵 판정 |
| 핵 합치기 | 경계 < 0.8 × 내부 AND 접촉 < 0.6 × 단축이면 따로 | 붙은 핵 구분 |
| green+ | 밝은 green ≥ 세포질의 2 % | 발현 세포 |
| 미토 분할 | 0.5 × 세포 Otsu, saddle 0.75, neck 1.0 | 조각 나누기 |
| puncta | max(Otsu, 중앙값 + 6 MAD) | 밝은 점 |
