# Mito Analyzer 사용 설명서 (v1.1)

English version: [manual_en.md](manual_en.md)

이 설명서는 Mito Analyzer를 처음 쓰는 사람이 설치부터 결과 해석까지 혼자 따라 할 수 있도록 쓴 것입니다.
각 단계의 계산 방법은 [algorithm.md](algorithm.md)에, 결과 표의 지표 뜻은 [metrics.md](metrics.md)에 자세히 있습니다.
설명서의 그림은 합성(가짜) 예시 이미지로 만든 화면입니다.

## 목차
1. 이 프로그램이 하는 일
2. 설치
3. 이미지 준비
4. 처음 분석해 보기 (따라 하기)
5. 화면 구성: Run 탭
6. 옵션 하나하나
6-1. 세포 ROI 검토와 수정
7. 임계값 수동 조정과 미리보기
8. 픽셀 크기 (µm/px)
9. 그룹 만들기
10. 결과 보기: 샘플별 탭
11. Analysis 탭: 여러 샘플 통합 분석
12. 출력 파일
13. 결과를 해석하는 법
14. 명령줄(터미널)로 실행
15. 문제 해결 (FAQ)

---

## 1. 이 프로그램이 하는 일
세 채널(빨강·초록·파랑) 형광 이미지에서 **세포 하나하나**를 찾아, 각 세포의
- **미토콘드리아 형태** (얼마나 길고 연결되어 있는지, 얼마나 잘게 쪼개졌는지; Fiji MiNA와 같은 방식)
- **관심 단백질(green)이 미토콘드리아에 얼마나 있는지** (미토 위 green 양, 비율, 밝은 점(puncta))

를 측정하고, 둘 사이의 관계(상관·회귀)와 **실험군 사이의 차이**를 보여 줍니다.

| 채널 | 염색 예 | 프로그램에서의 역할 |
|---|---|---|
| Red (빨강) | MitoTracker 등 미토콘드리아 | 미토 형태 측정, 세포 경계 찾기 |
| Green (초록) | 관심 단백질 (예: mNeonGreen 융합) | 미토 위 단백질 양, green+/− 세포 구분 |
| Blue (파랑) | DAPI / Hoechst 핵 | 세포 하나 = 핵 하나, 세포의 중심 |

## 2. 설치
세 가지 방법 중 하나를 고르세요. **분석만 할 사람은 2-1**, **코드를 고치며 쓸 사람은 2-3**이 편합니다.

### 2-1. DMG 파일로 설치 (가장 쉬움, Python 필요 없음)
1. `Mito Analyzer-1.1.1-apple-silicon.dmg` 파일을 받아 더블클릭합니다.
2. 열린 창에서 **Mito Analyzer** 아이콘을 **Applications** 폴더 아이콘 위로 끌어다 놓습니다.
3. 처음 열 때 "확인되지 않은 개발자" 경고가 뜨면 (Apple 서명이 없는 앱이라 그렇습니다)
   - 응용 프로그램 폴더에서 앱을 **우클릭 → 열기 → 열기**, 또는
   - **시스템 설정 → 개인정보 보호 및 보안** 맨 아래 **"그래도 열기"** 를 누른 뒤 다시 엽니다.
   한 번만 하면 다음부터는 그냥 열립니다.
- Apple Silicon(M1 이후) Mac 전용입니다. Intel Mac은 2-2나 2-3을 쓰세요.

### 2-2. 터미널 한 줄 설치 (자동 빌드)
터미널(응용 프로그램 → 유틸리티 → 터미널)을 열고 아래 한 줄을 붙여 넣고 Enter:
```bash
curl -fsSL https://raw.githubusercontent.com/windowrainYOON/mito-tracking/main/install.sh | bash
```
코드를 받고, Python이 없으면 관리자 권한 없이 사용자 폴더에 설치하고, 앱을 만들어 Dock에 고정합니다 (처음 5분 안팎). 같은 줄을 다시 실행하면 최신 버전으로 업데이트됩니다.
- 설치 위치: `~/Applications/MitoAnalyzer` (바꾸려면 `MITO_DIR=경로` 를 앞에 붙임)
- 아직 main 브랜치에 반영되기 전이라면 개발 버전을 받으세요: `curl … | MITO_REF=develop bash`

### 2-3. 직접 빌드 (코드를 고치며 쓸 때)
```bash
git clone https://github.com/windowrainYOON/mito-tracking.git
cd mito-tracking/application
./build_mac.sh --add-to-dock
```
- Python 3.11 이상이 필요합니다 (3.12 이상 권장). 없으면 스크립트가 설치 방법을 알려 줍니다.
- 코드를 고친 뒤 `./build_mac.sh` 를 다시 실행하면 같은 자리에 새로 빌드되어 Dock 아이콘이 그대로 새 버전을 엽니다.
- 옵션: `--clean` (환경을 처음부터 다시), `--locked` (검증된 패키지 버전으로), `--dmg` (다른 사람에게 줄 DMG도 만들기)
- Windows / Linux: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt` 후 `.venv/bin/python run_app.py`

## 3. 이미지 준비
- **같은 시야**의 세 채널을 **각각 TIFF 파일**로 저장합니다 (ImageJ/Fiji의 RGB 내보내기 또는 단일 채널 8/16-bit).
- **파일 이름 규칙**: 채널을 뺀 앞부분이 같고, 끝이 `Red` / `Green` / `Blue` 로 끝나야 세 파일이 한 세트로 묶입니다. 채널 번호(`Ch1` 등)는 있어도 되고 없어도 됩니다.
  ```
  Processed_New01_Ch1_Red.tif
  Processed_New01_Ch2_Green.tif
  Processed_New01_Ch3_Blue.tif
  ```
- **폴더 구성 (권장)**: 실험 조건마다 폴더 하나. 폴더를 통째로 넣으면 하위 폴더까지 모두 찾고, 폴더 이름이 *Dataset* 과 기본 *Group* 이 됩니다.
  ```
  내 실험/
    Control/        ← 이미지 세트 여러 개
    PRD 24h/
      Processed/New-01_results/…_Red.tif …   (하위 폴더도 괜찮음)
  ```
- **스케일 바**(흰색 글자·막대)는 자동으로 지우고 계산합니다.
- **픽셀 크기**: Fiji에서 Image → Properties로 µm 단위가 들어간 TIFF면 자동으로 읽습니다. 없으면 8절을 보세요.
- 한 이미지에 세포가 5–20개 정도, 세포 전체가 화면 안에 들어오는 시야가 좋습니다. 가장자리에 걸린 세포는 분석에서 빠질 수 있습니다.

### 3-1. Zeiss CZI 파일
ZEN에서 저장한 `.czi` 파일은 한 시야의 모든 채널을 담고 있어서, **CZI 파일 하나가 이미지 세트 하나**입니다. 파일 이름 규칙은 필요 없습니다.
- TIFF와 똑같이 넣으면 됩니다 (Add files…, Add folder…, 또는 표에 끌어다 놓기). 한 폴더에 TIFF와 CZI가 섞여 있어도 됩니다.
- **채널 역할**은 파일에 저장된 채널 이름으로 자동 추정합니다. DNA 염색(DAPI, Hoechst 등)은 핵, 방출 파장이 가장 긴 채널은 미토콘드리아, 나머지는 관심 단백질(POI)로 잡습니다. 추정 결과는 표의 Red / Green / Blue 칸에 `파일.czi ch0 AF568-T1` 처럼 보입니다. **한 번은 꼭 확인하세요.**
- 여러 파일의 역할을 한꺼번에 바꾸려면: 줄을 선택하고 (아무것도 선택하지 않으면 모든 CZI 줄) **Channel roles…** 를 누르거나 CZI 줄의 Red / Green / Blue 칸을 더블클릭한 뒤, 채널마다 *Mito*, *Protein (POI)*, *Nucleus*, *Not used* 중 하나를 고릅니다. 선택은 채널 번호 기준으로 고른 모든 파일에 적용되고, 나중에 넣는 CZI 중 채널 이름이 같은 파일에도 같은 역할이 자동으로 들어갑니다.

  ![채널 역할](manual_img/12_czi_roles.png)
- 픽셀 크기는 CZI 메타데이터에서 읽습니다 (µm/px 칸 툴팁과 settings 시트에 *CZI* 로 표시).
- z-stack은 최대 강도 투영(MIP)으로 분석하고, scene이나 시간점이 여러 개면 첫 번째를 씁니다. 타일(mosaic)은 이어 붙입니다. 8-bit보다 깊은 이미지는 비트 수마다 정해진 하나의 비율(예: 12-bit ÷ 16)로 0–255에 맞추고, 이미지마다 따로 늘리지 않아서 한 배치 안의 밝기를 그대로 비교할 수 있습니다.
- 원본 CZI는 같은 시야의 *Processed* TIFF보다 보통 어둡습니다. 자동 임계값은 이미지마다 맞춰지지만, processed TIFF에서 정한 수동 임계값을 그대로 쓰지 말고 미리보기로 한 번 확인하세요.

## 4. 처음 분석해 보기 (따라 하기)
![시작 화면](manual_img/01_start.png)

1. **앱을 엽니다.** 위쪽에 **Run** 과 **Analysis (all samples)** 두 탭이 있습니다. 처음에는 Run 탭입니다.
2. **이미지를 넣습니다.** **Add folder…** 를 눌러 실험 폴더(예: `Control`)를 고르거나, Finder에서 폴더를 표로 끌어다 놓습니다. 여러 조건이면 폴더마다 반복합니다.
   ![이미지 세트 표](manual_img/02_image_sets.png)
   - 표의 한 줄 = 이미지 세트 하나(빨강·초록·파랑 세 파일).
   - **µm/px** 칸에 빨간 `?` 가 보이면 그 파일에는 픽셀 크기 정보가 없다는 뜻입니다 (8절). Run all 때 물어봅니다.
3. **그룹을 확인합니다.** *Group* 열이 비교할 실험군입니다 (기본 = 폴더 이름). 바꾸려면 칸을 더블클릭해 고치거나, 여러 줄을 선택하고 **Set group…** 을 누릅니다.
4. **출력 폴더를 정합니다.** *Output folder* 에 결과를 저장할 폴더를 고릅니다. 비워 두면 입력 폴더들의 공통 상위 폴더 안에 `dataset` 폴더를 만듭니다.
   - 원본 데이터 폴더와 다른 곳을 고르는 것을 권장합니다.
5. **옵션은 처음엔 그대로 둡니다.** 기본값이 대부분의 이미지에 맞도록 정해져 있습니다.
6. **Run all** 을 누릅니다. 진행 상황이 아래 *Log* 와 표의 *Status* 칸에 나옵니다. 중간에 멈추려면 같은 버튼(**Stop**)을 누르면 지금 이미지까지만 끝내고 멈춥니다.
   - 먼저 모든 이미지의 **세포 ROI만** 찾은 뒤(이미지당 10–15초) **ROI 검토 창**이 열립니다. 세포를 확인하고 필요하면 고친 다음 **Confirm ROIs and continue analysis** 를 누르면 나머지 분석(이미지당 20초 안팎)이 이어집니다. 자세한 사용법은 6-1절.
   - 검토 없이 끝까지 자동으로 돌리려면 Options의 *Review and edit the cell ROIs before the analysis* 를 끕니다.
7. **결과를 봅니다.** 끝난 줄을 클릭하면 오른쪽 탭에 그 이미지의 결과가 나옵니다 (10절).
   ![분석이 끝난 화면](manual_img/05_after_run.png)
8. **먼저 Cell ROIs 탭에서 세포가 제대로 잡혔는지 확인합니다.** 빨간 굵은 선 = 분석에 쓰인 세포, 회색 점선 = 가장자리에 닿아 빠진 세포.
   ![Cell ROIs](manual_img/06_cell_rois.png)
   잘못 잡힌 세포가 많으면 7절(임계값)과 15절(문제 해결)을 보세요.
9. **실험군을 비교합니다.** 위쪽 **Analysis (all samples)** 탭으로 갑니다. 방금 돌린 결과 폴더가 자동으로 불러와져 있습니다 (11절).

## 5. 화면 구성: Run 탭
| 영역 | 내용 |
|---|---|
| **Image sets** 표 | Group · Dataset · Preset · Sample · µm/px · Red · Green · Blue · Thresholds · Status |
| 표 아래 버튼 | **Add files…** (파일 여러 개 선택), **Add folder…** (폴더, 하위 폴더까지), **Set group…** (선택한 줄의 그룹 정하기), **Channel roles…** (CZI 채널 역할, 3-1), **Remove** (선택한 줄 빼기), **Clear** (모두 빼기) |
| **Output folder** | 결과 저장 위치. 저장 구조: `<출력>/<dataset>/<preset>-<dataset>/<sample>/` |
| **Options** | 분석 설정 (6절). 스크롤해서 아래까지 보세요 |
| **Run all / Stop** | 표의 모든 세트를 차례로 분석 / 지금 세트까지 끝내고 멈춤 |
| **Open output folder** | Finder에서 결과 폴더 열기 |
| **Log** | 진행 상황, 경고 |
| 오른쪽 탭 | 선택한 이미지의 결과 (10절) |

표에서 바꿀 수 있는 칸
- **Group**: 비교할 실험군 이름 (9절)
- **Dataset**: 출력 폴더 이름의 일부 (기본 = 넣은 폴더 이름)
- **Preset**: 출력 폴더 구분용 이름 (기본 = 미토 분할 방식 `split`/`otsu`). 같은 데이터를 다른 설정으로 돌릴 때 이름을 바꾸면 결과가 섞이지 않습니다.
- **Sample**: 이미지 이름 (기본 = 파일 이름에서 채널 부분을 뺀 것)
- **µm/px**: 픽셀 크기. 더블클릭해서 직접 입력 (8절)
- **Red/Green/Blue**: 더블클릭하면 그 채널 파일을 다른 파일로 바꿀 수 있습니다 (CZI 줄은 채널 역할 바꾸기, 3-1). Red = 미토콘드리아, Green = 관심 단백질, Blue = 핵.
- **Thresholds**: 이 이미지에만 쓰는 수동 임계값 (7절). 미리보기 창에서 정합니다.

## 6. 옵션 하나하나
![Options](manual_img/03_options.png)

| 옵션 | 기본값 | 뜻 | 언제 바꾸나 |
|---|---|---|---|
| Review and edit the cell ROIs before the analysis | 켜짐 | 세포 ROI를 찾은 뒤 멈추고 검토 창을 엶 (6-1절) | 수정 없이 자동으로 끝까지 돌릴 때 끔 |
| Mito segmentation | Split objects | 미토콘드리아를 나누는 방식. *Split objects*: 국소 임계값 + 붙은 미토를 어둡고 좁은 접촉점에서 나눔 (조각화 지표 포함). *MiNA classic*: 세포마다 Otsu 하나 (Fiji MiNA와 같은 방식) | Fiji MiNA 결과와 직접 비교할 때만 MiNA classic |
| Exclude cells that share a weak border | 꺼짐 | 경계가 뚜렷하지 않은 이웃 세포 쌍(핵 두 개처럼 보이는 세포)을 분석에서 뺌 | 세포가 빽빽하게 붙어 경계가 불확실한 데이터에서 보수적으로 볼 때 |
| Keep cells whose tip touches the border | 켜짐 | 핵은 안쪽인데 끝부분만 가장자리에 닿은 세포를, 미토가 끊기는 선에서 잘라 분석에 포함 (잘린 부분은 빗금으로 표시) | 잘린 세포를 아예 빼고 싶으면 끔 |
| Also analyse cells touching the image border | 꺼짐 | 가장자리에 닿은 세포도 분석 | 보통은 켜지 마세요 (잘린 세포가 섞임) |
| Minimum cell area | 250 µm² | 이보다 작은 세포 ROI는 버림 | 세포가 아주 작은 세포주면 낮춤 |
| Weak-border threshold | 0.000 | "약한 경계" 판정 기준 (위의 Exclude 옵션과 함께) | 거의 바꿀 일 없음 |
| Pixel size (0 = auto) | 0 | 0이면 파일에서 읽음. 값을 넣으면 **모든 이미지**에 그 값 사용 | 모든 파일의 픽셀 정보가 틀렸거나 없을 때 |
| Green puncta threshold | 1.00 × | 자동 puncta 임계값에 곱하는 값. 1보다 크면 더 밝은 점만, 작으면 더 많이 | puncta가 너무 많이/적게 잡힐 때 (또는 7절 수동 값) |
| Minimum mito object | 0.05 µm² | 이보다 작은 미토 조각은 미토 조각 표에서 뺌 | 잡음 점이 많을 때 높임 |
| Nucleus detection | Texture + merge | 핵 찾는 방식. 염색질 무늬로 진짜 핵만 고르고, 휜 핵은 하나로, 붙은 핵은 둘로. *Texture, no merge*, *Global Otsu (old)* 도 선택 가능 | 기본값 권장. 핵 염색이 아주 깨끗하면 Otsu도 됨 |
| Dim (out-of-focus) nuclei also get their own cell | 켜짐 | 초점이 흐린 핵도 자기 세포를 가짐 (옆 세포에 미토가 넘어가지 않게) | 거의 항상 켬 |
| Green+ cell: bright green ≥ | 2.0 % | 세포질 중 밝은 green이 이 비율 이상이면 green+ (발현 세포) | 발현이 약하면 낮춤, 배경이 밝으면 높임 |
| Thresholds | 모두 auto | 5가지 임계값을 자동 또는 수동으로 (7절) | 자동 결과가 맞지 않을 때 |

## 6-1. 세포 ROI 검토와 수정
Options의 **Review and edit the cell ROIs before the analysis** (기본: 켜짐)가 켜져 있으면, Run all은 세포 ROI를 찾은 뒤 잠깐 멈추고 아래 창을 엽니다. 여기서 확인·수정한 ROI로 나머지 분석이 진행됩니다.

![ROI 검토 창](manual_img/12_roi_review.png)

화면
- 왼쪽: 이미지 세트 목록. `(분석 세포 / 전체 세포)`, ✓ = 열어 봄, ✎ = 수정함
- 가운데: 미토(빨강)·green(초록)·핵(파랑) 합성 이미지 위에 세포 경계(흰색), 핵 경계(하늘색), 번호
  - **회색** 경계·번호 = 이미지 가장자리에 닿아 분석에서 빠질 세포, `*` = 손으로 그린 세포, **노란** 굵은 선 = 선택한 세포
  - 위 체크박스로 채널, 핵 경계, 번호를 끄고 켤 수 있습니다
- 확대/이동: 마우스 휠, 또는 돋보기·손 아이콘 (집 아이콘 = 전체 보기)

도구 (괄호 = 단축키)
| 도구 | 사용법 |
|---|---|
| **Select (V)** | 세포를 클릭해 선택. 아래에 넓이와 씨앗(핵 / 흐린 핵 / 손으로 그림)이 표시됩니다 |
| **Add (A)** | 선택한 세포에 더할 영역을 마우스로 둘러 그립니다. 다른 세포의 영역이면 그 세포에서 가져옵니다. 세포와 이어져 있어야 합니다 |
| **Subtract (S)** | 선택한 세포에서 뺄 영역을 둘러 그립니다 (선택이 없으면 그 영역에 걸친 모든 세포에서). 세포 경계에 걸치게 그리세요 (세포 안에 구멍은 만들 수 없습니다) |
| **New cell (N)** | 놓친 세포의 윤곽을 그립니다 (다른 세포가 없는 곳만) |
| **Delete cell (Delete)** | 선택한 세포 삭제 (가짜 세포, 분석에서 뺄 세포) |
| **Merge… (M)** | 선택한 세포를 누른 뒤 Merge… → 합칠 세포 클릭. 하나의 세포가 둘로 나뉘었을 때. 서로 닿아 있어야 합니다 |
| **Undo / Redo** | ⌘Z / ⇧⌘Z |
| **Reset image** | 이 이미지를 처음 검출된 ROI로 되돌림 |
| **Esc** | 선택 해제, Select 도구로 |

자주 쓰는 수정
- **세포 두 개가 한 ROI** → Subtract로 한쪽을 잘라내고, 잘라낸 자리에 New cell로 새 세포를 그립니다.
- **한 세포가 둘로 나뉨** → 한쪽 선택 → Merge… → 다른 쪽 클릭.
- **배경이 ROI에 붙음** → 그 세포를 선택하고 Subtract로 배경 부분을 둘러 그립니다.
- **가짜 세포** → 선택 → Delete.
- **세포가 빠짐** → New cell로 윤곽을 그립니다 (번호에 `*`, 표의 *seed* 열 = `manual`).

마치기
- **◀ Previous / Next ▶** 로 모든 이미지를 확인한 뒤 **Confirm ROIs and continue analysis**. 열어 보지 않은 이미지가 있으면 그대로 진행할지 묻습니다.
- **Cancel** 은 분석을 시작하지 않고 멈춥니다 (아무것도 저장되지 않음).
- 확정한 ROI는 각 샘플 폴더에 `<sample>_roi_review.npz` 로 저장되고, settings 시트에 `rois_reviewed = True` 로 기록됩니다. 같은 이미지를 같은 출력 폴더로 다시 돌리면 **"이전에 검토한 ROI에서 시작할까요?"** 라고 묻습니다 (Yes = 이전 수정 그대로 불러옴, No = 새로 검출).
- 세포를 지우거나 그린 경우에도 나머지 계산(미토, green, 통계)은 확정한 ROI 기준으로 똑같이 이루어집니다.

## 7. 임계값 수동 조정과 미리보기
자동 임계값은 이미지마다 계산됩니다. 결과가 이상하면 직접 정할 수 있습니다.

| 임계값 | 정하는 것 | 너무 낮으면 | 너무 높으면 |
|---|---|---|---|
| Nuclei | 핵 후보 | 세포질 덩어리도 핵 → 가짜 세포 | 흐린 핵을 놓침 → 두 세포가 한 ROI |
| Cell area | 세포 영역 (미토 밀도 기준, 0–1) | 배경까지 세포에 포함 | 세포가 작게 잡힘 |
| Mitochondria | 미토 마스크 | 배경이 미토로 잡힘, 조각이 붙음 | 흐린 미토를 놓침, 과하게 쪼개짐 |
| Green+ (bright green) | 발현 세포 판정용 밝은 green | green− 세포가 green+ 로 | 약한 발현 세포를 놓침 |
| Green puncta | 밝은 green 점 | 잡음이 점으로 | 실제 점을 놓침 |

### 사용 방법
1. Options → Thresholds 아래 **Preview / adjust thresholds…** 를 누릅니다.
   ![임계값 미리보기](manual_img/04_threshold_preview.png)
2. 위에서 **Image set**(이미지)과 **Threshold**(종류)를 고릅니다.
3. **Manual** 을 켜고 슬라이더나 숫자 칸으로 값을 바꿉니다. 색칠된 부분 = 그 값 이상인 픽셀이 바로 다시 그려집니다. 옆에 그 이미지의 **자동 값**이 나옵니다.
   - 돋보기·손 아이콘으로 확대/이동, *Show* 에서 마스크/윤곽선/원본 보기.
4. 값을 정했으면
   - **Apply to all images**: 모든 이미지에 같은 값 (Options의 Thresholds 칸이 Manual로 바뀜)
   - **This image only**: 이 이미지에만 (표의 *Thresholds* 칸에 표시). 개별 값이 일괄 값보다 우선합니다.
   - **Clear this image's value**: 이 이미지의 개별 값 지우기
5. **Run all** 로 다시 분석합니다.

팁
- 여러 이미지를 넘겨 보며 **모든 이미지에 무난한 값**을 고르는 것이 좋습니다. 이미지마다 다른 값을 쓰면 실험군 비교가 공정하지 않을 수 있습니다.
- 사용한 값과 자동 값은 결과 `_results.xlsx` 의 *settings* 시트에 남습니다.
- 각 임계값이 정확히 어디에 쓰이는지: [algorithm.md 6절](algorithm.md#6-임계값-수동-조정과-미리보기)

## 8. 픽셀 크기 (µm/px)
모든 길이(µm)·넓이(µm²) 값이 픽셀 크기에 달려 있습니다.
- 파일에 정보가 있으면 자동으로 읽어 표의 **µm/px** 칸에 보여 줍니다 (마우스를 올리면 출처 표시).
- 없거나 이상한 값(예: 72 dpi 화면 해상도)이면 **빨간 `?`** 가 표시됩니다.
  - **Run all** 을 누르면 그 세트들의 값을 한 번 묻습니다. 현미경 획득 설정(예: 대물렌즈, zoom, 픽셀 수)에서 확인한 µm/px를 입력하세요.
  - 또는 `?` 칸을 더블클릭해 세트마다 직접 입력합니다.
- Options의 *Pixel size* 에 값을 넣으면 파일 정보와 상관없이 모든 이미지에 그 값을 씁니다.
- 확인: 결과 `_results.xlsx` → *settings* 시트의 `pixel_size_um` 과 `pixel_size_source` (TIFF / CZI / entered / assumed).
- Fiji에서 파일에 저장하려면: Image → Properties… 에서 Unit = micron, Pixel width/height 입력 후 TIFF로 저장.

## 9. 그룹 만들기
그룹 = 비교하고 싶은 실험군 (예: Control, PRD 24h, Drug A …).
- **처음 넣을 때**: 폴더 이름이 기본 그룹입니다. 조건마다 폴더를 나눠 두면 따로 정할 필요가 없습니다.
- **Run 탭에서**: *Group* 칸 더블클릭, 또는 여러 줄 선택 → **Set group…**
- **분석이 끝난 뒤**: Analysis 탭의 **Groups…** 버튼. 표에서 Group 칸을 고치고 **Save** 하면 다시 분석하지 않고 바로 반영됩니다.
  ![Groups 창](manual_img/11_groups_dialog.png)
- 그룹 정보는 결과 폴더의 `groups.csv` 에 저장됩니다.

## 10. 결과 보기: 샘플별 탭 (Run 탭 오른쪽)
표에서 끝난 줄을 클릭하면 그 이미지 결과가 나옵니다.

| 탭 | 내용 |
|---|---|
| **Cells** | 세포마다 모든 지표 (열 = 세포, `(+)`/`(−)` = green+/−) |
| **Mito objects** | 미토 조각 하나하나의 형태와 그 위 green. 위에서 green+/− 와 세포로 거를 수 있음 |
| **Green puncta** | green 밝은 점 하나하나 (넓이, 밝기, 미토 위인지) |
| **Correlation** | 이 이미지 안에서 미토 지표 × green 지표 상관 히트맵. 칸을 누르면 산점도 |
| **Green on mito** | 미토 조각(분홍 윤곽) 위 green, 미토 위 점(노랑) / 밖 점(하늘색) |
| **Green+ / − cells** | 밝은 green(노랑)과 green+ 세포(빨간 윤곽) / green− (회색) |
| **MiNA overlay**, **Cells (zoom)** | 미토 마스크(분홍), 골격(초록), 끝점(노랑), 분기점(파랑) |
| **Cell ROIs** | 세포 경계와 번호. 빨간 굵은 선 = 분석, 회색 점선 = 제외(edge), 빗금 = 가장자리에서 잘라낸 부분 |

![Green+ / − cells](manual_img/07_green_cells.png)

**확인 순서 (권장)**: Cell ROIs (세포가 맞게 나뉘었나) → Green+ / − cells (발현 세포 판정이 맞나) → MiNA overlay (미토가 잘 잡혔나) → 표와 그래프.

## 11. Analysis 탭: 여러 샘플 통합 분석
배치가 끝나면 그 결과 폴더가 자동으로 불러와집니다. 예전 결과는 **Browse…** 로 결과 폴더를 고르고 **Load**.

위쪽 버튼
- **Groups…**: 샘플의 그룹 바꾸기 (9절)
- **Export by group…**: 그룹별 표와 그룹 비교표를 다시 만들기 (12절). 대조군(control)을 고를 수 있습니다.

### 11-1. Regression / correlation
![Regression / correlation](manual_img/08_analysis_regression.png)
- 왼쪽 히트맵: 모든 (미토 지표 × green 지표) 쌍의 상관계수. 빨강 = 양의 상관, 파랑 = 음의 상관, `–` = 계산 불가 (세포 3개 미만이거나 값이 모두 같음).
- 칸을 누르면 오른쪽에 산점도 + 회귀선(95 % 신뢰 띠), r, R², p, Spearman ρ. 평균(또는 중앙값)으로 4분면을 나눠 각 분면 비율을 표시합니다.
- 위에서 *Level* (세포 / 미토 조각), *Cells* (전체 / green+ / green−), *Group*, *Dataset*, *Sample* 로 거르고, *Colour by* 로 색을 나눕니다. *log X / log Y* 는 로그 축.
- **Export tables…**: 지금 거른 세포·미토 표와 회귀표 저장.

### 11-2. Groups: mean / median (실험군 비교)
![Groups: mean / median](manual_img/09_groups_mean.png)
- *Metric* 에서 지표를 고르면 그룹별 점(세포 또는 이미지)과 평균 ± SD / SEM / 95 % CI 또는 중앙값 ± IQR 막대(*Show*).
- *Unit*: **cells** (세포 하나 = 점 하나) 또는 **images** (이미지 평균 = 점 하나). 결론은 images에서도 확인하세요 (13절).
- *Control*: 대조군을 고르면 각 그룹을 대조군과 비교하고 (Mann-Whitney, Holm 보정) 그래프에 `*` p<0.05, `**` p<0.01, `***` p<0.001.
- *Groups shown*: 보여줄 그룹 선택. *Table: this metric only*: 표에 지금 지표만.
- 아래 표: 지표 × 그룹마다 n, 평균, SD, SEM, 95 % CI, 중앙값, 사분위, 검정 p 값. 줄을 누르면 그 지표가 그래프에 나옵니다.
- *Export…* → **Plotted values (CSV)…**: 지금 그래프에 찍힌 점만 내보냅니다. 보이는 그룹과 현재 Level / Cells / Unit 설정을 따르고, 열은 group, dataset, preset, sample, cell(또는 mito), green_status(Unit = images이면 n_rows)와 그래프의 지표 하나뿐입니다. 파일 이름 기본값은 `<지표>_<단위>.csv`. **Statistics table (CSV)…** 는 아래 표를 내보냅니다.
- 그래프는 항상 4:3 비율로 창 크기에 맞춰집니다. 그룹이 많으면 이름의 공통 앞부분(예: `Condition_`)을 빼고 가운데를 `…` 로 줄여 글자가 겹치지 않게 합니다. 전체 이름은 아래 표에 있습니다.

### 11-3. Groups: correlation
![Groups: correlation](manual_img/10_groups_corr.png)
- *Mode* **A vs B**: 그룹 A, B의 상관 히트맵과 차이(B − A). 차이가 유의하면 `*`/`**` (Fisher z 검정). ◀ ▶ 로 B를 차례로 바꿉니다.
- *Mode* **All groups**: 그룹마다 작은 히트맵 + "그룹마다 상관이 다른 칸" 지도 (Cochran's Q).
- 그래프 창은 4:3 비율을 유지합니다. A vs B는 위에 A·B, 아래에 Δ와 열(X) 지표 목록. 그룹이 많으면 작은 히트맵을 4:3에 맞는 격자로 배치합니다.
- 아래에 차이가 큰 쌍 목록과 "몇 칸 중 몇 칸이 p<0.05 (우연 기대 수)" 가 나옵니다.

### 11-4. Cells / Mito objects (all samples)
모든 샘플의 세포·미토 표 (그룹·dataset·sample 열 포함). 열 머리를 누르면 정렬됩니다.

## 12. 출력 파일
```
<출력 폴더>/
  groups.csv                         샘플 → 그룹
  groups/                            그룹별 출력 (배치 후 자동)
    groups_overview.csv              그룹별 이미지·세포 수
    group_comparison.xlsx            지표 × 그룹 비교표 (전체/green+/green−, 세포/이미지/미토 단위)
    <그룹>/<그룹>_cells.csv, _mito.csv, _per_image.csv, _correlations.csv, _results.xlsx …
  <dataset>/<preset>-<dataset>/
    <preset>-<dataset>_all_cells.csv …   그 폴더 샘플을 합친 표
    <sample>/                        이미지 하나의 결과
      <sample>_results.xlsx          모든 표 (green_pos / green_neg / all 시트, settings 시트)
      <sample>_per_cell.csv          세포 표   (+ _green_pos / _green_neg 버전)
      <sample>_per_mito.csv          미토 조각 표
      <sample>_green_puncta.csv      green 점 표
      <sample>_correlations.csv      지표 쌍 상관표
      <sample>_RoiSet.zip            Fiji ROI Manager용 세포·핵 ROI (_RoiSet_filtered.zip = 분석 세포만)
      <sample>_mito_RoiSet.zip, _green_puncta_RoiSet.zip
      <sample>_summary.png, _green_cells.png, _green_on_mito.png, _mina_overlay.png …  확인용 그림
```
- Excel에서 보려면 `_results.xlsx` 또는 `groups/group_comparison.xlsx` 를 여세요.
- Fiji에서 ROI 확인: 원본 이미지를 열고 `_RoiSet.zip` 을 Fiji 창에 끌어다 놓으면 ROI Manager에 들어갑니다.
- *settings* 시트: 픽셀 크기, 사용한 임계값과 자동 값, 모든 옵션. 결과를 재현할 때 필요합니다.

## 13. 결과를 해석하는 법
- 지표 하나하나의 뜻과 "조각화가 심하면 어떤 지표가 오르고 내리는지": [metrics.md](metrics.md)
- **세포 수와 이미지 수**: 세포 단위 p 값은 같은 이미지 세포들이 서로 비슷해서 실제보다 작게 나오기 쉽습니다. **Unit = images** 에서도 같은 결론이 나오는지 확인하세요. 그룹당 이미지 5장 이상을 권장합니다.
- **미토 조각 단위**: 같은 세포의 조각들은 독립이 아니라 p 값이 매우 낙관적입니다. 경향을 보는 용도로만 쓰세요.
- **여러 검정**: 히트맵의 240칸을 한꺼번에 보면 약 5 %(12칸)는 우연히 p<0.05입니다. 한 칸보다 **같은 방향으로 묶여 나오는 패턴**을 보세요.
- **발현량**: green 세기 지표(평균, integrated)는 발현량에 크게 좌우됩니다. 그룹 사이 발현량이 다르면 비율 지표(fraction on mito, enrichment, Pearson, Manders)로 비교하는 것이 공정합니다. green+ 세포만 고르면 발현 세포끼리 비교됩니다.
- **결론 전에**: Cell ROIs 그림으로 세포 분할을 한 번 확인하세요.

## 14. 명령줄(터미널)로 실행
화면 없이 많은 데이터를 돌리거나 서버에서 쓸 때.
```bash
# 앱 안의 실행 파일
"/Applications/Mito Analyzer.app/Contents/MacOS/MitoAnalyzer" --cli --batch "내 실험/Control" -o ~/results --group Control
"/Applications/Mito Analyzer.app/Contents/MacOS/MitoAnalyzer" --cli --batch "내 실험/PRD 24h" -o ~/results --group PRD --control Control
# 소스에서
cd mito-tracking/application && .venv/bin/python run_app.py --cli --batch 폴더 -o 출력폴더
# CZI 파일 (파일 하나, 또는 --batch 로 폴더). 추정한 채널 역할이 틀리면 직접 지정
"/Applications/Mito Analyzer.app/Contents/MacOS/MitoAnalyzer" --cli --batch CZI폴더 -o ~/results --czi-channels mito=0,protein=1,nucleus=2
```
자주 쓰는 옵션: `--group 이름`, `--control 대조군`, `--pixel-size-um 0.099`, `--thr-nuclei`, `--cell-fg-level`, `--thr-mito`, `--thr-green-bright`, `--thr-puncta` (0 = 자동), `--mito-method otsu`, `--no-edge-trim`, `--include-edge-cells`, `--green-pos-percent 2`. 전체 목록: `--cli --help`.
같은 `-o` 폴더에 여러 번 돌리면 그룹이 합쳐지고, 마지막 실행 때 `groups/` 비교표가 다시 만들어집니다.

## 15. 문제 해결 (FAQ)
**앱이 열리지 않아요 / "확인되지 않은 개발자"** → 2-1의 우클릭 → 열기, 또는 시스템 설정 → 개인정보 보호 및 보안 → 그래도 열기.

**파일을 넣었는데 표에 안 나와요** → 파일 이름이 `…Red`, `…Green`, `…Blue` 로 끝나는 세트가 세 개 모두 있어야 합니다. 세트가 안 되는 파일은 목록으로 알려 줍니다. `Merged` 같은 다른 TIFF는 무시됩니다. CZI 파일은 채널이 3개 이상이고 읽을 수 있어야 합니다.

**CZI 파일이 엉뚱한 채널로 분석돼요** → 표의 Red / Green / Blue 칸을 확인하고 **Channel roles…** 로 고치세요 (3-1).

**µm/px 에 `?` 가 떠요** → 8절. 값을 입력하거나 Run all 때 물어볼 때 입력하세요.

**세포 두 개가 한 ROI로 잡혀요** → ROI 검토 창(6-1절)에서 Subtract + New cell로 바로 고칠 수 있습니다. 여러 이미지에서 반복되면 한쪽 핵을 못 찾은 경우가 대부분입니다. 미리보기에서 *Nuclei* 를 낮춰 보고, *Dim nuclei* 옵션이 켜져 있는지 확인하세요.

**세포가 아닌 곳에 ROI가 생겨요 (가짜 세포)** → 세포질의 매끈한 파란 덩어리가 핵으로 잡힌 경우. *Nuclei* 를 높이거나 Nucleus detection이 *Texture + merge* 인지 확인하세요.

**ROI가 배경까지 넓어요 / 너무 작아요** → *Cell area* 를 높이거나(넓을 때) 낮춥니다(작을 때). 0.12가 기본.

**분석된 세포가 너무 적어요** → 대부분 가장자리에 닿아 빠진 세포입니다. *Keep cells whose tip touches the border* 가 켜져 있는지 확인하고, 세포가 화면 안에 들어오게 촬영하는 것이 가장 좋습니다.

**green+ 세포가 하나도 없어요** → *Green+ cell: bright green ≥* 를 낮추거나 미리보기에서 *Green+ (bright green)* 를 낮춥니다. Green+ / − cells 탭에서 노란 영역을 확인하세요.

**표에 nan 이나 – 가 보여요** → 계산할 수 없는 값입니다 (세포 3개 미만, 값이 모두 같음, 미토가 없는 세포 등). 화면에 이유가 함께 표시됩니다.

**같은 데이터를 다른 설정으로 비교하고 싶어요** → 표의 *Preset* 이름을 바꿔(예: `thr30`) 돌리면 결과가 다른 폴더에 저장됩니다.

**업데이트는 어떻게?** → DMG: 새 DMG로 다시 설치. 한 줄 설치: 같은 줄 다시 실행. 직접 빌드: `git pull` 후 `./build_mac.sh`.

**버그를 알리고 싶어요** → 문제가 된 이미지의 `_results.xlsx` (settings 시트)와 Log 내용을 함께 보내 주세요.
