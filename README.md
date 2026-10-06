# mito-tracking

**Mito Analyzer**: a desktop app that measures, for each cell, mitochondrial morphology (MiNA) and how much of a protein of interest (green) sits on the mitochondria, then looks at how the two relate. The code is in [`application/`](application/), with the full description in [`application/README.md`](application/README.md).

📄 **Step-by-step algorithm explanation with an example image (PDF):** [`application/docs/Mito_Analyzer_algorithm.pdf`](application/docs/Mito_Analyzer_algorithm.pdf)

## 설치 방법

### macOS 앱 (.app)
```bash
git clone https://github.com/windowrainYOON/mito-tracking.git
cd mito-tracking/application
./build_mac.sh --add-to-dock     # application/Mito Analyzer.app 을 만들고 Dock 에 고정
open "Mito Analyzer.app"
```
- Python 3.11 이상이 필요합니다 (`brew install python@3.12` 또는 python.org 설치본).
- 서명하지 않은 로컬 빌드라서 macOS 가 막으면 앱을 우클릭 → 열기를 한 번 해 주세요.
- 업데이트: `git pull` 후 `./build_mac.sh` 를 다시 실행하면 같은 경로에 다시 빌드됩니다.

### 소스에서 바로 실행 (macOS / Windows / Linux)
```bash
cd mito-tracking/application
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python run_app.py                                              # GUI
.venv/bin/python run_app.py --cli RED.tif GREEN.tif BLUE.tif -o OUTDIR    # 이미지 세트 하나, 화면 없이
.venv/bin/python run_app.py --cli --batch FOLDER [FOLDER…] -o OUTDIR      # 폴더 안 모든 세트 (하위 폴더까지)
```

### 입력 이미지
같은 시야의 3채널 TIFF (ImageJ RGB 내보내기 또는 단일 채널 8-bit). 파일 이름으로 세트를 묶습니다: `…_Ch1_Red.tif` (미토콘드리아), `…_Ch2_Green.tif` (관심 단백질), `…_Ch3_Blue.tif` (핵).

## 알고리즘 요약
0. **밝기 자동 보정**: green 의 배경(하위 1 %)과 세포질 기준 밝기(75 %)를 기준 세트 단위로 맞춥니다. 크기 상수는 모두 µm 기준이라 배율·gain·노출이 달라도 같게 동작합니다.
1. **세포 ROI** (형태 + 밝기)
   - 핵 하나가 세포 하나입니다. 초점이 맞은 핵과 흐린 핵 모두 씨앗이 됩니다.
   - 핵 찾기: 세포질에 퍼진 흐린 blue와 진짜 핵을 염색질 무늬(texture)로 구분합니다. 배경을 빼고 낮은 임계값으로 후보를 잡은 뒤, 무늬가 약한(매끈한) 영역은 버리고, 같은 핵의 조각은 합치고(어둡고 좁은 접촉선이 있는 붙은 핵은 따로), 핵마다 자기 임계값으로 윤곽을 다듬습니다. 옵션 *Nucleus detection* 에서 예전 방식(Global Otsu)도 고를 수 있습니다.
   - 세포 경계는 미토콘드리아가 끊어지는 어두운 선을 따라갑니다 (red 밀도 지도에 Sato 어두운 능선 필터). green 자가형광의 골짜기도 보조로 씁니다.
   - compact watershed 로 각 세포가 핵을 중심으로 고르게 자랍니다.
   - 핵은 안쪽에 있고 세포 끝부분만 이미지 가장자리에 닿으면, 미토콘드리아가 끊어지는 선을 따라 그 끝부분을 잘라내고 분석합니다. 그래도 가장자리에 닿는 세포는 분석에서 뺍니다.
2. **green+ / green− 분류**: 밝은 green 이 세포질의 2 % 이상을 덮으면 green+ 입니다. ROI 모양은 바꾸지 않습니다.
3. **미토콘드리아 분할 + MiNA**
   - 세포마다 rolling-ball 배경 제거 → 국소 임계값 → 봉우리 기반 watershed 로 조각을 나누고, 보수적으로 다시 합칩니다.
   - 골격화해서 가지, 네트워크, 길이, 도넛 수와 조각화 지표를 계산합니다.
4. **미토콘드리아 위의 green**
   - 세포 단위: on/off 평균, 농축비, Pearson, Manders.
   - green puncta: top-hat 후 이미지당 하나의 임계값을 씁니다.
   - 미토콘드리아 조각 단위: 형태와 그 위의 green.
5. **출력과 통계**
   - 모든 지표 쌍의 회귀(r, R², p)와 Spearman ρ 를 구하고, 산점도를 사분면으로 나눠 보여 줍니다.
   - green+ 와 green− 는 각각 따로 계산해 출력합니다.
   - *Analysis* 탭에서 결과 폴더 전체를 모아 통합 분석할 수 있습니다.

설명 PDF 는 `cd application && .venv/bin/python tools/make_algorithm_pdf.py RED.tif GREEN.tif BLUE.tif -o docs/Mito_Analyzer_algorithm.pdf` 로 현재 코드 기준으로 다시 만들 수 있습니다.
