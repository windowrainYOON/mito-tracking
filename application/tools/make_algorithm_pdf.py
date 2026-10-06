"""Step-by-step algorithm explanation PDF, drawn from one example image set.

    python tools/make_algorithm_pdf.py RED.tif GREEN.tif BLUE.tif -o docs/Mito_Analyzer_algorithm.pdf

Runs the real pipeline on the set (into a temporary folder) and recomputes the intermediate images of
every step with the same functions, so the figures always match the current code.
"""
import argparse
import os
import sys
import tempfile

import numpy as np
from scipy import ndimage as ndi
from skimage import filters, morphology, segmentation

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib  # noqa: E402
matplotlib.use('Agg')
from matplotlib import font_manager, pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

from mito_app import __version__, cell_roi, green as green_q, green_cells, mina, mito_objects, nuclei, pipeline, stats  # noqa: E402

KO_FONTS = ('AppleGothic', 'Apple SD Gothic Neo', 'Malgun Gothic', 'NanumGothic', 'Noto Sans CJK KR',
            'WenQuanYi Zen Hei')


def set_font():
    have = {f.name for f in font_manager.fontManager.ttflist}
    for f in KO_FONTS:
        if f in have:
            plt.rcParams['font.family'] = f
            break
    plt.rcParams['axes.unicode_minus'] = False


def norm(x, p=99.5):
    x = x.astype(float)
    return np.clip(x / max(np.percentile(x, p), 1e-9), 0, 1)


def rgb(r=None, g=None, b=None):
    shape = next(c for c in (r, g, b) if c is not None).shape
    z = np.zeros(shape)
    return np.dstack([z if c is None else c for c in (r, g, b)])


def outline(img, lab, color=(0, 1, 1), mode='inner'):
    out = img.copy()
    out[segmentation.find_boundaries(lab, mode=mode)] = color
    return out


def page(pdf, title, text, panels, cols=None):
    """One A4-landscape page: title, explanation on the left, image panels (title, image[, cmap]) on the right."""
    fig = plt.figure(figsize=(11.69, 8.27))
    fig.text(0.04, 0.95, title, fontsize=17, weight='bold', va='top')
    fig.text(0.03, 0.88, text, fontsize=8.6, va='top', linespacing=1.55, wrap=True)
    n = len(panels)
    cols = cols or (1 if n == 1 else 2 if n <= 4 else 3)
    rows = int(np.ceil(n / cols))
    x0, w, y0, h = 0.42, 0.57, 0.04, 0.86
    cw, ch = w / cols, h / rows
    for i, pnl in enumerate(panels):
        r_, c_ = divmod(i, cols)
        ax = fig.add_axes([x0 + c_ * cw + 0.005, y0 + (rows - 1 - r_) * ch + 0.005, cw - 0.01, ch - 0.035])
        if callable(pnl[1]):
            pnl[1](ax)
        else:
            ax.imshow(pnl[1], cmap=pnl[2] if len(pnl) > 2 else None, interpolation='nearest')
            ax.set_axis_off()
        if pnl[0]:
            ax.set_title(pnl[0], fontsize=9)
    pdf.savefig(fig)
    plt.close(fig)


def label_cells(ax, lab, names=None, color='white'):
    for l in np.unique(lab[lab > 0]):
        cy, cx = ndi.center_of_mass(lab == l)
        ax.text(cx, cy, names.get(l, str(l)) if names else str(l), color=color, fontsize=7, ha='center', va='center',
                weight='bold', bbox=dict(facecolor='black', alpha=0.5, pad=1, edgecolor='none'))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('red'); ap.add_argument('green'); ap.add_argument('blue')
    ap.add_argument('-o', '--out', default='docs/Mito_Analyzer_algorithm.pdf')
    A = ap.parse_args(argv)
    set_font()
    tmp = tempfile.mkdtemp()
    res_run = pipeline.run(A.red, A.green, A.blue, tmp, log=lambda s: None)
    rows = {r['cell']: r for r in res_run['rows']}
    ok_ids = [int(r['roi'][4:]) for r in res_run['rois'] if r['status'] == 'ok']

    # ---- recompute intermediates with the pipeline's own functions ----
    a, g, px = cell_roi.load_green(A.green)
    px = px if px != 1.0 else cell_roi.REF_PX_UM
    k = float(np.clip(cell_roi.REF_PX_UM / px, 0.2, 5))
    nuc_img, _ = cell_roi.load_nuclei(A.blue, k)
    nuc0, dim, ninfo = nuclei.detect(nuc_img, k, pipeline.Params().nuclei_method)
    red, _ = mina.load_channel(A.red, 0)
    green, _ = mina.load_channel(A.green, 1)
    gn, nrm = cell_roi.normalize_green(g, k)
    n_focus = int(nuc0.max())
    nuc = np.where(dim > 0, dim, nuc0)
    seg = cell_roi.segment(gn, nuclei=nuc, k=k)
    morph = cell_roi.segment_morph(red, nuc, seg['landscape'], px, k)
    lab = morph['labels']
    md, ridge, fg = morph['mito_density'], morph['ridge'], morph['fg']
    cost = (1 - md) + ridge + seg['landscape']
    gstatus, ginfo = green_cells.classify(gn, lab, nuc, k)
    R, Gc, B = norm(red), norm(green), norm(nuc_img)
    um = lambda v: v / px  # noqa: E731

    # example cell for the zoomed steps: the largest analysed one
    ex = max(ok_ids, key=lambda c: (lab == c).sum())
    ys, xs = np.nonzero(lab == ex)
    pad = int(um(3))
    sl = (slice(max(0, ys.min() - pad), ys.max() + pad), slice(max(0, xs.min() - pad), xs.max() + pad))
    cell = lab == ex
    p0 = dict(mito_objects.P, sigma=mito_objects.P['sigma'] * k)
    sm = mito_objects.preprocess(red, px, p0)
    q = mito_objects.scaled_params(sm, np.isin(lab, ok_ids), k, p0)
    mfg, mlab, msep = mito_objects.segment_cell(sm, cell, px, q)
    msk = morphology.skeletonize(msep > 0, method='lee')
    plab, pthr, tophat = green_q.detect_puncta(np.clip(green.astype(float) - nrm['offset'], 0, None),
                                               np.isin(lab, ok_ids), k)
    exname = f'cell{ex:02d}'

    os.makedirs(os.path.dirname(os.path.abspath(A.out)), exist_ok=True)
    with PdfPages(A.out) as pdf:
        # 0. overview
        page(pdf, f'Mito Analyzer {__version__} — 알고리즘 단계별 설명',
             '세포 하나하나의 미토콘드리아 형태(MiNA)와 관심 단백질(green)이 미토콘드리아 위에\n'
             '얼마나 있는지를 측정하고, 둘의 관계를 회귀분석으로 보는 프로그램입니다.\n\n'
             '입력: 같은 시야의 3채널 이미지\n'
             '  - Red  = 미토콘드리아 (MitoTracker 등)\n'
             '  - Green = 관심 단백질 (POI)\n'
             '  - Blue = 핵 (DAPI/Hoechst)\n\n'
             '처리 순서\n'
             '  0. 밝기 자동 보정 (이미지마다 gain·노출이 달라도 같은 기준)\n'
             '  1. 세포 ROI: 핵 하나 = 세포 하나, 미토콘드리아가\n'
             '     끊어지는 선을 경계로 사용\n'
             '  2. green+ / green- 세포 분류 (ROI 모양은 바꾸지 않음)\n'
             '  3. 세포마다 미토콘드리아 분할 + MiNA 형태 분석\n'
             '  4. 미토콘드리아 위의 green 신호, green puncta\n'
             '  5. 표·이미지 출력, 상관/회귀분석, 여러 샘플 통합 분석\n\n'
             f'예시 이미지: {os.path.basename(A.red)} 외 2채널\n'
             f'{red.shape[1]}×{red.shape[0]} px, 픽셀 {px:.4f} µm\n'
             f'모든 크기 상수는 µm 기준이라 배율이 달라도 같은 크기로 동작합니다.',
             [('Red: 미토콘드리아', rgb(R, R * 0.55)), ('Green: 관심 단백질', rgb(None, Gc)),
              ('Blue: 핵', rgb(None, None, B)), ('Merge', rgb(R, Gc, B))])

        # 1. auto-levels
        sm4 = ndi.gaussian_filter(g, 4 * k)

        def hist(ax):
            ax.hist(sm4.ravel(), bins=200, color='#4a9a4a')
            for (v, t), fy in zip(((np.percentile(sm4, 1), '1%: 배경'), (np.percentile(sm4, 75), '75%: 기준 밝기')), (0.92, 0.80)):
                ax.axvline(v, color='k', ls='--', lw=1)
                ax.text(v, fy, ' ' + t, fontsize=8, transform=ax.get_xaxis_transform())
            ax.set_yscale('log'); ax.set_xlabel('green (smoothed)'); ax.tick_params(labelsize=7)
        page(pdf, '0단계 — 밝기 자동 보정 (auto-levels)',
             'green 이미지를 σ = 4 px(기준 픽셀)로 부드럽게 한 뒤\n'
             '  - 하위 1 % 값 = 세포 사이 어두운 틈 = 배경(offset)\n'
             '  - 75 % 값 - 배경 = 세포질 기준 밝기(unit)\n'
             '로 보고, 기준 세트와 같은 단위가 되도록 선형 변환합니다.\n\n'
             f'이 이미지: 배경 {nrm["offset"]:.2f}, gain ×{nrm["gain"]:.2f}\n\n'
             '이 덕분에 노출·gain·bit depth가 다른 이미지에서도 같은\n'
             '임계값과 상수가 그대로 쓰입니다. 크기 상수는 픽셀 크기\n'
             '비율 k = 0.099 µm / 픽셀 크기로 함께 조정됩니다.\n\n'
             'green의 약한 자가형광은 세포 경계의 보조 정보로도 쓰입니다\n'
             '(1단계의 비용 지도).',
             [('원본 green (표시용 자동 대비)', Gc, 'gray'), ('보정된 green', np.clip(gn / 20, 0, 1), 'gray'),
              ('밝기 분포와 기준점', hist)])

        # 2. nuclei
        nimg = rgb(B * 0.3, B * 0.3, B)
        nimg = outline(nimg, ninfo['rejected'], (0, 0.9, 1))
        nimg = outline(nimg, nuc0, (1, 1, 0)); nimg = outline(nimg, dim, (1, 0.3, 1))
        sflat = nuclei.flatten(filters.gaussian(nuc_img, 2 * k, preserve_range=True), k, px=px)
        tex = np.zeros(nuc_img.shape)
        for l, (t, _) in ninfo['texture'].items():
            tex[ninfo['regions'] == l] = t
        page(pdf, '1-1단계 — 핵 찾기 (세포의 씨앗)',
             '세포 하나에는 핵이 하나 있고, 세포질은 핵을 둘러싸고 있다는\n'
             '형태를 그대로 이용합니다. 핵이 각 세포 ROI의 씨앗(seed)입니다.\n\n'
             '핵 염색에는 세포질에 퍼진 흐린 blue, 잡음, 핵 안의 어두운 구멍이\n'
             '섞여 있습니다. 진짜 핵은 염색질 무늬(texture)가 있고, 세포질의\n'
             'blue는 매끈합니다. 이 차이로 핵을 고릅니다.\n\n'
             '  1. σ = 2 px로 부드럽게, 25 µm 배경 제거 → Otsu 값의 50 %로 후보\n'
             '  2. 거리 지도 watershed로 후보를 조각으로 나눔\n'
             '  3. 조각마다 무늬 값(잡음을 뺀 표준편차 / 중앙값) 계산:\n'
             '     ≥ 0.22 이거나, ≥ 0.17 이면서 이 이미지의 무늬 있는 핵만큼\n'
             '     밝으면(중앙값 ≥ 0.7 × 기준) 남김. 나머지는 버림 (하늘색)\n'
             '  4. 남은 이웃 조각은 합침. 접촉선이 어둡고(경계 < 0.8 × 내부)\n'
             '     좁을 때만 따로 둠 → 휜 핵은 하나로, 붙은 두 핵은 둘로\n'
             '  5. 핵마다 자기 임계값(핵 중심과 주변의 중간)으로 윤곽을 다시\n'
             '     잡고, 가장자리 홈을 메움. 둥글지 않으면(solidity < 0.8) 버림\n'
             '  6. 밝기가 기준의 50 % 이상이고 충분히 크면 초점 맞은 핵(노랑),\n'
             '     아니면 흐린 핵(분홍, 최소 약 40 µm²). 흐린 핵도 자기 세포를\n'
             '     가집니다 (표의 seed 열 = "dim nucleus").\n\n'
             f'이 이미지: 초점 맞은 핵 {n_focus}개, 흐린 핵 {int(nuc.max()) - n_focus}개, '
             f'버린 매끈한 영역 {ninfo["n_rejected"]}개',
             [('A. 배경 제거한 blue', norm(sflat), 'gray'),
              ('B. 후보 영역의 무늬 값 (밝을수록 무늬 강함)', np.clip(tex / 0.4, 0, 1), 'magma'),
              ('C. 결과 (노랑: 초점 핵, 분홍: 흐린 핵, 하늘색: 버린 영역)', nimg)])

        # 3. mito density / ridges / fg
        page(pdf, '1-2단계 — 미토콘드리아 밀도와 "끊어지는 선"',
             '세포와 세포 사이에는 미토콘드리아가 없는 가는 선이 있습니다.\n'
             'red 이미지에서 이 선을 찾아 경계로 씁니다.\n\n'
             'A. 미토콘드리아 밀도: red를 σ = 1.5 µm로 부드럽게 →\n'
             '   상위 1 % 값으로 나눔 (0-1). 미토콘드리아가 빽빽한 곳일수록 밝음.\n\n'
             'B. 끊어지는 선: 밀도 지도에 Sato 필터(어두운 능선, 폭 1.5-2.5 µm)\n'
             '   → 주변보다 가늘고 어두운 선(= 세포 사이 틈)이 밝게 나옴.\n\n'
             'C. 전경(세포가 있는 곳): red를 σ = 3 µm로 부드럽게 한 값이\n'
             '   상위 1 % 값의 12 % 이상인 곳 + 핵. 작은 틈은 닫고(2 µm),\n'
             '   200 µm² 이하 구멍은 채웁니다. 배경은 어떤 세포에도 안 들어감.',
             [('A. 미토콘드리아 밀도', md, 'magma'), ('B. 끊어지는 선 (Sato, 어두운 능선)', np.clip(ridge, 0, 1), 'viridis'),
              ('C. 전경 (흰색) + 핵', rgb(fg * 0.8 + (nuc > 0) * 0.2, fg * 0.8, fg * 0.8 + (nuc > 0) * 0.2)),
              ('red 원본', R, 'gray')])

        # 4. cost + watershed
        def roi_panel(ax):
            img = outline(rgb(R, R * 0.55, B * 0.8), lab)
            ax.imshow(img); ax.set_axis_off()
            label_cells(ax, lab, {l: f'{l}' + ('*' if l > n_focus else '') + ('' if l in ok_ids else ' (edge)')
                                  for l in np.unique(lab[lab > 0])})
        page(pdf, '1-3단계 — 세포 ROI: 핵에서 고르게 퍼지는 watershed',
             '비용 지도 = (1 - 미토콘드리아 밀도) + 끊어지는 선 + green 경계 지도\n'
             '  - 미토콘드리아가 빽빽한 곳은 비용이 낮아 쉽게 넘어가고,\n'
             '  - 끊어지는 선 위는 비용이 높아 경계가 그 선에 놓입니다.\n'
             '  - green 경계 지도: 보정된 green의 어두운 골짜기(Sato, 큰 폭)와\n'
             '    green 밀도가 낮은 곳.\n\n'
             '모든 핵(초점 + 흐린 핵)에서 동시에 watershed를 시작하고,\n'
             'compactness(0.003)를 주어 핵에서 멀어질수록 비용이 조금씩\n'
             '늘게 합니다 → 각 세포가 핵을 중심으로 고르게(마름모꼴에 가깝게)\n'
             '자라고, 한 세포가 길게 뻗어 이웃을 삼키지 않습니다.\n\n'
             '마지막으로 ROI를 매끄럽게 다듬고(열림 연산, 핵은 유지),\n'
             '이미지 가장자리에 닿는 세포는 "edge"로 분석에서 뺍니다\n'
             '(옵션으로 포함 가능). * = 흐린 핵에서 자란 세포.\n\n'
             f'이 이미지: 세포 ROI {int(lab.max())}개, 분석 대상 {len(ok_ids)}개',
             [('비용 지도 (밝을수록 넘기 어려움)', np.clip(cost / np.percentile(cost, 99), 0, 1), 'inferno'),
              ('최종 세포 ROI (청록), 번호', roi_panel)], cols=2)

        # 5. green+/-
        def green_panel(ax):
            v = np.clip(gn / max(np.percentile(gn, 99.5), 1e-9), 0, 1)
            img = rgb(v * 0.15, v, v * 0.15)
            img[ginfo['mask']] = (1, 0.95, 0.2)
            ax.imshow(img); ax.set_axis_off()
            for l in np.unique(lab[lab > 0]):
                pos = gstatus[l] == 'positive'
                ax.contour(lab == l, [0.5], colors='#ff5050' if pos else '#c8c8c8', linewidths=1.6 if pos else 0.8)
            label_cells(ax, lab, {l: f"{l}{'+' if gstatus[l] == 'positive' else '-'} "
                                     f"{100 * ginfo['fraction'][l]:.0f}%" for l in np.unique(lab[lab > 0])})
        page(pdf, '2단계 — green+ / green- 세포 분류',
             '관심 단백질은 보통 세포 단위로 켜지거나 꺼집니다.\n\n'
             '밝은 green 마스크 (노랑)\n'
             '  - 보정된 green을 σ = 1.5 px로 부드럽게 → log\n'
             '  - Otsu 두 번: (세포질 vs 신호) → (약한 신호 vs 강한 신호)\n'
             '  - 단, 세포질 기준 밝기의 3배 미만은 절대 밝다고 보지 않음\n'
             '    → green을 발현하는 세포가 없는 이미지에서는 마스크가 거의 없음\n\n'
             '세포 점수 = 세포질(세포 - 핵) 중 밝은 green이 덮는 비율\n'
             '  - 점수 ≥ 2 % (옵션) → green+ (빨간 윤곽), 아니면 green- (회색)\n\n'
             '중요: 이 단계는 ROI 모양을 바꾸지 않습니다. 분류만 하고,\n'
             'green+와 green- 세포는 이후 모든 측정·상관분석을 각각 따로\n'
             '계산해서 _green_pos / _green_neg 파일로 출력합니다.\n\n'
             f'이 이미지: 밝은 green 기준 ≥ {ginfo["threshold"]:.1f} (보정 단위)',
             [('밝은 green (노랑), 세포 번호 + 점수', green_panel)], cols=1)

        # 6. mito segmentation in one cell
        cm = cell[sl]
        objs = np.zeros(cm.shape + (3,))
        rng = np.random.default_rng(1)
        colors = rng.uniform(0.3, 1, (int(msep.max()) + 1, 3)); colors[0] = 0
        objs = colors[msep[sl]]
        sk_img = rgb(norm(red)[sl] * 0.5, norm(red)[sl] * 0.5, norm(red)[sl] * 0.5)
        sk_img[(msep > 0)[sl]] = (0.6, 0.2, 0.6)
        sk_img[msk[sl]] = (0.2, 1, 0.2)
        r_ = rows.get(exname, {})
        page(pdf, f'3단계 — 미토콘드리아 분할과 MiNA (예: {exname})',
             '분석 대상 세포마다 red 채널에서 미토콘드리아를 나눕니다.\n'
             '(기본 "Split objects", 옵션으로 MiNA classic = 세포별 Otsu 하나)\n\n'
             '  1. 배경 제거: rolling ball (반지름 1.5 µm) → Gaussian\n'
             '  2. 전경: 국소 평균(1.1 µm 창)보다 밝고, 세포 Otsu의 50 % 이상.\n'
             '     아주 작은 구멍은 채움(가짜 고리 방지), 점 같은 잡음 제거\n'
             '  3. 과분할: 밝은 봉우리마다 씨앗(h-maxima) → watershed\n'
             '  4. 보수적 재합침: 두 조각이 떨어져 있으려면 접점이 어둡고\n'
             '     (봉우리의 75 % 미만) 동시에 목처럼 좁아야 함\n'
             '  5. 맞닿은 조각 사이에 1 px 틈 → MiNA가 따로 셉니다\n\n'
             'MiNA: 골격화(skeleton) → 가지(branch), 연결점, 끝점, 네트워크,\n'
             '가지 길이, 도넛(고리) 수, footprint(면적) 등. 조각화 지표:\n'
             '100 µm²당 조각 수, 면적가중 평균 크기, form factor 등.\n\n'
             + (f'{exname}: footprint {r_.get("mitochondrial_footprint_um2", float("nan")):.1f} µm², '
                f'네트워크 {r_.get("n_networks", "-")}, 가지 {r_.get("n_branches", "-")}' if r_ else ''),
             [('red (세포 안)', np.where(cm, norm(red)[sl], norm(red)[sl] * 0.25), 'gray'),
              ('배경 제거 후', norm(sm)[sl] * cm, 'gray'),
              ('전경 마스크', mfg[sl], 'gray'), ('분할된 미토콘드리아 (색 = 조각)', objs),
              ('MiNA 골격 (초록) / footprint (보라)', sk_img)], cols=3)

        # 7. green on mito
        gz = norm(green)[sl]
        gimg = rgb(norm(red)[sl] * 0.6, gz, norm(red)[sl] * 0.6)
        gimg[segmentation.find_boundaries(msep[sl] > 0)] = (1, 0.2, 1)
        on = ndi.binary_dilation(msep > 0)
        pimg = rgb(gz * 0.4, gz * 0.4, gz * 0.4)
        pimg[((plab > 0) & on)[sl]] = (1, 1, 0.2)
        pimg[((plab > 0) & ~on)[sl]] = (0.2, 1, 1)
        page(pdf, f'4단계 — 미토콘드리아 위의 green (예: {exname})',
             'green에서 배경(0단계 offset)을 뺀 값으로 측정합니다.\n\n'
             '세포 단위\n'
             '  - 미토콘드리아 위 / 밖 green 평균, 농축비(on/off),\n'
             '    전체 green 중 미토콘드리아 위 비율, 적분값\n'
             '  - red-green Pearson, Manders 계수\n\n'
             'green puncta (점상 신호)\n'
             '  - green을 σ = 1 px로 부드럽게 → white top-hat (반지름 8 px)으로\n'
             '    퍼진 세포질 신호 제거\n'
             '  - 이미지당 임계값 하나 = max(Otsu, 중앙값 + 6 MAD) × 민감도\n'
             '    → 한 이미지 안의 어두운/밝은 세포를 같은 기준으로 비교\n'
             '  - 노랑 = 미토콘드리아 위 puncta, 청록 = 밖\n\n'
             '미토콘드리아 조각 단위\n'
             '  - 조각마다 길이, 가지 수, 종횡비, 면적, solidity, form factor\n'
             '  - 그 위의 green 평균/적분/최대, 겹치는 puncta 수와 면적\n\n'
             f'이 이미지: puncta 임계값 {pthr:.1f} (top-hat 단위)',
             [('green (초록) + red, 미토콘드리아 윤곽 (분홍)', gimg),
              ('green top-hat', norm(tophat)[sl], 'gray'), ('puncta: 노랑 = 미토 위, 청록 = 밖', pimg)], cols=2)

        # 8. statistics
        mrows = [m for m in res_run['mito_rows']]
        x = np.array([pipeline._num(m.get('length_um')) for m in mrows], float)
        x = np.log10(np.where(x > 0, x, np.nan))
        y = np.array([pipeline._num(m.get('green_mean')) for m in mrows], float)
        okm = np.isfinite(x) & np.isfinite(y)
        x, y = x[okm], y[okm]
        reg = stats.regression(x, y)

        def scatter(ax):
            ax.scatter(x, y, s=4, alpha=0.35, color='#3060c0', lw=0)
            mx, my = np.mean(x), np.mean(y)
            ax.axvline(mx, color='gray', lw=0.8, ls='--'); ax.axhline(my, color='gray', lw=0.8, ls='--')
            if np.isfinite(reg['slope']):
                xx = np.linspace(x.min(), x.max(), 50)
                ax.plot(xx, reg['slope'] * xx + reg['intercept'], color='#c03030', lw=1.5)
            for qx, qy, t in ((1, 1, 'I'), (0, 1, 'II'), (0, 0, 'III'), (1, 0, 'IV')):
                m_ = ((x >= mx) == bool(qx)) & ((y >= my) == bool(qy))
                ax.text(0.97 if qx else 0.03, 0.95 if qy else 0.05, f'{t}: {100 * m_.mean():.0f}%',
                        transform=ax.transAxes, ha='right' if qx else 'left', va='top' if qy else 'bottom', fontsize=8)
            ax.set_xlabel('log10 미토콘드리아 조각 길이 (µm)', fontsize=8); ax.set_ylabel('조각 위 green 평균', fontsize=8)
            ax.tick_params(labelsize=7)
            ax.set_title(f"예: 조각 길이 × green, r = {reg['r']:.3f}, p = {reg['p']:.2g}, n = {reg['n']}", fontsize=9)
        page(pdf, '5단계 — 출력과 통계',
             '상관/회귀분석\n'
             '  - 모든 (미토콘드리아 지표 × green 지표) 쌍에 대해 최소제곱 회귀선,\n'
             '    Pearson r, R², p, Spearman ρ를 세포 단위·조각 단위로 계산\n'
             '  - 전체 / green+ 만 / green- 만 각각 따로\n'
             '  - 산점도는 X·Y 평균(또는 중앙값)에서 사분면으로 나눠 각 사분면\n'
             '    비율을 표시 → 음의 상관(II·IV 사분면)도 바로 보입니다\n'
             '  - 점이 3개 미만이거나 값이 모두 같으면 nan 대신 이유를 표시\n\n'
             '여러 샘플 통합 분석 (Analysis 탭)\n'
             '  - 결과 폴더를 지정하면 하위 폴더의 _per_cell.csv / _per_mito.csv를\n'
             '    모두 모아 dataset·preset·sample 열을 붙여 함께 분석\n'
             '  - green+/-, dataset, sample로 거르고 색 구분, 회귀표 내보내기\n\n'
             '샘플마다 출력\n'
             '  - _results.xlsx (green_pos / green_neg / all 시트, settings)\n'
             '  - _per_cell / _per_mito / _green_puncta / _correlations CSV\n'
             '    (각각 _green_pos, _green_neg 버전 포함)\n'
             '  - Fiji용 ROI zip (세포, 미토콘드리아, puncta), 오버레이 이미지\n\n'
             '참고: 조각 단위 점들은 같은 세포 안에서 서로 독립이 아니므로\n'
             'p 값은 세포 단위 결과와 함께 해석하는 것이 안전합니다.',
             [('', scatter)], cols=1)
    print(f'Wrote {A.out}')


if __name__ == '__main__':
    main()
