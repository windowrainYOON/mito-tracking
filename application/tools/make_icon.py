"""Draw the app icon (cell with mitochondria and green puncta) -> resources/icon.png + icon.icns.
Run: python tools/make_icon.py"""
import math, os, random

import numpy as np
from PIL import Image, ImageDraw, ImageFilter

S = 4096  # supersampled canvas, downscaled to 1024 at the end
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'resources')


def squircle_mask(size, margin, n=5.0):
    y, x = np.mgrid[0:size, 0:size].astype(float)
    r = size / 2 - margin
    u, v = np.abs(x - size / 2) / r, np.abs(y - size / 2) / r
    return ((u ** n + v ** n) <= 1).astype(np.uint8) * 255


def capsule(draw, cx, cy, length, width, angle, fill, outline=None, ow=0):
    """Rounded, slightly bent mitochondrion drawn as one polygon (clean edges at any size)."""
    ts = np.linspace(-0.5, 0.5, 60)
    bend = 0.18 * length * (ts * ts - 0.25)
    xs = cx + ts * length * math.cos(angle) - bend * math.sin(angle)
    ys = cy + ts * length * math.sin(angle) + bend * math.cos(angle)
    pts = list(zip(xs, ys))

    def outline_poly(half):
        dx, dy = np.gradient(xs), np.gradient(ys)
        n = np.hypot(dx, dy); nx, ny = -dy / n, dx / n
        left = list(zip(xs + nx * half, ys + ny * half))
        right = list(zip(xs - nx * half, ys - ny * half))
        a1 = math.atan2(dy[-1], dx[-1]); a0 = math.atan2(dy[0], dx[0])
        cap1 = [(xs[-1] + half * math.cos(a1 + math.pi / 2 - t), ys[-1] + half * math.sin(a1 + math.pi / 2 - t))
                for t in np.linspace(0, math.pi, 24)]
        cap0 = [(xs[0] + half * math.cos(a0 - math.pi / 2 - t), ys[0] + half * math.sin(a0 - math.pi / 2 - t))
                for t in np.linspace(0, math.pi, 24)]
        return left + cap1 + right[::-1] + cap0

    if outline:
        draw.polygon(outline_poly(width / 2 + ow), fill=outline)
    draw.polygon(outline_poly(width / 2), fill=fill)
    return pts


def main():
    os.makedirs(OUT, exist_ok=True)
    # background: deep navy radial gradient
    y, x = np.mgrid[0:S, 0:S].astype(float)
    d = np.hypot(x - S * 0.42, y - S * 0.36) / (S * 0.85)
    top, bot = np.array([34, 52, 98]), np.array([8, 12, 32])
    bg = (top * (1 - d[..., None]).clip(0, 1) + bot * d[..., None].clip(0, 1)).astype(np.uint8)
    img = Image.fromarray(bg).convert('RGBA')

    # cell body: soft translucent membrane
    cell = Image.new('RGBA', (S, S), (0, 0, 0, 0))
    cd = ImageDraw.Draw(cell)
    cx, cy, rx, ry = S * 0.5, S * 0.52, S * 0.36, S * 0.32
    cd.ellipse([cx - rx, cy - ry, cx + rx, cy + ry], fill=(120, 170, 255, 38), outline=(170, 210, 255, 150), width=int(S * 0.012))
    img = Image.alpha_composite(img, cell.filter(ImageFilter.GaussianBlur(S * 0.002)))

    # nucleus: blue with glow
    nuc = Image.new('RGBA', (S, S), (0, 0, 0, 0))
    nd = ImageDraw.Draw(nuc)
    nx, ny, nr = S * 0.43, S * 0.47, S * 0.12
    for i in range(40):  # radial shading, lighter towards the upper left
        f = 1 - i / 40
        c = (int(55 + 70 * (1 - f)), int(100 + 70 * (1 - f)), 255, 255)
        ox, oy = -nr * 0.25 * (1 - f), -nr * 0.25 * (1 - f)
        nd.ellipse([nx + ox - nr * 1.15 * f, ny + oy - nr * f, nx + ox + nr * 1.15 * f, ny + oy + nr * f], fill=c)
    glow = nuc.filter(ImageFilter.GaussianBlur(S * 0.03))
    img = Image.alpha_composite(img, glow)
    img = Image.alpha_composite(img, nuc.filter(ImageFilter.GaussianBlur(S * 0.004)))

    # mitochondria: magenta capsules with cristae, around the nucleus
    mito = Image.new('RGBA', (S, S), (0, 0, 0, 0))
    md = ImageDraw.Draw(mito)
    rnd = random.Random(7)
    specs = [(0.64, 0.33, 0.20, -0.5), (0.72, 0.52, 0.22, 1.35), (0.60, 0.70, 0.19, 0.35), (0.36, 0.73, 0.17, -0.25),
             (0.25, 0.55, 0.15, 1.6), (0.33, 0.30, 0.14, 0.55), (0.52, 0.27, 0.10, 0.1), (0.78, 0.70, 0.10, -1.0)]
    puncta_sites = []
    for fx, fy, fl, ang in specs:
        w = S * 0.045
        pts = capsule(md, S * fx, S * fy, S * fl, w, ang, fill=(235, 70, 170, 255), outline=(255, 150, 215, 255), ow=S * 0.006)
        for t in np.linspace(0.15, 0.85, 6):  # cristae
            p = pts[int(t * (len(pts) - 1))]; q = pts[min(int(t * (len(pts) - 1)) + 1, len(pts) - 1)]
            a = math.atan2(q[1] - p[1], q[0] - p[0]) + math.pi / 2
            h = w * 0.26
            md.line([(p[0] - h * math.cos(a), p[1] - h * math.sin(a)), (p[0] + h * math.cos(a), p[1] + h * math.sin(a))],
                    fill=(150, 30, 110, 255), width=int(S * 0.006))
        puncta_sites += [pts[rnd.randint(5, len(pts) - 6)] for _ in range(2 if fl > 0.15 else 1)]
    shadow = mito.filter(ImageFilter.GaussianBlur(S * 0.012))
    shadow = Image.fromarray((np.array(shadow) * np.array([0, 0, 0, 0.55])).astype(np.uint8))
    img = Image.alpha_composite(img, shadow)
    img = Image.alpha_composite(img, mito)

    # green POI puncta sitting on mitochondria
    pun = Image.new('RGBA', (S, S), (0, 0, 0, 0))
    pd = ImageDraw.Draw(pun)
    for (px, py) in puncta_sites:
        r = S * rnd.uniform(0.016, 0.023)
        pd.ellipse([px - r, py - r, px + r, py + r], fill=(90, 255, 120, 255))
    img = Image.alpha_composite(img, pun.filter(ImageFilter.GaussianBlur(S * 0.018)))
    img = Image.alpha_composite(img, pun.filter(ImageFilter.GaussianBlur(S * 0.002)))

    # glossy top highlight
    hl = Image.new('RGBA', (S, S), (0, 0, 0, 0))
    ImageDraw.Draw(hl).ellipse([-S * 0.3, -S * 0.75, S * 1.3, S * 0.42], fill=(255, 255, 255, 16))
    img = Image.alpha_composite(img, hl.filter(ImageFilter.GaussianBlur(S * 0.04)))

    # macOS-style squircle with margin, plus a subtle drop shadow
    margin = int(S * 0.09)
    mask = Image.fromarray(squircle_mask(S, margin))
    shape = Image.new('RGBA', (S, S), (0, 0, 0, 0)); shape.paste(img, (0, 0), mask)
    sh = Image.new('RGBA', (S, S), (0, 0, 0, 0)); sh.paste((0, 0, 0, 110), (0, int(S * 0.012)), mask)
    final = Image.alpha_composite(sh.filter(ImageFilter.GaussianBlur(S * 0.012)), shape)
    final = final.resize((1024, 1024), Image.LANCZOS)
    final.save(os.path.join(OUT, 'icon.png'))
    final.save(os.path.join(OUT, 'icon.icns'), sizes=[(16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)])
    final.save(os.path.join(OUT, 'icon.ico'), sizes=[(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
    print('wrote', OUT)


if __name__ == '__main__':
    main()
