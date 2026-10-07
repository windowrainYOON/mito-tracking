"""Image input: TIFF files and Zeiss CZI files.

A CZI file holds all channels of a field in one file. The pipeline works on one image per channel, so a CZI
channel is addressed as '<file>.czi::ch<N>' (N = 0-based channel index in the file); `imread` and
`pixel_size` accept such references as well as plain TIFF paths.

The CZI reader is self-contained (ZISRAW segment format: file header -> subblock directory -> subblocks);
uncompressed subblocks need nothing else, JPEG XR / zstd compressed ones need the optional imagecodecs package.
Per channel it takes the first scene and time point, stitches mosaic tiles, and makes a maximum-intensity
projection when the file has several z planes. Images of more than 8 bits are scaled to 0-255 with one fixed
factor per bit depth (e.g. 12-bit / 16), never per image, so values stay comparable across a batch and the
TIFF-tuned thresholds keep their meaning."""
import os
import re
import struct
import xml.etree.ElementTree as ET
from functools import lru_cache

import numpy as np
import tifffile

CZI_EXT = '.czi'
REF_RX = re.compile(r'^(.*\.czi)::ch(\d+)$', re.I)
# CZI pixel types: (numpy dtype, samples per pixel)
PIXEL_TYPES = {0: ('<u1', 1), 1: ('<u2', 1), 2: ('<f4', 1), 3: ('<u1', 3), 4: ('<u2', 3), 9: ('<u1', 4),
               12: ('<u4', 1)}
# channel roles in the pipeline, keyed as the image-set dicts are ('red' = mitochondria, ...)
ROLES = (('red', 'Mito'), ('green', 'Protein (POI)'), ('blue', 'Nucleus'))


def is_czi(path):
    return str(path).lower().endswith(CZI_EXT)


def channel_ref(path, ch):
    return f'{path}::ch{int(ch)}'


def split_ref(path):
    """'<file>.czi::ch2' -> ('<file>.czi', 2); anything else -> (path, None)."""
    m = REF_RX.match(str(path))
    return (m.group(1), int(m.group(2))) if m else (path, None)


def display_name(path):
    """Short label for a channel file or CZI channel reference ('New-01.czi ch2 DAPI-T3')."""
    f, ch = split_ref(path)
    if ch is None:
        return os.path.basename(path)
    chans = czi_info(f)['channels']
    label = chans[ch]['name'] if ch < len(chans) else ''
    return f'{os.path.basename(f)} ch{ch}' + (f' {label}' if label else '')


def imread(path):
    """Image array of a TIFF file or of one CZI channel (2-D uint8 for CZI)."""
    f, ch = split_ref(path)
    if ch is None:
        return tifffile.imread(path)
    st = os.stat(f)
    return _czi_channel(f, ch, st.st_mtime_ns, st.st_size).copy()


def pixel_size(path):
    """(pixel size in um, source) for a CZI file / channel reference; (None, reason) when it has none."""
    f, _ = split_ref(path)
    try:
        px = czi_info(f)['px_um']
    except Exception as e:
        return None, f'cannot read CZI ({type(e).__name__}: {e})'
    return (px, 'CZI Scaling') if px else (None, 'no Scaling in the CZI metadata')


def czi_info(path):
    """{'channels': [{'name', 'fluor', 'emission'}], 'px_um', 'size': (h, w), 'bits', 'n_z', 'n_t', 'n_scenes'}"""
    st = os.stat(path)
    return _czi_info(path, st.st_mtime_ns, st.st_size)


def guess_roles(channels):
    """Default role -> channel index from the channel names / dyes / emission wavelengths:
    a DNA stain (DAPI, Hoechst, ...) or the shortest emission is the nucleus, then the longest emission is
    mitochondria and the remaining one the protein (the red / green / blue order of the TIFF exports)."""
    n = len(channels)
    if n < 3:
        return {}
    dna = re.compile(r'dapi|hoechst|draq|sytox|nuc|dna', re.I)
    idx = list(range(n))
    em = [c.get('emission') or 0 for c in channels]
    hits = [i for i in idx if dna.search(f"{channels[i].get('name', '')} {channels[i].get('fluor', '')}")]
    if hits:
        blue = hits[0]
    elif all(em):
        blue = min(idx, key=lambda i: em[i])
    else:
        blue = n - 1
    rest = [i for i in idx if i != blue]
    if all(em[i] for i in rest):
        rest.sort(key=lambda i: -em[i])
    return {'red': rest[0], 'green': rest[1], 'blue': blue}


def parse_roles(text):
    """'mito=0,protein=1,nucleus=2' (or red=/green=/blue=) -> {'red': 0, 'green': 1, 'blue': 2}."""
    alias = {'mito': 'red', 'red': 'red', 'protein': 'green', 'poi': 'green', 'green': 'green',
             'nucleus': 'blue', 'nuclei': 'blue', 'dapi': 'blue', 'blue': 'blue'}
    out = {}
    for part in filter(None, (s.strip() for s in text.split(','))):
        k, _, v = part.partition('=')
        if k.strip().lower() not in alias or not v.strip().isdigit():
            raise ValueError(f"bad channel role '{part}' (use e.g. mito=0,protein=1,nucleus=2)")
        out[alias[k.strip().lower()]] = int(v)
    if set(out) != {'red', 'green', 'blue'}:
        raise ValueError('give all three channel roles: mito=…,protein=…,nucleus=…')
    if len(set(out.values())) != 3:
        raise ValueError('the three roles need three different channels')
    return out


def czi_set(path, roles, dataset, root):
    """Image-set dict (as pipeline.find_sets returns) for one CZI file with role -> channel index."""
    s = {k: channel_ref(path, roles[k]) for k in ('red', 'green', 'blue')}
    return dict(s, dataset=dataset, name=os.path.splitext(os.path.basename(path))[0], root=root,
                czi=path, roles=dict(roles))


# ---- CZI file format -------------------------------------------------------------------------------

def _segment(fh, pos, want=None):
    fh.seek(pos)
    sid, _alloc, used = struct.unpack('<16sqq', fh.read(32))
    sid = sid.rstrip(b'\0').decode('ascii', 'replace')
    if want and sid != want:
        raise ValueError(f'not a CZI file or damaged ({sid!r} at {pos}, expected {want})')
    return sid, used


def _dir_entry(buf, o):
    """DirectoryEntryDV at offset o -> (entry dict, offset after it)."""
    if buf[o:o + 2] != b'DV':
        raise ValueError('unsupported CZI directory entry')
    ptype, fpos, _part, comp, pyr = struct.unpack_from('<iqiiB', buf, o + 2)
    ndim, = struct.unpack_from('<i', buf, o + 28)
    dims, p = {}, o + 32
    for _ in range(ndim):
        name, start, size, _coord, stored = struct.unpack_from('<4siifi', buf, p)
        dims[name.rstrip(b'\0').decode()] = (start, size, stored)
        p += 20
    return dict(ptype=ptype, pos=fpos, comp=comp, pyramid=pyr, dims=dims), p


def _header(fh):
    _segment(fh, 0, 'ZISRAWFILE')
    fh.seek(32)
    raw = fh.read(80)
    dir_pos, meta_pos = struct.unpack_from('<qq', raw, 4 * 4 + 32 + 4)
    return dir_pos, meta_pos


def _directory(fh, dir_pos):
    _segment(fh, dir_pos, 'ZISRAWDIRECTORY')
    n, = struct.unpack('<i', fh.read(4))
    fh.read(124)
    buf = fh.read(n * (32 + 20 * 12) + 4096)
    entries, o = [], 0
    for _ in range(n):
        if o + 32 > len(buf):  # many dimensions: read more
            buf += fh.read(1 << 20)
        e, o = _dir_entry(buf, o)
        entries.append(e)
    return entries


@lru_cache(maxsize=64)
def _czi_info(path, _mtime, _size):
    with open(path, 'rb') as fh:
        dir_pos, meta_pos = _header(fh)
        entries = _directory(fh, dir_pos)
        xml = ''
        if meta_pos:
            _segment(fh, meta_pos, 'ZISRAWMETADATA')
            xml_size, = struct.unpack('<i', fh.read(4))
            fh.seek(meta_pos + 32 + 256)
            xml = fh.read(xml_size).decode('utf-8', 'replace')
    root = ET.fromstring(xml) if xml else ET.Element('none')

    def txt(el, tag):
        x = el.find(tag)
        return x.text.strip() if x is not None and x.text else ''

    n_c = max((e['dims'].get('C', (0, 1, 1))[0] + 1 for e in entries), default=1)
    chans = []
    dims_ch = root.findall('.//Information/Image/Dimensions/Channels/Channel')
    for i in range(n_c):
        el = dims_ch[i] if i < len(dims_ch) else ET.Element('Channel')
        em = txt(el, 'EmissionWavelength')
        try:
            em = float(em)
        except ValueError:
            em = 0.0
        chans.append(dict(name=el.get('Name') or f'Channel {i}', fluor=txt(el, 'Fluor'), emission=em))
    px = None
    for d in root.findall('.//Scaling/Items/Distance'):
        if d.get('Id') == 'X':
            try:
                v = float(txt(d, 'Value')) * 1e6  # metres -> um
                px = v if 0.005 <= v <= 10 else None
            except ValueError:
                pass
    bits = txt(root, './/Information/Image/ComponentBitCount')
    full = [e for e in entries if not e['pyramid']]
    ys = [(e['dims']['Y'][0], e['dims']['Y'][0] + e['dims']['Y'][1]) for e in full if 'Y' in e['dims']]
    xs = [(e['dims']['X'][0], e['dims']['X'][0] + e['dims']['X'][1]) for e in full if 'X' in e['dims']]

    def count(k):
        return len({e['dims'][k][0] for e in full if k in e['dims']}) or 1

    return dict(channels=chans, px_um=px, bits=int(bits) if bits.isdigit() else None,
                size=(max(b for _, b in ys) - min(a for a, _ in ys), max(b for _, b in xs) - min(a for a, _ in xs)),
                n_z=count('Z'), n_t=count('T'), n_scenes=count('S'), entries=full, ptypes={e['ptype'] for e in full})


def _decode(raw, comp, dtype, shape):
    if comp == 0:
        return np.frombuffer(raw, dtype).reshape(shape)
    try:
        import imagecodecs
    except ImportError:
        raise ValueError('this CZI file is compressed; reading it needs the imagecodecs package '
                         '(pip install imagecodecs)') from None
    if comp == 4:  # JPEG XR
        return np.asarray(imagecodecs.jpegxr_decode(raw)).reshape(shape)
    if comp == 1:
        return np.asarray(imagecodecs.jpeg8_decode(raw)).reshape(shape)
    if comp == 2:
        return np.frombuffer(imagecodecs.lzw_decode(raw), dtype).reshape(shape)
    if comp == 5:  # zstd, no header
        return np.frombuffer(imagecodecs.zstd_decode(raw), dtype).reshape(shape)
    if comp == 6:  # zstd with header; optional hi/lo byte packing of 16-bit data
        n = raw[0]
        hilo = n >= 3 and raw[1] == 1 and raw[2] & 1
        out = np.frombuffer(imagecodecs.zstd_decode(bytes(raw[n:])), np.uint8)
        if hilo:
            half = out.size // 2
            out = np.stack([out[:half], out[half:]], -1).ravel()
        return out.view(dtype).reshape(shape)
    raise ValueError(f'unsupported CZI compression {comp}')


@lru_cache(maxsize=8)
def _czi_channel(path, ch, _mtime, _size):
    info = _czi_info(path, _mtime, _size)
    full = info['entries']
    if not full:
        raise ValueError(f'{os.path.basename(path)}: no image data')
    first = {k: min(e['dims'][k][0] for e in full if k in e['dims']) for k in ('T', 'S', 'B', 'H', 'V', 'I', 'R')
             if any(k in e['dims'] for e in full)}
    sel = [e for e in full if e['dims'].get('C', (0,))[0] == ch
           and all(e['dims'].get(k, (v,))[0] == v for k, v in first.items())]
    if not sel:
        raise ValueError(f'{os.path.basename(path)} has no channel {ch} (it has {len(info["channels"])})')
    y0 = min(e['dims']['Y'][0] for e in full); x0 = min(e['dims']['X'][0] for e in full)
    h, w = info['size']
    planes = {}
    with open(path, 'rb') as fh:
        for e in sel:
            if e['ptype'] not in PIXEL_TYPES:
                raise ValueError(f'unsupported CZI pixel type {e["ptype"]}')
            dtype, spp = PIXEL_TYPES[e['ptype']]
            _segment(fh, e['pos'], 'ZISRAWSUBBLOCK')
            meta_size, _att, data_size = struct.unpack('<iiq', fh.read(16))
            ent, _ = _dir_entry(fh.read(4096), 0)
            hdr = max(256, 16 + 32 + 20 * len(ent['dims']))
            fh.seek(e['pos'] + 32 + hdr + meta_size)
            sy, sx = e['dims']['Y'][2], e['dims']['X'][2]
            tile = _decode(fh.read(data_size), e['comp'], dtype, (sy, sx, spp) if spp > 1 else (sy, sx))
            if tile.ndim == 3:  # BGR(A) colour: use the luminance-like mean of the colour samples
                tile = tile[..., :3].mean(-1)
            z = e['dims'].get('Z', (0,))[0]
            plane = planes.get(z)
            if plane is None:
                plane = planes[z] = np.zeros((h, w), tile.dtype)
            ys, xs = e['dims']['Y'][0] - y0, e['dims']['X'][0] - x0
            plane[ys:ys + sy, xs:xs + sx] = tile[:h - ys, :w - xs]
    img = np.max(np.stack(list(planes.values())), 0) if len(planes) > 1 else next(iter(planes.values()))
    if img.dtype != np.uint8:
        top = float(2 ** info['bits'] - 1) if info['bits'] else (65535.0 if img.dtype == np.uint16 else
                                                                  max(float(img.max()), 1.0))
        img = np.clip(np.round(img.astype(float) * (255.0 / top)), 0, 255).astype(np.uint8)
    return img
