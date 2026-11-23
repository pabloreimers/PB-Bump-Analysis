"""LPsP morphology in the whole fly brain, from the public FlyWire connectome (v783).

Draws the two LPsP neurons (one per side) over the FlyWire brain-surface outline, with
only the PB and the LALs shaded -- the neuropils where LPsP actually arborizes. An
in-volume test of each LPsP mesh against the FlyWire neuropil meshes gives:
    PB ~35-39%, own-side LAL ~13%, WED / SPS / IPS 2-10% each,
    ~30% outside every neuropil (soma + cell-body fibre + connecting neurites)
so "PB + LAL" is the honest simplification; the lateral arbor spills ventrally from the
LAL into the posterior slope, which is where the name comes from.

Why FlyWire and not hemibrain: in hemibrain v1.2.1 both LPsP neurites leave the imaged
volume ventrally, so the somata and posterior-slope arbors are missing. FlyWire is a
complete brain. (Hemibrain is still the better source for PB glomerulus-level synapse
counts; see ~/Documents/GitHub/hemibrain.)

Neurotransmitter (FlyWire public annotations, Supplemental_file1_neuron_annotations.tsv):
known_nt = dopamine + glutamate; ACh-, GABA-, 5-HT-, TA-, OA-negative (Alekseyenko et al.
2013 immuno; Reimers et al. intersection; Wolff et al. 2024 EASI-FISH). The image-based
top_nt prediction is low-confidence and disagrees between the two cells (ACh 0.48 R,
glutamate 0.45 L) -- don't quote it.

Coordinates: FlyWire nm frame, checked from neuropil centroids:
    x: fly's LEFT (low) -> RIGHT (high)     (ME_R at high x)
    y: DORSAL (low)     -> VENTRAL (high)   (GNG at high y)
    z: ANTERIOR (low)   -> POSTERIOR (high) (EB < PB)
View conventions: frontal = seen from the front (fly's right on page left); posterior =
from behind (fly's right on page right, as in our 2p imaging); dorsal = from above with
anterior up (fly's right on page right); lateral = from the fly's right side.

Rendering: neuron meshes are dense (~1M vertices), so they're drawn by rasterising the
projected vertices; neuropil/brain meshes are simplified, so their triangles are
rasterised and outlined. Output is PNG (raster) -- the figure scripts import it as an
image; scale bars are in the image, so they survive resizing in Illustrator.

Data: everything is public, no FlyWire/CAVE login needed. Downloaded once into
CACHE (data/flywire_lpsp/, gitignored):
    LPsP meshes     precomputed://gs://flywire_v141_m783 (cloud-volume, https)
    neuropil meshes fafbseg.flywire.get_neuropil_volumes
    brain surface   flybrains' bundled FLYWIRE.ply

Python env (outside MATLAB). h5py doesn't build on Python 3.14 without pkg-config, so
flybrains is installed --no-deps (only its bundled meshes are used):
    python3 -m venv ~/venvs/flywire
    ~/venvs/flywire/bin/pip install navis fafbseg cloud-volume scikit-image trimesh
    ~/venvs/flywire/bin/pip install --no-deps flybrains
Run:
    ~/venvs/flywire/bin/python ugly_figures/scripts/fig_lpsp_morphology_flywire.py [out_dir]
"""
import os, glob, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from scipy import ndimage
from skimage.draw import polygon as skpoly
from skimage import measure

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
CACHE = os.path.join(REPO, 'data', 'flywire_lpsp')
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, 'ugly_figures', 'exports', 'V2', 'fig1_anatomy')
os.makedirs(OUT, exist_ok=True)

NEURONS = {720575940617672226: ('LPsP R', '#c0007a'), 720575940624452902: ('LPsP L', '#0a8a8a')}
SHADED = ['PB', 'LAL_L', 'LAL_R']

# horizontal, vertical page coordinates (both increase right/up), axis-key labels
VIEWS = {
    'frontal':   (lambda p: -p[:, 0], lambda p: -p[:, 1], ('R', 'L'), ('V', 'D')),
    'posterior': (lambda p: p[:, 0],  lambda p: -p[:, 1], ('L', 'R'), ('V', 'D')),
    'dorsal':    (lambda p: p[:, 0],  lambda p: -p[:, 2], ('L', 'R'), ('P', 'A')),
    'lateral':   (lambda p: -p[:, 2], lambda p: -p[:, 1], ('P', 'A'), ('V', 'D')),
}


# ----------------------------------------------------------------------------- data
def fetch_cache():
    """Download (once) the LPsP meshes, neuropil meshes and brain surface into CACHE."""
    os.makedirs(os.path.join(CACHE, 'neuropils'), exist_ok=True)
    missing = [i for i in NEURONS if not os.path.exists(os.path.join(CACHE, f'lpsp_{i}.npz'))]
    if missing:
        import cloudvolume as cv
        vol = cv.CloudVolume('precomputed://gs://flywire_v141_m783', use_https=True, progress=False)
        meshes = vol.mesh.get(missing)
        for i in missing:
            np.savez_compressed(os.path.join(CACHE, f'lpsp_{i}.npz'), v=meshes[i].vertices, f=meshes[i].faces)
    if not glob.glob(os.path.join(CACHE, 'neuropils', '*.npz')):
        from fafbseg import flywire
        names = list(flywire.get_neuropil_volumes(None))
        for n, v in zip(names, flywire.get_neuropil_volumes(names)):
            np.savez_compressed(os.path.join(CACHE, 'neuropils', f'{n}.npz'),
                                v=np.asarray(v.vertices), f=np.asarray(v.faces))
    if not os.path.exists(os.path.join(CACHE, 'brain_FLYWIRE.npz')):
        import trimesh, flybrains
        m = trimesh.load(os.path.join(os.path.dirname(flybrains.__file__), 'meshes', 'FLYWIRE.ply'))
        np.savez_compressed(os.path.join(CACHE, 'brain_FLYWIRE.npz'), v=np.asarray(m.vertices), f=np.asarray(m.faces))


fetch_cache()
neuropils = {os.path.basename(f)[:-4]: np.load(f) for f in glob.glob(os.path.join(CACHE, 'neuropils', '*.npz'))}
neurons = {i: np.load(os.path.join(CACHE, f'lpsp_{i}.npz'))['v'] for i in NEURONS}
BRAIN = np.load(os.path.join(CACHE, 'brain_FLYWIRE.npz'))


# ----------------------------------------------------------------------------- drawing
def raster_faces(v2, f, ext, res):
    x0, x1, y0, y1 = ext
    img = np.zeros((int((y1 - y0) / res) + 1, int((x1 - x0) / res) + 1), bool)
    px = (v2[:, 0] - x0) / res; py = (v2[:, 1] - y0) / res
    for tri in f:
        rr, cc = skpoly(py[tri], px[tri], img.shape)
        img[rr, cc] = True
    return ndimage.binary_closing(img, iterations=2)


def contours(img, ext, res):
    return [np.c_[c[:, 1] * res + ext[0], c[:, 0] * res + ext[2]]
            for c in measure.find_contours(img.astype(float), 0.5) if len(c) > 10]


def neuron_img(v2, ext, res, grow=1):
    x0, x1, y0, y1 = ext
    img = np.zeros((int((y1 - y0) / res) + 1, int((x1 - x0) / res) + 1), bool)
    ix = ((v2[:, 0] - x0) / res).astype(int); iy = ((v2[:, 1] - y0) / res).astype(int)
    ok = (ix >= 0) & (ix < img.shape[1]) & (iy >= 0) & (iy < img.shape[0])
    img[iy[ok], ix[ok]] = True
    return ndimage.binary_dilation(img, iterations=grow) if grow else img


def draw(view, ax, which=None, ext=None, res=None, labels=True, grow=1):
    H, V, hl, vl = VIEWS[view]
    which = list(NEURONS) if which is None else which
    if ext is None:
        allv = BRAIN['v']; pad = 8000
        ext = (H(allv).min() - pad, H(allv).max() + pad, V(allv).min() - pad, V(allv).max() + pad)
    res = res or (ext[1] - ext[0]) / 2400
    # whole-brain outline: the FlyWire brain surface, NOT a union of neuropils (the union
    # is lumpy where the optic-lobe neuropils don't tile)
    v2 = np.c_[H(BRAIN['v']), V(BRAIN['v'])]
    brain = ndimage.binary_fill_holes(raster_faces(v2, BRAIN['f'], ext, res))
    for c in contours(brain, ext, res):
        ax.fill(c[:, 0], c[:, 1], color='#f2f2f2', lw=0, zorder=0)
        ax.plot(c[:, 0], c[:, 1], color='#b5b5b5', lw=0.6, zorder=1)
    for r in SHADED:
        d = neuropils[r]; v2 = np.c_[H(d['v']), V(d['v'])]
        for c in contours(raster_faces(v2, d['f'], ext, res), ext, res):
            ax.fill(c[:, 0], c[:, 1], color='#9a9a9a', alpha=0.3, lw=0, zorder=2)
            ax.plot(c[:, 0], c[:, 1], color='#8a8a8a', lw=0.5, zorder=2)
        if labels and r != 'LAL_R':  # one LAL label is enough; placed just above the neuropil
            ax.text(v2[:, 0].mean(), v2[:, 1].max(), r.split('_')[0], fontsize=7, color='#444444',
                    ha='center', va='bottom', zorder=7, fontweight='bold')
    rgba = np.zeros(brain.shape + (4,))
    for i in which:  # later entries drawn on top
        v = neurons[i]
        m = neuron_img(np.c_[H(v), V(v)], ext, res, grow)
        rgba[m, :3] = to_rgb(NEURONS[i][1]); rgba[m, 3] = 0.95
    ax.imshow(rgba, extent=(ext[0] - res / 2, ext[0] + (brain.shape[1] - .5) * res,
                            ext[2] - res / 2, ext[2] + (brain.shape[0] - .5) * res),
              origin='lower', interpolation='nearest', zorder=5)
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    ax.set_aspect('equal'); ax.axis('off')
    for i in which:
        ax.plot([], [], color=NEURONS[i][1], lw=2, label=NEURONS[i][0])
    ax.legend(loc='upper right', frameon=False, fontsize=7)
    # scale bar + orientation key (baked into the PNG, so they survive resizing)
    w = ext[1] - ext[0]; h = ext[3] - ext[2]
    sb = 100000 if w > 400000 else (20000 if w > 80000 else 10000)
    xs = ext[1] - 0.04 * w
    ax.plot([xs - sb, xs], [ext[2] + 0.04 * h] * 2, 'k', lw=2, zorder=8)
    ax.text(xs - sb / 2, ext[2] + 0.06 * h, f'{sb / 1000:.0f} µm', ha='center', va='bottom', fontsize=7)
    ox, oy, L = ext[0] + 0.07 * w, ext[2] + 0.1 * h, 0.035 * w
    for a, b in [((ox - L, oy), (ox + L, oy)), ((ox, oy - L), (ox, oy + L))]:
        ax.annotate('', b, a, arrowprops=dict(arrowstyle='<->', lw=0.8))
    ax.text(ox - L * 1.4, oy, hl[0], ha='right', va='center', fontsize=7)
    ax.text(ox + L * 1.4, oy, hl[1], ha='left', va='center', fontsize=7)
    ax.text(ox, oy - L * 1.4, vl[0], ha='center', va='top', fontsize=7)
    ax.text(ox, oy + L * 1.4, vl[1], ha='center', va='bottom', fontsize=7)


# ----------------------------------------------------------------------------- figures
if __name__ == '__main__':
    # whole brain, every view (titled, for browsing)
    for v in VIEWS:
        fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
        draw(v, ax)
        ax.set_title(f'LPsP in the whole brain (FlyWire v783) — {v} view', fontsize=9)
        fig.savefig(os.path.join(OUT, f'fw_lpsp_brain_{v}.png'), bbox_inches='tight', facecolor='white')
        plt.close(fig)

    # untitled posterior whole-brain panel: the one Figure_1_V2.m imports
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    draw('posterior', ax)
    fig.savefig(os.path.join(OUT, 'fig1_lpsp_morphology_posterior.png'), bbox_inches='tight',
                facecolor='white', pad_inches=0.02)
    plt.close(fig)

    # central-brain crops: each neuron alone, then both
    for v in ['posterior', 'dorsal', 'lateral']:
        H, V, *_ = VIEWS[v]
        allp = np.vstack(list(neurons.values()) + [neuropils[r]['v'] for r in SHADED])
        m = 10000
        ext = (H(allp).min() - m, H(allp).max() + m, V(allp).min() - m, V(allp).max() + m)
        fig, axs = plt.subplots(1, 3, figsize=(13, 4.6), dpi=300)
        for ax, sel in zip(axs, [[720575940617672226], [720575940624452902], list(NEURONS)]):
            draw(v, ax, which=sel, ext=ext, res=(ext[1] - ext[0]) / 1400, grow=0)
        fig.suptitle(f'LPsP, {v} view (FlyWire v783)', fontsize=10)
        fig.savefig(os.path.join(OUT, f'fw_lpsp_closeup_{v}.png'), bbox_inches='tight', facecolor='white')
        plt.close(fig)
    written = sorted(f for f in os.listdir(OUT) if f.startswith(('fw_', 'fig1_')))
    print('wrote', written)

    # shared copy in Dropbox, next to the V2 figure PDFs (skipped if Dropbox isn't on this machine)
    import shutil
    lpsp_dropbox = os.path.expanduser('~/Dropbox-HarvardUniversity/Pablo Reimers/LPsP')
    if os.path.isdir(lpsp_dropbox):
        dst = os.path.join(lpsp_dropbox, 'V2', 'fig1_anatomy')
        os.makedirs(dst, exist_ok=True)
        for f in written:
            shutil.copy2(os.path.join(OUT, f), dst)
        print('copied to', dst)
