"""
Write a CXI v1.6 file for publication from a merged hits file plus the per-event
metadata collected on maxwell by get_event_metadata.py.

    python make_publication_cxi.py <merged.cxi> <event_metadata.h5> <output.cxi>
        --mask Ery/recon_3D_nosym/mask.h5 --edge-mask Ery/edge_mask.h5 [--geom-run 600]

Inputs
    merged.cxi          frames, hit scores, per-run backgrounds and manual selection
                        (merge_cxi_files.py + add_background_cxi.py + local edits)
    event_metadata.h5   run, timestamp, pulseId, per-pulse XGM/LITFRM/hit finder values,
                        per-train control values, recomputed background weighting
                        (get_event_metadata.py, run on maxwell)
    geom/r{run}.geom    per-run AGIPD geometry (reference geometry + quadrant motors)

Layout (CXIDB spec, NeXus-compatible). Per-event arrays are indexed by
experiment_identifier; per-pixel arrays by module_identifier:y:x, where a
"module" is one of the 128 AGIPD tiles (see Detector tiles below). N is the
number of events, Nb the number of runs (one background per run).

    /
      cxi_version = 160                                            int64 (160 = v1.6)
      entry_1/                                                     NXentry
        experiment_identifier   (N,)            string             p007927_r{run:04d}_t{trainId}_c{cellId}
        title                   scalar          string             "p007927 {sample} hits, runs {first}-{last}"
        experiment_description  scalar          string             --description
        program_name            scalar          string             "xfel7927/offline/make_publication_cxi.py"
        start_time              scalar          string             ISO 8601, start of the first run
        end_time                scalar          string             ISO 8601, last train in the file
        sample_1/                                                  NXsample
          name                  scalar          string             "Erythrocruorin"
        instrument_1/                                              NXinstrument
          name                  scalar          string             "SPB"
          source_1/                                                NXsource
            name                scalar          string             "European XFEL SASE1"
            energy              (N,)            float32  J         photon energy, undulator (per train)
            pulse_energy        (N,)            float32  J         XGM SPB_XTD9 energy of this frame's pulse (via LITFRM)
            pulse_energy_sigma  (N,)            float32  J         XGM uncertainty of pulse_energy
            experiment_identifier -> /entry_1/experiment_identifier
          xgm_1/                                                   NXcollection, SPB_XTD9_XGM (downstream)
            description         scalar          string
            intensity           (N,)            float32  J         pulse energy of this frame's pulse
            intensity_sigma     (N,)            float32  J
            x, y                (N,)            float32  m         beam position at the XGM
            x_sigma, y_sigma    (N,)            float32  m
            experiment_identifier -> /entry_1/experiment_identifier
          xgm_2/                                                   NXcollection, SA1_XTD2_XGM (upstream), as xgm_1
          attenuator_1/                                            NXattenuator, SA1_XTD2_ATT
            type                scalar          string             "SA1_XTD2_ATT"
            attenuator_transmission (N,)        float32  dimensionless  (per train)
            experiment_identifier -> /entry_1/experiment_identifier
          attenuator_2/                                            NXattenuator, SPB_XTD9_ATT, as attenuator_1
          electrospray/                                            NXcollection, injector (all per train)
            CO2_capillary       (N,)            float32  ?         flow controller reading (measureCapacity)
            CO2_chamber         (N,)            float32  ?            "
            DP2                 (N,)            float32  ?            "
            He_chamber          (N,)            float32  ?            "
            N2_capillary        (N,)            float32  ?            "
            N2_chamber          (N,)            float32  ?            "
            liquidjet_He        (N,)            float32  ?            "
            voltage             (N,)            float32  V         electrospray high voltage
            current             (N,)            float32  A         electrospray current
            injector_x, injector_y, injector_z (N,) float32 m       injector motor positions
            chamber_pressure    (N,)            float32  ?         SPB_IRU_VAC/GAUGE/GAUGE_FR_6 (likely mbar)
            experiment_identifier -> /entry_1/experiment_identifier
          detector_1/                                              NXdetector
            description         scalar          string             "AGIPD 1M"
            distance            scalar          float64  m         sample to detector distance (DET_DIST)
            x_pixel_size        scalar          float64  m
            y_pixel_size        scalar          float64  m
            module_identifier   (128,)          string             "AGIPD{mm}T{t}", tile t of module mm
            corner_position     (128, 3)        float32  m         corner of pixel (0, 0) of each tile
            basis_vectors       (128, 2, 3)     float32  m         pixel steps along y (slow) and x (fast)
            xyz_map             (3, 128, 64, 128) float32 m        pixel centre positions, same geometry
            quadrant            (4,)            string             "Q1".."Q4"
            quadrant_correction (4, 3)          float32  m         alternative geometry, see Geometry below
            mask                (128, 64, 128)  uint32             CXI mask bits, see Mask below
            data                (N, 128, 64, 128) uint8  counts    photons per pixel, NOT masked
            powder              (128, 64, 128)  float32  counts    mean of data over all events
            start_time          (N,)            string             train timestamp, ISO 8601 (ns)
            run                 (N,)            uint16             run number
            trainId             (N,)            uint64
            cellId              (N,)            uint16             AGIPD memory cell
            pulseId             (N,)            int64              pulse index in the train
            vds_index           (N,)            int64              frame index in the per-run VDS file
            data_white          (Nb, 128, 64, 128) float32 counts  per-run background, see Background below
            background_run      (Nb,)           uint16             run of each data_white
            background_index    (N,)            uint32             index into data_white for each event
            background_weighting (N,)           float32            scale of data_white for each event
            score/                                                 per-event scores
              hit_sigma         (N,)            float32            (photons in hit_finding_mask - train median)
                                                                   / (1.4826 * (median - 25th percentile));
                                                                   hits were selected with hit_sigma > 4
              photon_counts     (N,)            float32            photons in the frame (per-cell mask applied)
              lit_pixels        (N,)            float32            pixels with >= 1 photon (per-cell mask applied)
              background_counts (N,)            float32            sum_pix(background_weighting * data_white[background_index])
              facility_hitscore (N,)            float32            SPI_HITFINDER hitscore
              facility_hit_flag (N,)            int8               SPI_HITFINDER hitFlag (1 = hit)
              facility_miss_flag (N,)           int8               SPI_HITFINDER missFlag (1 = miss)
              facility_threshold_mu (N,)        float32            SPI_HITFINDER threshold mu (per train)
              facility_threshold_sigma (N,)     float32            SPI_HITFINDER threshold sigma (per train)
              facility_lit_pixels (N,)          float32            correction pipeline litPixels, summed over modules
              facility_total_intensity (N,)     float32            correction pipeline totalIntensity, summed over modules
              facility_unmasked_pixels (N,)     float32            correction pipeline unmaskedPixels, summed over modules
              experiment_identifier -> /entry_1/experiment_identifier
            tag                 (3,)            string             "is_crap", "is_good", "is_strong_hit"
            tags                (3, N)          int8               manual selection, 1 = tagged
            motors/                                                NXcollection, AGIPD motors (per train)
              Q1M1 .. Q4M2      (N,)            float32  m         quadrant motor positions
              z_stepper         (N,)            float32  m         detector z stage position
              experiment_identifier -> /entry_1/experiment_identifier
            note_1/                                                NXnote, geometry provenance
              file_name         scalar          string             reference CrystFEL geometry file name
              description       scalar          string             how the geometry was derived
              data              scalar          string             contents of the reference geometry file
            experiment_identifier -> /entry_1/experiment_identifier
        data_1/                                                    NXdata
          data -> /entry_1/instrument_1/detector_1/data
          experiment_identifier -> /entry_1/experiment_identifier
        process_1/                                                 NXprocess
          program               scalar          string             "make_publication_cxi.py"
          version               scalar          string             git describe of xfel7927
          date                  scalar          string             ISO 8601
          command               scalar          string             command line

Axes attributes used as dimension scales:
    every (N,) dataset                            axes = "experiment_identifier"
    data                                          axes = "experiment_identifier:module_identifier:y:x"
    mask, powder                                  axes = "module_identifier:y:x"
    xyz_map                                       axes = "coordinate:module_identifier:y:x"
    corner_position                               axes = "module_identifier:coordinate"
    basis_vectors                                 axes = "module_identifier:dimension:coordinate"
    quadrant_correction                           axes = "quadrant:coordinate"
    data_white                                    axes = "background_run:module_identifier:y:x"
    tags                                          axes = "tag:experiment_identifier"
Each named axis is a dataset (or link) in the same group, except the implicit
axes y, x, coordinate and dimension. Per-event datasets from the metadata file
also carry `source` and `key` attributes (Karabo device and key); where the
facility data records no unit, `units` is omitted and `units_note` says so
(marked "?" above).

Detector tiles
    Each AGIPD module (512 x 128 pixels) is 8 tiles of 64 x 128 pixels separated
    by 2 pixel gaps, so a single corner/basis per module would misplace pixels by
    up to 14 pixels. module_identifier "AGIPD{mm}T{t}" (index 8 * mm + t) is rows
    64t..64t+63 of module mm. Any per-pixel array reshapes (C order) to the
    facility layout:  data[d].reshape(16, 512, 128).
    Pixel centre of (tile m, row i, column j):
        corner_position[m] + (i + 0.5) * basis_vectors[m, 0] + (j + 0.5) * basis_vectors[m, 1]
    Quadrant q (0-3) is modules 4q..4q+3, i.e. tiles 32q..32q+31.

Geometry
    The facility reference geometry (note_1, agipd_p008316_r0024_v04.geom) moved
    to the recorded quadrant motor positions with extra_geom.motors.AGIPD_1MMotors
    (geom/r{run}.geom; the motors are constant for all runs in this file to 2 um),
    with z = distance and shifted by BEAM_CENTRE_SHIFT = (-500, -500) um in x, y.
    The shift is the mean beam centre correction found from the data, which
    agrees with an independent crystallographic refinement (O. Yefanov, p8216).
    quadrant_correction is an alternative, per-quadrant refinement fitted to the
    data and used in the 3D reconstructions; adding quadrant_correction[q] to the
    positions of the tiles in quadrant q reproduces that geometry. It is not
    supported by the independent refinement and is given for reproducibility.

Mask (bits set = pixel excluded; 0 = good)
    0x00000080  CXI_PIXEL_IS_BAD   excluded by the analysis mask (--mask)
    0x00010000  user               bad in every memory cell of every run (merged per-cell masks)
    0x00020000  user               ASIC edge pixel (--edge-mask)
    The definitions are also stored in mask.attrs["bit_definitions"].

Background
    data_white[b] is the per-pixel mean of the non-hit frames of run
    background_run[b] (scratch/powder/r{run}_powder_is_hit_False_per_pixel.h5),
    with zero pixels and pixels deviating by more than 20% from the
    (polarisation and solid angle corrected) radial profile replaced by that
    profile (Ery/fill_background_gaps.py). The estimated background of event d is
        B_d = background_weighting[d] * data_white[background_index[d]]
    with
        background_weighting[d] = a_t * e_d / <e>
        a_t  = mean photon count of the non-hit frames in the train of event d
               / total photon count of the run's non-hit mean (before filling)
               (-1 if the train has no non-hit frames)
        e_d  = XGM pulse energy of event d (source_1/pulse_energy)
        <e>  = mean of e over all frames of the run with e > 1 mJ
               (e_d / <e> is set to 1 for frames with e <= 1 mJ)
    score/background_counts[d] = sum over pixels of B_d.

Powder
    powder is the mean of data over all N events (unmasked, background included).
"""
import argparse
import datetime
import os
import subprocess
import sys

import numpy as np
import h5py
import extra_geom
from tqdm import tqdm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PROPOSAL = 7927
DET_DIST = 715e-3
REFERENCE_GEOM = 'geom/agipd_p008316_r0024_v04.geom'

# mean beam centre correction: pixel positions are shifted by this relative to the geometry file
BEAM_CENTRE_SHIFT = np.array([-500e-6, -500e-6, 0.])

# per-quadrant shifts (pixels) relative to the geometry file used in the 3D
# reconstructions (Ery/recon_3D_nosym/config.py, "round 4"), stored relative to
# BEAM_CENTRE_SHIFT as quadrant_correction
RECON_QUADRANT_SHIFTS_PX = np.array([[-2, -4], [0, -2], [-4, -2], [-4, -2]], dtype=float)

# CXI mask bits (cxi.h) and user bits
CXI_PIXEL_IS_BAD = 0x00000080
BIT_NOT_GOOD_PIXELS = 0x00010000
BIT_ASIC_EDGE = 0x00020000

GZ = dict(compression='gzip', compression_opts=1, shuffle=True)

NMODULES, NTILES = 16, 8
TILE_SHAPE = (64, 128)
PANEL_SHAPE = (NMODULES * NTILES,) + TILE_SHAPE

DET = 'entry_1/instrument_1/detector_1'


def to_tiles(a):
    """(..., 16, 512, 128) -> (..., 128, 64, 128)"""
    return a.reshape(a.shape[:-3] + PANEL_SHAPE)


def git_version():
    try:
        return subprocess.check_output(['git', '-C', ROOT, 'describe', '--always', '--dirty'], text=True).strip()
    except Exception:
        return 'unknown'


def link_expid(group):
    if 'experiment_identifier' not in group:
        group['experiment_identifier'] = h5py.SoftLink('/entry_1/experiment_identifier')


def per_event(group, name, data, units=None, **attrs):
    """write a per-event dataset with axes and the experiment_identifier link"""
    ds = group.create_dataset(name, data=data, **GZ)
    ds.attrs['axes'] = 'experiment_identifier'
    if units:
        ds.attrs['units'] = units
    for k, v in attrs.items():
        ds.attrs[k] = v
    link_expid(group)
    return ds


def copy_meta(group, name, meta, key, **attrs):
    """copy a metadata column, keeping its units/source/key attributes"""
    if key not in meta:
        print(f'warning: {key} not in metadata file, {group.name}/{name} skipped')
        return
    ds = meta[key]
    units = ds.attrs.get('units', '')
    known = units and not units.startswith('unknown')
    out = per_event(group, name, ds[()], units if known else None, **attrs)
    if units.startswith('unknown'):
        out.attrs['units_note'] = units + ' (not recorded in the facility data)'
    for a in ('source', 'key'):
        if a in ds.attrs:
            out.attrs[a] = ds.attrs[a]


def load_geometry(runs, geom_run):
    """corner_position, basis_vectors, xyz_map per tile for the (single) geometry of all runs"""
    fnam = os.path.join(ROOT, f'geom/r{geom_run:04d}.geom')
    geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(fnam)
    pos = geom.get_pixel_positions()

    # check that every run with a geometry file has the same geometry
    for r in runs:
        f = os.path.join(ROOT, f'geom/r{r:04d}.geom')
        if not os.path.exists(f):
            print(f'warning: no geometry file for run {r}, assuming r{geom_run:04d}')
            continue
        p = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(f).get_pixel_positions()
        d = np.abs(p - pos).max() / geom.pixel_size
        if d > 0.1:
            raise ValueError(f'geometry of run {r} differs from run {geom_run} by {d:.2f} pixels')

    pos[..., 2] = DET_DIST
    pos += BEAM_CENTRE_SHIFT
    pos = pos.reshape(PANEL_SHAPE + (3,))

    basis = np.empty((pos.shape[0], 2, 3), dtype=np.float32)
    basis[:, 0] = pos[:, 1, 0] - pos[:, 0, 0]     # slow scan (y)
    basis[:, 1] = pos[:, 0, 1] - pos[:, 0, 0]     # fast scan (x)
    corner = (pos[:, 0, 0] - 0.5 * basis[:, 0] - 0.5 * basis[:, 1]).astype(np.float32)

    # every pixel centre must be reproduced by the per-tile corner and basis vectors
    i, j = np.mgrid[:TILE_SHAPE[0], :TILE_SHAPE[1]] + 0.5
    model = corner[:, None, None] + i[..., None] * basis[:, None, None, 0] + j[..., None] * basis[:, None, None, 1]
    err = np.abs(model - pos).max() / geom.pixel_size
    assert err < 1e-3, f'per-tile geometry is not exact ({err} pixels)'

    xyz = np.transpose(pos, (3, 0, 1, 2)).astype(np.float32)

    quadrant_correction = np.zeros((4, 3), dtype=np.float32)
    quadrant_correction[:, :2] = RECON_QUADRANT_SHIFTS_PX * geom.pixel_size - BEAM_CENTRE_SHIFT[:2]
    return corner, basis, xyz, quadrant_correction, float(geom.pixel_size), fnam


def build_mask(good_pixels, mask_file, edge_file):
    mask = np.zeros(good_pixels.shape, dtype=np.uint32)
    defs = []
    if mask_file:
        with h5py.File(mask_file) as f:
            good = f['data'][()].astype(bool)
        mask[~good] |= CXI_PIXEL_IS_BAD
        name = os.path.join(*os.path.normpath(mask_file).split(os.sep)[-2:])
        defs.append(f'0x{CXI_PIXEL_IS_BAD:08x} CXI_PIXEL_IS_BAD: excluded by the analysis mask ({name})')
    mask[~good_pixels] |= BIT_NOT_GOOD_PIXELS
    defs.append(f'0x{BIT_NOT_GOOD_PIXELS:08x} user: bad in every memory cell of every run (merged per-cell masks)')
    if edge_file:
        with h5py.File(edge_file) as f:
            edge_good = f['data'][()].astype(bool)
        mask[~edge_good] |= BIT_ASIC_EDGE
        defs.append(f'0x{BIT_ASIC_EDGE:08x} user: ASIC edge pixel')
    return to_tiles(mask), defs


def main():
    parser = argparse.ArgumentParser(description='Write a publication CXI file from a merged hits file and per-event metadata')
    parser.add_argument('input', help='merged hits cxi file, e.g. Ery_all_hits_no_mask.cxi')
    parser.add_argument('metadata', help='per-event metadata from get_event_metadata.py')
    parser.add_argument('output', help='output cxi file')
    parser.add_argument('--mask', help='analysis good-pixel mask (h5 with /data, True = good)')
    parser.add_argument('--edge-mask', help='ASIC edge mask (h5 with /data, True = good)')
    parser.add_argument('--geom-run', type=int, default=600, help='run whose geometry file is used')
    parser.add_argument('--description', default='Single particle diffraction patterns of erythrocruorin '
                        'recorded with the AGIPD 1M detector at the SPB/SFX instrument of the European XFEL.')
    args = parser.parse_args()

    f = h5py.File(args.input, 'r')
    meta = h5py.File(args.metadata, 'r')

    trainId = f['entry_1/trainId'][()]
    cellId = f['entry_1/cellId'][()]
    N = len(trainId)
    run = meta['run'][()]
    if meta.attrs.get('input', '').split('/')[-1] != os.path.basename(args.input):
        print(f'warning: metadata file was made for {meta.attrs.get("input")}')
    if len(run) != N:
        sys.exit(f'metadata has {len(run)} events, input has {N}')
    if np.any(run == 0):
        sys.exit(f'{np.sum(run == 0)} events have no run number in {args.metadata}')
    runs = np.unique(run)

    expid = np.array([f'p{PROPOSAL:06d}_r{r:04d}_t{t}_c{c}' for r, t, c in zip(run, trainId, cellId)],
                     dtype=object)
    assert len(np.unique(expid)) == N, 'experiment identifiers are not unique'

    corner, basis, xyz, quadrant_correction, pixel_size, geom_file = load_geometry(runs, args.geom_run)

    good_pixels = f[f'{DET}/good_pixels'][()].astype(bool)
    mask, mask_defs = build_mask(good_pixels, args.mask, args.edge_mask)

    # per-run background (data_white): background_index -> run
    bindex = f['entry_1/background_index'][()]
    background = f[f'{DET}/background'][()]
    background_run = np.zeros(background.shape[0], dtype=np.uint16)
    for i in range(len(background_run)):
        r = np.unique(run[bindex == i])
        assert len(r) == 1, f'background {i} is used by runs {r}'
        background_run[i] = r[0]

    if 'pulse/background_weighting' not in meta:
        sys.exit(f'{args.metadata} has no recomputed pulse/background_weighting')
    bweight = meta['pulse/background_weighting'][()]
    background_counts = bweight * background.reshape(background.shape[0], -1).sum(1)[bindex]

    timestamps = meta['timestamp'][()]
    run_start = dict(zip(meta['runs/run'][()], meta['runs/start_time'].asstr()[()]))

    with h5py.File(args.output, 'w') as g:
        g['cxi_version'] = 160

        entry = g.create_group('entry_1')
        entry.attrs['NX_class'] = 'NXentry'
        entry.create_dataset('experiment_identifier', data=expid, dtype=h5py.string_dtype(), **GZ)
        entry['title'] = f'p{PROPOSAL:06d} {f["entry_1/sample_1/name"].asstr()[()]} hits, runs {runs[0]}-{runs[-1]}'
        entry['experiment_description'] = args.description
        entry['program_name'] = 'xfel7927/offline/make_publication_cxi.py'
        entry['start_time'] = min(run_start[r] for r in runs)
        entry['end_time'] = max(t.decode() for t in timestamps)

        sample = entry.create_group('sample_1')
        sample.attrs['NX_class'] = 'NXsample'
        sample['name'] = f['entry_1/sample_1/name'].asstr()[()]

        instrument = entry.create_group('instrument_1')
        instrument.attrs['NX_class'] = 'NXinstrument'
        instrument['name'] = 'SPB'

        # source
        source = instrument.create_group('source_1')
        source.attrs['NX_class'] = 'NXsource'
        source['name'] = 'European XFEL SASE1'
        copy_meta(source, 'energy', meta, 'train/undulator_energy',
                  description='photon energy, undulator setting (per train)')
        copy_meta(source, 'pulse_energy', meta, 'pulse/LITFRM_energyPerFrame',
                  description='XGM (SPB_XTD9) pulse energy of the pulse that produced this frame')
        copy_meta(source, 'pulse_energy_sigma', meta, 'pulse/LITFRM_energySigma')

        # XGMs
        for i, (name, desc) in enumerate((('XGM_SPB_XTD9', 'SPB_XTD9_XGM, downstream gas monitor'),
                                         ('XGM_SA1_XTD2', 'SA1_XTD2_XGM, upstream gas monitor')), 1):
            xgm = instrument.create_group(f'xgm_{i}')
            xgm.attrs['NX_class'] = 'NXcollection'
            xgm['description'] = desc
            for out, key in (('intensity', 'intensityTD'), ('intensity_sigma', 'intensitySigmaTD'),
                             ('x', 'xTD'), ('y', 'yTD'), ('x_sigma', 'xSigmaTD'), ('y_sigma', 'ySigmaTD')):
                copy_meta(xgm, out, meta, f'pulse/{name}_{key}')

        # attenuators
        for i, name in enumerate(('SA1_XTD2', 'SPB_XTD9'), 1):
            att = instrument.create_group(f'attenuator_{i}')
            att.attrs['NX_class'] = 'NXattenuator'
            att['type'] = f'{name}_ATT'
            copy_meta(att, 'attenuator_transmission', meta, f'train/attenuator/transmission_{name}')

        # electrospray injector
        es = instrument.create_group('electrospray')
        es.attrs['NX_class'] = 'NXcollection'
        for key in sorted(meta['train/electrospray']):
            copy_meta(es, key, meta, f'train/electrospray/{key}')
        for ax in 'xyz':
            copy_meta(es, f'injector_{ax}', meta, f'train/injector/{ax}')
        copy_meta(es, 'chamber_pressure', meta, 'train/chamber_pressure')

        # detector
        det = instrument.create_group('detector_1')
        det.attrs['NX_class'] = 'NXdetector'
        det['description'] = 'AGIPD 1M'
        for name, val in (('distance', DET_DIST), ('x_pixel_size', pixel_size), ('y_pixel_size', pixel_size)):
            det.create_dataset(name, data=val).attrs['units'] = 'm'
        det.create_dataset('module_identifier', data=[f'AGIPD{m:02d}T{t}' for m in range(NMODULES) for t in range(NTILES)],
                           dtype=h5py.string_dtype()).attrs['description'] = \
            'tile t = rows 64t..64t+63 of AGIPD module mm; per-pixel arrays reshape (C order) to (16, 512, 128)'

        ds = det.create_dataset('corner_position', data=corner)
        ds.attrs.update(units='m', axes='module_identifier:coordinate')
        ds = det.create_dataset('basis_vectors', data=basis)
        ds.attrs.update(units='m', axes='module_identifier:dimension:coordinate')
        ds = det.create_dataset('xyz_map', data=xyz, **GZ)
        ds.attrs.update(units='m', axes='coordinate:module_identifier:y:x',
                        description='pixel centre positions (same geometry as corner_position/basis_vectors)')
        det.create_dataset('quadrant', data=['Q1', 'Q2', 'Q3', 'Q4'], dtype=h5py.string_dtype()).attrs['description'] = \
            'quadrant q is modules 4q..4q+3, i.e. tiles (module_identifier) 32q..32q+31'
        ds = det.create_dataset('quadrant_correction', data=quadrant_correction)
        ds.attrs.update(units='m', axes='quadrant:coordinate',
                        description='alternative geometry used in the 3D reconstructions: add quadrant_correction[q] to '
                                    'the positions of the tiles in quadrant q; fitted to the data, not supported by an '
                                    'independent refinement')

        ds = det.create_dataset('mask', data=mask, **GZ)
        ds.attrs['axes'] = 'module_identifier:y:x'
        ds.attrs['bit_definitions'] = mask_defs

        # frames, reshaped into tiles; accumulate the powder on the way
        src = f[f'{DET}/data']
        dst = det.create_dataset('data', shape=(N,) + PANEL_SHAPE, dtype=src.dtype,
                                 chunks=(1,) + PANEL_SHAPE, compression='gzip', compression_opts=1)
        powder = np.zeros(PANEL_SHAPE, dtype=np.uint64)
        block = 256
        for i in tqdm(range(0, N, block), desc='copying frames'):
            frames = to_tiles(src[i:i + block])
            dst[i:i + block] = frames
            powder += frames.sum(0, dtype=np.uint64)
        dst.attrs.update(axes='experiment_identifier:module_identifier:y:x', signal=1, units='counts',
                         description='photons per pixel, not masked')
        ds = det.create_dataset('powder', data=(powder / N).astype(np.float32), **GZ)
        ds.attrs.update(axes='module_identifier:y:x', units='counts',
                        description=f'mean of data over all {N} events (unmasked, background included)')
        link_expid(det)

        per_event(det, 'start_time', timestamps.astype(object), description='train timestamp (ISO 8601)')
        per_event(det, 'run', run)
        per_event(det, 'trainId', trainId.astype(np.uint64))
        per_event(det, 'cellId', cellId.astype(np.uint16))
        copy_meta(det, 'pulseId', meta, 'pulseId')
        per_event(det, 'vds_index', f['entry_1/vds_index'][()], description='frame index in the per-run VDS file')

        # background
        ds = det.create_dataset('data_white', data=to_tiles(background), chunks=(1,) + mask.shape, **GZ)
        ds.attrs.update(axes='background_run:module_identifier:y:x', units='counts',
                        description='per-run mean of non-hit frames; zero pixels and pixels deviating by more than 20% '
                                    'from the radial profile replaced by the radial profile')
        det.create_dataset('background_run', data=background_run).attrs['description'] = 'run of each data_white'
        per_event(det, 'background_index', bindex, description='index into data_white for this event')
        per_event(det, 'background_weighting', bweight,
                  description='background of this event = background_weighting * data_white[background_index]; '
                              'per-train scale from the non-hit frames times the normalised XGM pulse energy')

        # scores
        score = det.create_group('score')
        for name, key, desc in (
                ('hit_sigma', f'{DET}/hit_sigma',
                 '(photons in hit_finding_mask - train median) / (1.4826 * (median - 25th percentile))'),
                ('photon_counts', f'{DET}/photon_counts', 'photons in the frame (per-cell mask applied)'),
                ('lit_pixels', f'{DET}/lit_pixels', 'pixels with at least one photon (per-cell mask applied)')):
            per_event(score, name, f[key][()], description=desc)
        per_event(score, 'background_counts', background_counts.astype(np.float32),
                  description='sum over pixels of background_weighting * data_white[background_index]')
        for name, key in (('facility_hitscore', 'pulse/hitfinder_hitscore'),
                          ('facility_hit_flag', 'pulse/hitfinder_hitFlag'),
                          ('facility_miss_flag', 'pulse/hitfinder_missFlag'),
                          ('facility_threshold_mu', 'train/hitfinder_threshold_mu'),
                          ('facility_threshold_sigma', 'train/hitfinder_threshold_sig'),
                          ('facility_lit_pixels', 'pulse/litpx_litPixels'),
                          ('facility_total_intensity', 'pulse/litpx_totalIntensity'),
                          ('facility_unmasked_pixels', 'pulse/litpx_unmaskedPixels')):
            copy_meta(score, name, meta, key)

        # manual selection as CXI tags
        if 'manual_selection' in f:
            names = sorted(f['manual_selection'])
            tags = np.array([f[f'manual_selection/{n}'][()] for n in names], dtype=np.int8)
            det.create_dataset('tag', data=names, dtype=h5py.string_dtype())
            ds = det.create_dataset('tags', data=tags, **GZ)
            ds.attrs['headings'] = names
            ds.attrs['axes'] = 'tag:experiment_identifier'
            ds.attrs['description'] = 'manual selection, 1 = tagged'

        # AGIPD motors
        motors = det.create_group('motors')
        motors.attrs['NX_class'] = 'NXcollection'
        for key in sorted(meta['train/agipd_motors']):
            copy_meta(motors, key, meta, f'train/agipd_motors/{key}')

        # geometry provenance
        note = det.create_group('note_1')
        note.attrs['NX_class'] = 'NXnote'
        note['file_name'] = os.path.basename(REFERENCE_GEOM)
        note['description'] = (
            'Detector geometry: facility reference geometry (this file) moved to the recorded AGIPD quadrant '
            f'motor positions with extra_geom.motors.AGIPD_1MMotors ({os.path.relpath(geom_file, ROOT)}); '
            'the motor positions are constant for all runs in this file. '
            f'Pixel positions are then shifted by {BEAM_CENTRE_SHIFT[:2].tolist()} m in x, y '
            '(mean beam centre correction from the data, consistent with an independent crystallographic '
            'refinement) and z is set to the detector distance. quadrant_correction gives the per-quadrant '
            'alternative used in the 3D reconstructions.')
        with open(os.path.join(ROOT, REFERENCE_GEOM)) as gf:
            note['data'] = gf.read()

        data_1 = entry.create_group('data_1')
        data_1.attrs['NX_class'] = 'NXdata'
        data_1['data'] = h5py.SoftLink(f'/{DET}/data')
        link_expid(data_1)

        process = entry.create_group('process_1')
        process.attrs['NX_class'] = 'NXprocess'
        process['program'] = 'make_publication_cxi.py'
        process['version'] = git_version()
        process['date'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        process['command'] = ' '.join(sys.argv)

    print(f'written {args.output}')


if __name__ == '__main__':
    main()
