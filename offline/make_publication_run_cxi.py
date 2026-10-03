"""
Write every frame of one run (e.g. a gas background run, no sample) to a CXI
v1.6 file for publication, in the same layout as make_publication_cxi.py.

    python make_publication_run_cxi.py <run> [-m <event_metadata.h5>] [-o <output.cxi>] [-n <nproc>]

Inputs (defaults under {PREFIX}scratch/)
    vds/r{run}.cxi                      frames (photons), trainId, cellId
    publication/r{run}_metadata.h5      get_event_metadata.py -i vds/r{run}.cxi --runs {run}
    det/r{run}_mask.h5                  per-cell masks of the run (union -> mask bit 0x00010000)
    Ery/recon_3D_nosym/mask.h5          analysis mask (mask bit CXI_PIXEL_IS_BAD)
    geom/r{run}.geom                    AGIPD geometry (reference geometry + quadrant motors)
    log/run_table.json                  sample name (gas condition) of the run

Output: publication/p007927_r{run:04d}_gas_background.cxi

Layout (see make_publication_cxi.py for the meaning of the shared entries,
Detector tiles, Geometry and Mask). N is the number of frames in the run's VDS.

    /
      cxi_version = 160                                            int64
      entry_1/                                                     NXentry
        experiment_identifier   (N,)            string             p007927_r{run:04d}_t{trainId}_c{cellId}
        title                   scalar          string             "p007927 run {run} gas background ({condition}), all frames"
        experiment_description  scalar          string             --description
        program_name            scalar          string             "xfel7927/offline/make_publication_run_cxi.py"
        start_time, end_time    scalar          string             ISO 8601, first and last event
        sample_1/                                                  NXsample
          name                  scalar          string             gas condition from the run table, e.g. "CO2_N2_He"
          description           scalar          string             "no sample (gas background) ..."
        instrument_1/                                              NXinstrument, as make_publication_cxi.py:
          name, source_1/, xgm_1/, xgm_2/, attenuator_1/, attenuator_2/, electrospray/
                                                                   (electrospray/ holds the gas flows per event)
          detector_1/                                              NXdetector, as make_publication_cxi.py:
            description, distance, x_pixel_size, y_pixel_size, module_identifier,
            corner_position, basis_vectors, xyz_map, quadrant, quadrant_correction, mask
            data                (N, 128, 64, 128) uint8  counts    photons per pixel (clipped to 0..255), NOT masked
            powder              (128, 64, 128)  float32  counts    mean of data over all N frames
            start_time          (N,)            string             approximate train time, ISO 8601 (+-1 s)
            run                 (N,)            uint16
            trainId             (N,)            uint64
            cellId              (N,)            uint16
            pulseId             (N,)            int64
            vds_index           (N,)            int64              frame index in the run's VDS file (0..N-1)
            score/
              photon_counts     (N,)            float32            photons in the frame over pixels with mask == 0
              lit_pixels        (N,)            float32            pixels with >= 1 photon and mask == 0
              facility_*        (N,)                               as make_publication_cxi.py
              experiment_identifier -> /entry_1/experiment_identifier
            motors/, note_1/                                       as make_publication_cxi.py
            experiment_identifier -> /entry_1/experiment_identifier
        data_1/                                                    NXdata
          data -> /entry_1/instrument_1/detector_1/data
          experiment_identifier -> /entry_1/experiment_identifier
        process_1/                                                 NXprocess (program, version, date, command)

Axes attributes as in make_publication_cxi.py (no data_white, background_*,
hit_sigma, background_counts or tags in these files).

Frames are read from the VDS in blocks by worker processes, which clip to uint8,
reshape to tiles, compress (deflate level 1) and accumulate the scores and the
powder; the main process writes the compressed chunks in order.
"""
import argparse
import io
import json
import multiprocessing as mp
import os
import zlib

import numpy as np
import h5py
from tqdm import tqdm

from constants import PREFIX
from make_publication_cxi import (PROPOSAL, DET, PANEL_SHAPE, GZ, to_tiles, per_event, link_expid,
                                  load_geometry, build_mask, edge_mask, approx_timestamps,
                                  write_header, write_beam_and_injector, write_detector_static,
                                  write_event_ids, write_facility_scores, write_motors_and_note,
                                  write_footer, write_powder, frames_attrs)

VDS_DATA = '/entry_1/instrument_1/detector_1/data'
PROGRAM = 'make_publication_run_cxi.py'

GAS_DESCRIPTION = 'no sample (gas background); gas condition "{}" as logged in the run table; the injector ' \
                  'gas flows of every event are in instrument_1/electrospray'

# worker state, set by init_worker
_vds = None
_good = None


def init_worker(vds_file, good):
    global _vds, _good
    _vds = h5py.File(vds_file, 'r')[VDS_DATA]
    _good = good


def process_block(span):
    """read frames start:stop, return their compressed chunks, scores and powder sum"""
    start, stop = span
    frames = to_tiles(np.clip(_vds[start:stop], 0, 255).astype(np.uint8))
    chunks = [zlib.compress(np.ascontiguousarray(fr).tobytes(), 1) for fr in frames]
    photons = (frames * _good).sum(axis=(1, 2, 3), dtype=np.uint64).astype(np.float32)
    lit = ((frames > 0) & _good).sum(axis=(1, 2, 3)).astype(np.float32)
    return start, chunks, photons, lit, frames.sum(0, dtype=np.uint64)


def sample_name(run):
    with open(f'{PREFIX}scratch/log/run_table.json') as f:
        table = json.load(f)
    for v in table.values():
        if isinstance(v, dict) and v.get('Run number') == run:
            return v['Sample']
    raise ValueError(f'run {run} not in the run table')


def main():
    parser = argparse.ArgumentParser(description='Write every frame of one run to a publication CXI file')
    parser.add_argument('run', type=int)
    parser.add_argument('-m', '--metadata', help='per-event metadata (default scratch/publication/r{run}_metadata.h5)')
    parser.add_argument('-o', '--output', help='output file (default scratch/publication/p007927_r{run}_gas_background.cxi)')
    parser.add_argument('--vds', help='VDS file (default scratch/vds/r{run}.cxi)')
    parser.add_argument('--mask', default=f'{PREFIX}scratch/Ery/recon_3D_nosym/mask.h5',
                        help='analysis good-pixel mask (h5 with /data, True = good)')
    parser.add_argument('--cell-mask', help='per-cell mask file (default scratch/det/r{run}_mask.h5)')
    parser.add_argument('-n', '--nproc', type=int, default=8, help='worker processes')
    parser.add_argument('--block', type=int, default=64, help='frames per worker task')
    parser.add_argument('--max-frames', type=int, help='only write the first frames (testing)')
    parser.add_argument('--description', default='Gas background (no sample) diffraction recorded with the AGIPD 1M '
                        'detector at the SPB/SFX instrument of the European XFEL, in the same conditions as the '
                        'single particle imaging runs of this proposal.')
    args = parser.parse_args()
    r = args.run
    out_dir = f'{PREFIX}scratch/publication'
    args.metadata = args.metadata or f'{out_dir}/r{r:04d}_metadata.h5'
    args.output = args.output or f'{out_dir}/p{PROPOSAL:06d}_r{r:04d}_gas_background.cxi'
    args.vds = args.vds or f'{PREFIX}scratch/vds/r{r:04d}.cxi'
    args.cell_mask = args.cell_mask or f'{PREFIX}scratch/det/r{r:04d}_mask.h5'
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    with h5py.File(args.vds) as f:
        trainId = f['entry_1/trainId'][()]
        cellId = f['entry_1/cellId'][:, 0]
        N = f[VDS_DATA].shape[0]
    N = min(N, args.max_frames or N)
    trainId, cellId = trainId[:N], cellId[:N]
    vds_index = np.arange(N, dtype=np.int64)

    meta = h5py.File(args.metadata, 'r')
    if meta['run'].shape[0] < N or np.any(meta['run'][:N] != r):
        raise ValueError(f'{args.metadata} does not cover the frames of run {r}')
    if meta['run'].shape[0] > N:
        meta = truncate(meta, N)
    run = meta['run'][()]

    expid = np.array([f'p{PROPOSAL:06d}_r{r:04d}_t{t}_c{c}' for t, c in zip(trainId, cellId)], dtype=object)

    *geometry, geom_file = load_geometry([r], r)

    with h5py.File(args.cell_mask) as f:
        good_pixels = f['entry_1/good_pixels'][()].astype(bool).any(axis=0)
    mask, mask_defs = build_mask(good_pixels, args.mask, edge_good=edge_mask())
    mask_defs = [d.replace('every run', 'this run') for d in mask_defs]
    good = mask == 0

    timestamps = approx_timestamps(run, trainId, meta['timestamp'].asstr()[()])
    condition = sample_name(r)

    with h5py.File(args.output, 'w') as g:
        entry, instrument = write_header(
            g, expid, f'p{PROPOSAL:06d} run {r} gas background ({condition}), all frames',
            args.description, PROGRAM, timestamps, condition, GAS_DESCRIPTION.format(condition))
        write_beam_and_injector(instrument, meta)
        det = write_detector_static(instrument, geometry, mask, mask_defs)

        dst = det.create_dataset('data', shape=(N,) + PANEL_SHAPE, dtype=np.uint8,
                                 chunks=(1,) + PANEL_SHAPE, compression='gzip', compression_opts=1)
        photons = np.zeros(N, dtype=np.float32)
        lit = np.zeros(N, dtype=np.float32)
        powder = np.zeros(PANEL_SHAPE, dtype=np.uint64)
        spans = [(i, min(i + args.block, N)) for i in range(0, N, args.block)]
        with mp.Pool(args.nproc, initializer=init_worker, initargs=(args.vds, good)) as pool:
            for start, chunks, p, l, s in tqdm(pool.imap(process_block, spans), total=len(spans), desc='frames'):
                for i, c in enumerate(chunks):
                    dst.id.write_direct_chunk((start + i, 0, 0, 0), c)
                photons[start:start + len(chunks)] = p
                lit[start:start + len(chunks)] = l
                powder += s
        frames_attrs(dst)
        write_powder(det, powder, N)
        link_expid(det)

        write_event_ids(det, meta, timestamps, run, trainId, cellId, vds_index)

        score = det.create_group('score')
        per_event(score, 'photon_counts', photons, description='photons in the frame over pixels with mask == 0')
        per_event(score, 'lit_pixels', lit, description='pixels with at least one photon and mask == 0')
        write_facility_scores(score, meta)

        write_motors_and_note(det, meta, geom_file)
        write_footer(entry, PROGRAM)

    print(f'written {args.output}')


def truncate(meta, N):
    """in-memory copy of the metadata file with the per-event datasets cut to the first N events"""
    M = meta['run'].shape[0]
    out = h5py.File(io.BytesIO(), 'w')
    for k, v in meta.attrs.items():
        out.attrs[k] = v

    def copy(name, obj):
        if isinstance(obj, h5py.Dataset):
            ds = out.create_dataset(name, data=obj[:N] if obj.shape[:1] == (M,) else obj[()], dtype=obj.dtype)
            for k, v in obj.attrs.items():
                ds.attrs[k] = v
    meta.visititems(copy)
    return out


if __name__ == '__main__':
    main()
