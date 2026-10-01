"""
Collect per-event metadata, missing from a merged hits cxi file, from the raw/proc data.

The merged file (e.g. Ery_all_hits_no_mask.cxi) only stores trainId/cellId per
event. For every event this looks up the run, and then every quantity we have,
at the finest resolution available (per pulse where recorded, otherwise the
per-train value of the event's train). Values are converted to SI where the
facility unit is known; each dataset carries `units`, `source`, `key` attributes.

Pulse alignment (verified on r0600):
    AGIPD frame with cellId c has pulseId = LITFRM detectorPulseId[c] (= 4c here)
    LITFRM energyPerFrame[c] is the XGM pulse energy for that frame
    XGM pulse index for frame c is xgmPulseId[cumsum(nPulsePerFrame)[c] - 1]
      (indexing the XGM by cellId or pulseId, as add_pulsedata.py and
       xfel10662/make_cxi_file.py do, picks the wrong pulse for ~99% of frames)

Output layout (all per-event arrays have shape (Nevents,)):
    run, pulseId, timestamp
    pulse/<name>          per-pulse quantities (XGM, LITFRM, hit finder, lit pixels)
    pulse/background_weighting  add_background_cxi.py weighting recomputed with LITFRM pulse
                          energies (pulse/background_weighting_old reproduces the old values)
    train/<name>          per-train control values (electrospray, motors, attenuators, ...)
    runs/{run, first_train, last_train, start_time}

Usage (on maxwell, after source ../source_this_at_euxfel; use an srun/sbatch node):
    python get_event_metadata.py --list 600     # print available injector sources/keys/units
    python get_event_metadata.py                # write sidecar for the default Ery file
    python get_event_metadata.py -i <cxi> -o <h5> -s <sample substring> [--runs 600 601]
"""
import argparse
import json
import os
import sys

import numpy as np
import h5py
import extra_data
import scipy.constants as sc
from tqdm import tqdm

from constants import PREFIX, EXP_ID, DET_NAME, NMODULES

LITFRM_SRC = 'SPB_IRU_AGIPD1M1/REDU/LITFRM:output'
HITFINDER_SRC = f'{DET_NAME}/REDU/SPI_HITFINDER:output'
CORR_SRC = DET_NAME + '/CORR/{}CH0:output'
DATA_SELECTOR = 'SPB_IRU_AGIPD1M1/MDL/DATA_SELECTOR'

# Per-pulse XGM keys, read at the XGM pulse index of each frame.
XGM_KEYS = ('data.intensityTD', 'data.intensitySigmaTD', 'data.xTD', 'data.yTD',
            'data.xSigmaTD', 'data.ySigmaTD')
XGM_SRCS = {
    'XGM_SPB_XTD9': 'SPB_XTD9_XGM/XGM/DOOCS:output',   # downstream, near sample
    'XGM_SA1_XTD2': 'SA1_XTD2_XGM/XGM/DOOCS:output',   # upstream
}

# Per-train control values, written under train/<name>.
# Electrospray / aerosol injector parameters per Johan Bielecki (SPB/SFX), 2026-04-17.
CONTROL_SRCS = {
    'electrospray/CO2_capillary': ('SPB_IRU_AEROSOL/FLOW/CO2_CAPILLARY', 'measureCapacity'),
    'electrospray/CO2_chamber':   ('SPB_IRU_AEROSOL/FLOW/CO2_CHAMBER',   'measureCapacity'),
    'electrospray/DP2':           ('SPB_IRU_AEROSOL/FLOW/DP_2',          'measureCapacity'),
    'electrospray/He_chamber':    ('SPB_IRU_AEROSOL/FLOW/HE_CHAMBER',    'measureCapacity'),
    'electrospray/N2_capillary':  ('SPB_IRU_AEROSOL/FLOW/N2_CAPILLARY',  'measureCapacity'),
    'electrospray/N2_chamber':    ('SPB_IRU_AEROSOL/FLOW/N2_CHAMBER',    'measureCapacity'),
    'electrospray/liquidjet_He':  ('SPB_IRU_LIQUIDJET/FLOW/HE',          'measureCapacity'),
    'electrospray/voltage':       ('SPB_EXP_HV/MDL/SHQ1',                'channel1.voltage'),
    'electrospray/current':       ('SPB_EXP_HV/MDL/SHQ1',                'channel1.current'),
    'chamber_pressure':           ('SPB_IRU_VAC/GAUGE/GAUGE_FR_6',       'value'),
    'injector/x':                 ('SPB_IRU_INJMOV/MOTOR/X',             'actualPosition'),
    'injector/y':                 ('SPB_IRU_INJMOV/MOTOR/Y',             'actualPosition'),
    'injector/z':                 ('SPB_IRU_INJMOV/MOTOR/Z',             'actualPosition'),
    'attenuator/transmission_SA1_XTD2': ('SA1_XTD2_ATT/MDL/MAIN',        'actual.transmission'),
    'attenuator/transmission_SPB_XTD9': ('SPB_XTD9_ATT/MDL/MAIN',        'actual.transmission'),
    'undulator_energy':           ('SPB_XTD2_UND/DOOCS/ENERGY',          'actualPosition'),
    'xgm_wavelength_used':        ('SPB_XTD9_XGM/XGM/DOOCS',             'pulseEnergy.wavelengthUsed'),
}
# AGIPD quadrant motors and detector z stage (via the slow data selector, as
# used by extra.components.AGIPD1MQuadrantMotors for the per-run geometry)
for q in range(1, 5):
    for m in range(1, 3):
        CONTROL_SRCS[f'agipd_motors/Q{q}M{m}'] = (DATA_SELECTOR, f'spbIruAgipd1MMotorQ{q}M{m}.actualPosition')
CONTROL_SRCS['agipd_motors/z_stepper'] = (DATA_SELECTOR, 'spbIruAgipd1MMotorZStepper.actualPosition')

# Units not recorded in the raw data (checked on r0600); confirm with the beamline.
# The flow channels report measureCapacity in the Bronkhorst controller's
# capacity unit (measure is % of full scale).
MISSING_UNITS = {
    'chamber_pressure': 'unknown (likely mbar)',
    'attenuator/transmission_SA1_XTD2': 'dimensionless',
    'attenuator/transmission_SPB_XTD9': 'dimensionless',
}

# facility unit -> (SI unit, scale)
TO_SI = {
    'μJ': ('J', 1e-6), 'uJ': ('J', 1e-6),
    'mm': ('m', 1e-3),
    'nm': ('m', 1e-9),
    'keV': ('J', 1e3 * sc.e),
}

# search patterns for --list
INJECTOR_PATTERNS = ('AEROSOL', 'LIQUIDJET', 'EXP_HV', 'INJ', 'ESI', 'NOZZLE', 'FLOW', 'PRESS', 'GAUGE')


def get_units(kd):
    """Unit string for a KeyData, e.g. 'mm', or '' if not recorded."""
    try:
        return str(kd.units or '')
    except Exception:
        return ''


def to_si(vals, units):
    if units in TO_SI:
        si, scale = TO_SI[units]
        return vals * scale, si
    return vals, units


def list_sources(run):
    dc = extra_data.open_run(int(EXP_ID), run, data='raw')
    srcs = sorted(s for s in dc.all_sources if any(p in s for p in INJECTOR_PATTERNS))
    for src in srcs:
        print(src)
        for key in sorted(dc[src].keys(inc_timestamps=False)):
            kd = dc[src, key]
            try:
                v = kd.ndarray()
                if v.dtype.kind not in 'fiub':
                    continue
                stats = f'min={np.min(v):.4g} max={np.max(v):.4g}'
            except Exception as e:
                stats = f'({e!r})'
            print(f'    {key:50s} [{get_units(kd)}] {stats}')


def get_runs(sample):
    with open(f'{PREFIX}scratch/log/run_table.json') as f:
        table = json.load(f)
    return sorted(v['Run number'] for v in table.values()
                  if isinstance(v, dict) and sample in str(v.get('Sample')))


def iso(t):
    return '' if np.isnat(t) else np.datetime_as_string(t, unit='ns') + 'Z'


def pair_key(train, pulse):
    return (np.asarray(train, dtype=np.int64) << 32) | np.asarray(pulse, dtype=np.int64)


def lookup(src_keys, target_keys):
    """index into src_keys for each target key, and a mask of which were found"""
    order = np.argsort(src_keys)
    pos = np.clip(np.searchsorted(src_keys[order], target_keys), 0, max(len(order) - 1, 0))
    found = (len(order) > 0) & (src_keys[order][pos] == target_keys)
    return order[pos], found


class Columns:
    """per-event output columns, filled run by run"""
    def __init__(self, Nevents):
        self.N = Nevents
        self.data = {}
        self.attrs = {}

    def set(self, name, idx, vals, units='', source='', key='', dtype=np.float32, fill=np.nan):
        if name not in self.data:
            self.data[name] = np.full(self.N, fill, dtype=dtype)
        self.data[name][idx] = vals
        self.attrs.setdefault(name, dict(units=units, source=source, key=key))


def per_train_rows(kd, tids):
    """row of kd.ndarray() for each train id in tids, and found mask"""
    coords = np.asarray(kd.train_id_coordinates(), dtype=np.int64)
    return lookup(coords, tids)


def process_run(dc, ev, cellId, trainId, out):
    """fill out for events ev (indices into the merged file) belonging to run dc"""
    t = trainId[ev]
    c = cellId[ev].astype(np.int64)

    # LITFRM: pulseId, XGM energy per frame and the XGM pulse index of each frame
    xgm_index = None
    pulseId = None
    if LITFRM_SRC in dc.all_sources:
        lf = dc[LITFRM_SRC]
        rows, ok = per_train_rows(lf['data.detectorPulseId'], t)
        dpid = lf['data.detectorPulseId'].ndarray()
        npf  = lf['data.nPulsePerFrame'].ndarray().astype(np.int64)
        xpid = lf['data.xgmPulseId'].ndarray().astype(np.int64)
        rank = np.cumsum(npf, axis=1) - 1

        r, cc, e = rows[ok], c[ok], ev[ok]
        pulseId = np.full(len(ev), -1, dtype=np.int64)
        pulseId[ok] = dpid[r, cc]
        out.set('pulseId', e, dpid[r, cc], dtype=np.int64, fill=-1, source=LITFRM_SRC, key='data.detectorPulseId')

        has_pulse = npf[r, cc] > 0
        xgm_index = np.full(len(ev), -1, dtype=np.int64)
        xgm_index[np.where(ok)[0][has_pulse]] = xpid[r[has_pulse], rank[r[has_pulse], cc[has_pulse]]]

        for key in ('data.energyPerFrame', 'data.energySigma'):
            kd = lf[key]
            vals, units = to_si(kd.ndarray()[r, cc].astype(float), get_units(kd))
            out.set(f'pulse/LITFRM_{key.split(".")[-1]}', e, vals, units, LITFRM_SRC, key)
    else:
        print(f'WARNING: {LITFRM_SRC} not in run', file=sys.stderr)

    # XGMs at the frame's XGM pulse index
    if xgm_index is not None:
        for name, src in XGM_SRCS.items():
            if src not in dc.all_sources:
                continue
            for key in XGM_KEYS:
                kd = dc[src, key]
                rows, ok = per_train_rows(kd, t)
                ok &= xgm_index >= 0
                arr = kd.ndarray()
                vals, units = to_si(arr[rows[ok], xgm_index[ok]].astype(float), get_units(kd))
                out.set(f'pulse/{name}_{key.split(".")[-1]}', ev[ok], vals, units, src, key)

    # facility hit finder, matched on (trainId, pulseId)
    if pulseId is not None and HITFINDER_SRC in dc.all_sources:
        hf = dc[HITFINDER_SRC]
        hkeys = pair_key(hf['data.trainId'].ndarray(), hf['data.pulseId'].ndarray())
        idx, ok = lookup(hkeys, pair_key(t, pulseId))
        for key, dtype in (('data.hitscore', np.float32), ('data.hitFlag', np.int8), ('data.missFlag', np.int8)):
            fill = np.nan if dtype is np.float32 else -1
            out.set(f'pulse/hitfinder_{key.split(".")[-1]}', ev[ok], hf[key].ndarray()[idx[ok]],
                    '', HITFINDER_SRC, key, dtype=dtype, fill=fill)
        for key in ('threshold.mu', 'threshold.sig'):
            rows, ok = per_train_rows(hf[key], t)
            out.set(f'train/hitfinder_{key.replace(".", "_")}', ev[ok], hf[key].ndarray()[rows[ok]],
                    '', HITFINDER_SRC, key)

    # lit pixels / total intensity from the correction pipeline, summed over modules
    if pulseId is not None:
        sums = {k: np.zeros(len(ev)) for k in ('litPixels', 'totalIntensity', 'unmaskedPixels')}
        nmod = np.zeros(len(ev), dtype=int)
        for m in range(NMODULES):
            src = CORR_SRC.format(m)
            if src not in dc.all_sources or 'litpx.litPixels' not in dc[src].keys():
                continue
            s = dc[src]
            idx, ok = lookup(pair_key(s['litpx.trainId'].ndarray(), s['litpx.pulseId'].ndarray()),
                             pair_key(t, pulseId))
            for k in sums:
                sums[k][ok] += s[f'litpx.{k}'].ndarray()[idx[ok]]
            nmod += ok
        ok = nmod > 0
        for k, v in sums.items():
            out.set(f'pulse/litpx_{k}', ev[ok], v[ok], '', CORR_SRC.format('*'), f'litpx.{k}')
        out.set('pulse/litpx_nmodules', ev, nmod, '', CORR_SRC.format('*'), 'modules summed', dtype=np.int8, fill=0)

    # per-train control values
    for name, (src, key) in CONTROL_SRCS.items():
        if src not in dc.all_sources or key not in dc[src].keys(inc_timestamps=False):
            continue
        kd = dc[src, key]
        rows, ok = per_train_rows(kd, t)
        units = get_units(kd) or MISSING_UNITS.get(name, 'unknown')
        vals, units = to_si(kd.ndarray().ravel()[rows[ok]].astype(float), units)
        out.set(f'train/{name}', ev[ok], vals, units, src, key)


# minimum reliable pulse energy reading (J), as in add_background_cxi.py
EMIN = 1e-3


def background_weighting(run, dc, ev, vds_index, trainId, cellId, out):
    """Recompute add_background_cxi.py's per-frame background weighting
        b_d = a_t(d) * e_d / <e>
    with e_d the XGM pulse energy of frame d taken via LITFRM, instead of the
    mis-indexed events-file pulse_energy. a_t (per train, from the photon counts
    of the misses and the per-run non-hit powder) is unchanged. The old weighting
    is recomputed too, as a check against the merged file."""
    events_file = f'{PREFIX}scratch/events/r{run:04d}_events.h5'
    back_file = f'{PREFIX}scratch/powder/r{run:04d}_powder_is_hit_False_per_pixel.h5'
    if not (os.path.exists(events_file) and os.path.exists(back_file)):
        print(f'WARNING: no events or background file for run {run}, background weighting skipped', file=sys.stderr)
        return
    with h5py.File(events_file) as f:
        tid = f['trainId'][()].astype(np.int64)
        cid = f['cellId'][()]
        cid = (cid[:, 0] if cid.ndim == 2 else cid).astype(np.int64)
        misses = ~f['is_hit'][()]
        photon_counts = f['total_intens'][()]
        e_old = f['pulse_energy'][()]
    with h5py.File(back_file) as f:
        back_counts = np.sum(f['data'][()])

    v = vds_index[ev]
    assert np.array_equal(tid[v], trainId[ev]) and np.array_equal(cid[v], cellId[ev]), \
        f'events file of run {run} does not match the merged file vds_index'

    # a_t: mean miss photon counts in the train / mean background counts (-1 if no misses)
    ut, inv = np.unique(tid, return_inverse=True)
    s = np.bincount(inv, photon_counts * misses)
    n = np.bincount(inv, misses)
    a_t = np.where(n > 0, s / (back_counts * np.clip(n, 1, None)), -1.)
    a_d = a_t[inv]

    def normalise(e):
        e = np.array(e, dtype=float)
        m = e > EMIN
        e[m] /= np.mean(e[m])
        e[~m] = 1
        return e

    # pulse energy of every frame in the run via LITFRM
    e_new = np.full(len(tid), np.nan)
    kd = dc[LITFRM_SRC, 'data.energyPerFrame']
    rows, ok = per_train_rows(kd, tid)
    vals, units = to_si(kd.ndarray().astype(float), get_units(kd))
    assert units == 'J', units
    e_new[ok] = vals[rows[ok], cid[ok]]

    src = f'{events_file} + {back_file} + {LITFRM_SRC}'
    out.set('pulse/background_weighting', ev, (normalise(e_new) * a_d)[v], '', src,
            'a_t * e_d / <e>, e from LITFRM energyPerFrame')
    out.set('pulse/background_weighting_old', ev, (normalise(e_old) * a_d)[v], '', src,
            'a_t * e_d / <e>, e from events pulse_energy (as add_background_cxi.py)')
    out.set('pulse/background_train_factor', ev, a_d[v], '', src, 'a_t')


def main():
    parser = argparse.ArgumentParser(description='Collect per-event run number, timestamp, pulse-resolved and per-train metadata for a merged cxi file')
    parser.add_argument('--list', type=int, metavar='RUN', help='list injector sources and keys for RUN and exit')
    parser.add_argument('-i', '--input', default=f'{PREFIX}scratch/saved_hits/Ery_all_hits_no_mask.cxi')
    parser.add_argument('-o', '--output', default=f'{PREFIX}scratch/saved_hits/Ery_event_metadata.h5')
    parser.add_argument('-s', '--sample', default='Ery', help='substring of the run table sample name')
    parser.add_argument('--runs', type=int, nargs='+', help='only process these runs (for testing)')
    args = parser.parse_args()

    if args.list is not None:
        list_sources(args.list)
        return

    with h5py.File(args.input) as f:
        trainId = f['entry_1/trainId'][()].astype(np.int64)
        cellId  = f['entry_1/cellId'][()].astype(np.int64)
        vds_index = f['entry_1/vds_index'][()].astype(np.int64)
    Nevents = len(trainId)
    print(f'{Nevents} events in {args.input}')

    runs = args.runs or get_runs(args.sample)
    print(f'{len(runs)} candidate runs')

    out = Columns(Nevents)
    event_run  = np.zeros(Nevents, dtype=np.uint16)
    event_time = np.full(Nevents, '', dtype=object)
    run_rows = []

    for run in tqdm(runs):
        try:
            dc = extra_data.open_run(int(EXP_ID), run, data='all')
        except Exception as e:
            print(f'WARNING: could not open run {run}: {e!r}', file=sys.stderr)
            continue

        tids = np.asarray(dc.train_ids, dtype=np.int64)
        ev = np.where((trainId >= tids[0]) & (trainId <= tids[-1]))[0]
        if len(ev) == 0:
            continue

        ts = dc.train_timestamps()
        run_rows.append((run, tids[0], tids[-1], iso(ts[0])))
        event_run[ev] = run

        pos, ok = lookup(tids, trainId[ev])
        event_time[ev[ok]] = [iso(x) for x in ts[pos[ok]]]

        # only the trains containing events
        if LITFRM_SRC in dc.all_sources:
            background_weighting(run, dc, ev, vds_index, trainId, cellId, out)

        # only the trains containing events
        dc = dc.select_trains(extra_data.by_id[np.unique(trainId[ev])])
        process_run(dc, ev, cellId, trainId, out)

    missing = np.sum(event_run == 0)
    if missing:
        print(f'WARNING: {missing} events not matched to any run', file=sys.stderr)

    with h5py.File(args.output, 'w') as f:
        f.attrs['input'] = args.input
        f['run'] = event_run
        f.create_dataset('timestamp', data=list(event_time.astype(str)), dtype=h5py.string_dtype())
        for name, v in out.data.items():
            ds = f.create_dataset(name, data=v, compression='gzip')
            for k, a in out.attrs[name].items():
                ds.attrs[k] = a
        r = list(zip(*run_rows))
        f['runs/run']         = np.array(r[0], dtype=np.uint16)
        f['runs/first_train'] = np.array(r[1], dtype=np.uint64)
        f['runs/last_train']  = np.array(r[2], dtype=np.uint64)
        f.create_dataset('runs/start_time', data=list(r[3]), dtype=h5py.string_dtype())

    print(f'{len(run_rows)} runs contribute events; written {args.output}')


if __name__ == '__main__':
    main()
