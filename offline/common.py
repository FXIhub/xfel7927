import pathlib
import os

from constants import EXP_ID 

# .../xfel7927/
root = pathlib.Path(__file__).parent.parent.resolve()

def get_geom(run_no):
    return make_geom(run_no)

def make_geom(run_no):
    out_fnam = root.joinpath(f'geom/r{run_no:04}.geom').resolve()
    ref_fnam = root.joinpath('geom/agipd_p008316_r0024_v04.geom').resolve()
    
    if not os.path.exists(out_fnam) :
        from extra_geom import AGIPD_1MGeometry
        from extra_geom.motors import AGIPD_1MMotors
        from extra.components import AGIPD1MQuadrantMotors
        from extra_data import open_run

        ref_geom = AGIPD_1MGeometry.from_crystfel_geom(ref_fnam)
        run = open_run(int(EXP_ID), run_no)
        motors = AGIPD1MQuadrantMotors(run)
        tracker = AGIPD_1MMotors(ref_geom)
        geom = tracker.geom_at(motors.most_frequent_positions())
        
        geom.write_crystfel_geom(out_fnam)

    return out_fnam




