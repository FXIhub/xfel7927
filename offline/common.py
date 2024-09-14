import pathlib

# .../xfel7927/
root = pathlib.Path(__file__).parent.parent.resolve()

def get_geom(run):
    return root.joinpath('geom/motor_p7076_from_4462.geom').resolve()
