import pathlib

# .../xfel7927/
root = pathlib.Path(__file__).parent.parent.resolve()

def get_geom(run):
    return root.joinpath('geom/r0035_powder.geom').resolve()
