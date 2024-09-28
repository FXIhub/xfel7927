# Graphical Mask Maker for DSSC detector 
Based on CsPadMaskMaker.

This provides a geometry corrected view of the hexagonal pixels of the DSSC detector, each pixel is placed in x,y space then the surrounding area is 'filled' to provide an approximate representation of the hexagonal grid with a small number of pixels. Masked pixels are in data-space, so when masking a pixel on screen this will lookup the corresponding data pixel and then map that change to all on-screen pixels representing a single hexagonal pixel.

The boolean mask is saved to:
    mask.h5
        entry_1/good_pixels


## Example
exfel-python has all of the dependencies built in:
```
$ module load exfel exfel-python
# python maskmaker.py /gpfs/exfel/exp/SPB/202405/p007927/scratch/saved_hits/r0065_hits.cxi/entry_1/instrument_1/detector_1/background -g /gpfs/exfel/exp/SPB/202405/p007927/scratch/saved_hits/r0065_hits.cxi/entry_1/instrument_1/detector_1/xyz_map
```
