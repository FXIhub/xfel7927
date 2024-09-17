# xfel7927
Data analysis for the SPB/SFX experiment "xfel7927"

## Files
```text
Experiment root  : /gpfs/exfel/exp/SPB/202405/p007927
AGIPD Data       : /gpfs/exfel/exp/SPB/202405/p007927/raw/r*/RAW*AGIPD*.h5
Meta  Data       : /gpfs/exfel/exp/SPB/202405/p007927/raw/r*/RAW*DA*.h5
AGIPD (corrected): /gpfs/exfel/exp/SPB/202405/p007927/proc/r*/COR*AGIPD*.h5
Analysis scripts : /gpfs/exfel/exp/SPB/202405/p007927/usr/Shared/<user>/
Analysis files   : /gpfs/exfel/exp/SPB/202405/p007927/scratch/<user>/
```

## Pipeline
```text
XFEL Data (/raw)
   └► facility calibrations (/proc/r*/CORR*)           
   └► VDS files (/scratch/vds/r*.cxi)
      └► event info (/scratch/events/r*_events.h5)
         └► cxi files for hits (/scratch/hits/r*_hits.cxi)
            ├► static EMC (/scratch/emc_static/r*_emc_static.h5)
            └► 2D EMC (/scratch/emc_2D/r*_emc_2D.h5)
               └► classification (write to event info)
                  └► 3D EMC (/scratch/emc_3D/r*_emc_3D.h5)
```

## Scripts
These are located in `offline/`
```text
XFEL Data (/raw)
   └► facility calibrations (automatically triggered by EuXFEL)           
      └► VDS files (submit_vds.sh)
         └► event info (submit_events.sh)
            └► cxi files for hits (submit_cxi.sh)
               ├► static EMC (...)
               └► 2D EMC (...)
                  └► classification (...)
                     └► 3D EMC (...)
```
