# 2D-EMC: Expand Maximise Compress with inplane rotations

## Branches
- simple: uses single core python on minimal working examples with hard coded parameters. 
- main: simple + mpi, input and output parameters and structured code. 
- optimise: optimised with opencl, code becomes unreadable.
- optimise-fast: further optimised and with many comprimises to make the code fast, even if the results are dodgy.

## Notation
- I[t, mi, mj] class t at pixel mi, mj
- K[d, i] photon counts for frame d at pixel i
- W[t, r, i] estimated photon counts for class t, rotation angle r and pixel i
- w[d] fluence estimates
- B[l, i] estimated photon counts for background class l
- b[d, l] estimated background weights
- T[d, t, r, i] = w[d] W[t, r, i] + sum_l b[d, l] B[l, i] estimated model for frame
- P[d, t, r] = R[d, t, r] / sum_tr R[d, t, r] probability matrix 

- R[d, t, r] likelihoods
    R = Prod_i T^K e^-T 
    log R = sum_i [K log(T)  - T]

- E expectation value 
    E = sum_dtr P sum_i log( K logT - T ) 
- (mi, mj) = R_r[i] rotation matrix for rotation angle indexed by r, maps dector pixel i to model pixels (mi, mj)
    R(t[r]) = [[cos(t) -sin(t)], [sin(t) cos(t)]]
    t[r]    = 2 pi r / M 
    (u, v)  = R . (x[i], y[i])
    mi[u]    = i0 + u / dx
    mj[v]    = j0 + v / dx
- dx is the model square voxel side length
    
- M = number of in-plane angles
- C[i] = solid angle and polarisation factor 
    C[i]     = Omega[i] P[i]
    P[i]     = 1 - (x[i] / radius[i])^2 (for polarisation along x-axis)
    Omega[i] = pixel_area[i] z[i] / radius[i]^3 (solid angle correction)

## Expand
- W[t, r, i] = C[i] I[t, R_r[i].mi, R_r[i].mj]

## Maximise
Line search to find points where gradients are zero:

- gW[t, r, i] = sum_d w[d] P[d, t, r] (K[d, i] / T[d, t, r, i] - 1)
    = sum_d P[d, t, r] K[d, i] w[d] / T[d, t, r, i]  - \sum_d w[d] P[d, t, r] 

- gw[d] = sum_tri W[t, r, i] P[d, t, r] (K[d, i] / T[d, t, r, i] - 1)
    = sum_tri P[d, t, r] K[d, i] W[t, r, i] / T[d, t, r, i]  - \sum_tri W[t, r, i] P[d, t, r] 

- gb[d, l] = sum_tri B[l, i] P[d, t, r] (K[d, i] / T[d, t, r, i] - 1)
    = sum_tri P[d, t, r] K[d, i] B[l, i] / T[d, t, r, i]  - \sum_tri B[l, i] P[d, t, r] 
    = sum_tri P[d, t, r] K[d, i] B[l, i] / T[d, t, r, i]  - \sum_i B[l, i] 

## Compress
- sum_ri (I[t, mi[i], mj[j]] += W[t, r, i] / C[i])
- sum_ri (O[t, mi[i], mj[j]] += i)
- I /= O

## Notes
### update W
updating tomograms per class, rotation, looping over frames, kernel call over pixel chunks:
    - each work item operates on one pixel 
    - broadcast P[d], w[d], b[d], coallesed read of K[i], B[i]
    - each of these is read once per frame, class, rotation, pixel and iter (1e15)
    - that is 4 floats = 3600 TB 
    - and 1 char (K) = 914 TB (although a char is probably just as slow as a float)
    - for a total of 15.5 PB of reading
    - 15.5 PB / (1.5 terabyetes/s A100) = 3 hours
    - calculate T, 3 flops per iter = 3 petaFLOP
    - + the rest = 9 petaFLOP
    - 9 petaFLOP / (19.5 teraFLOPS A100) = 7.9 minutes

amaizing if that 7.9 minutes is correct. I'm guessing I am way below that though factor of 10-100 1.3 13 hours.

On paper, it's the global reads I have to reduce, but perhaps caching makes this better already.

### update w
update w per frame, one work group per frame w work group reduce for summing:
    - each work item loops over a subset of pixels
    - main issue is that all W are recalculated for each frame
    - that's at least 4 petaFlOP 
    - could reduce this by looping over frames and pixels within a work group
    - This is essentially the same thing that happens in the logR calc, compare 
      each tomogram to each frame, but multiplied by iter
    
    
updating tomograms per class, rotation, looping over frames, kernel call over pixel chunks:
