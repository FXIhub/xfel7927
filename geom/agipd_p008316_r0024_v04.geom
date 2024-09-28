; AGIPD-1M geometry file written by EXtra-geom 1.12.0
; You may need to edit this file to add:
; - data and mask locations in the file
; - mask_good & mask_bad values to interpret the mask
; - adu_per_eV & photon_energy
; - clen (detector distance)
;
; See: http://www.desy.de/~twhite/crystfel/manual-crystfel_geometry.html

;XGEOM MOTORS=4,2
;XGEOM MOTOR_Q1=56.399078369140625,83.81175231933594
;XGEOM MOTOR_Q2=-15.18103313446045,-72.31230163574219
;XGEOM MOTOR_Q3=44.998077392578125,697.5655517578125
;XGEOM MOTOR_Q4=-10.368175506591797,-30.78021812438965

data = /entry_1/instrument_1/detector_1/data

mask0_data = /entry_1/data_1/mask ;
mask0_goodbits = 0 ;
mask0_badbits = 65535 ;

mask1_file = /gpfs/exfel/exp/SPB/202405/p007927/usr/Shared/mask/mask_hvoff_20240807.h5
mask1_data = /entry_1/data_1/mask ;
mask1_goodbits = 0 ;
mask1_badbits = 1 ;

mask2_file = /gpfs/exfel/exp/SPB/202405/p007927/usr/Shared/mask/mask_p900463_extra_asic.h5
mask2_data = /entry_1/data_1/mask ;
mask2_goodbits = 0 ;
mask2_badbits = 1 ;

dim0 = %
res = 5000.0 ; pixels per metre

; Beam energy in eV
photon_energy = 11699.0

; Camera length, aka detector distance
clen = 0.12322

; Analogue Digital Units per eV
adu_per_eV = 1.0

rigid_group_p0 = p0a0,p0a1,p0a2,p0a3,p0a4,p0a5,p0a6,p0a7
rigid_group_p1 = p1a0,p1a1,p1a2,p1a3,p1a4,p1a5,p1a6,p1a7
rigid_group_p2 = p2a0,p2a1,p2a2,p2a3,p2a4,p2a5,p2a6,p2a7
rigid_group_p3 = p3a0,p3a1,p3a2,p3a3,p3a4,p3a5,p3a6,p3a7
rigid_group_p4 = p4a0,p4a1,p4a2,p4a3,p4a4,p4a5,p4a6,p4a7
rigid_group_p5 = p5a0,p5a1,p5a2,p5a3,p5a4,p5a5,p5a6,p5a7
rigid_group_p6 = p6a0,p6a1,p6a2,p6a3,p6a4,p6a5,p6a6,p6a7
rigid_group_p7 = p7a0,p7a1,p7a2,p7a3,p7a4,p7a5,p7a6,p7a7
rigid_group_p8 = p8a0,p8a1,p8a2,p8a3,p8a4,p8a5,p8a6,p8a7
rigid_group_p9 = p9a0,p9a1,p9a2,p9a3,p9a4,p9a5,p9a6,p9a7
rigid_group_p10 = p10a0,p10a1,p10a2,p10a3,p10a4,p10a5,p10a6,p10a7
rigid_group_p11 = p11a0,p11a1,p11a2,p11a3,p11a4,p11a5,p11a6,p11a7
rigid_group_p12 = p12a0,p12a1,p12a2,p12a3,p12a4,p12a5,p12a6,p12a7
rigid_group_p13 = p13a0,p13a1,p13a2,p13a3,p13a4,p13a5,p13a6,p13a7
rigid_group_p14 = p14a0,p14a1,p14a2,p14a3,p14a4,p14a5,p14a6,p14a7
rigid_group_p15 = p15a0,p15a1,p15a2,p15a3,p15a4,p15a5,p15a6,p15a7
rigid_group_collection_modules = p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15

rigid_group_q0 = p0a0,p0a1,p0a2,p0a3,p0a4,p0a5,p0a6,p0a7,p1a0,p1a1,p1a2,p1a3,p1a4,p1a5,p1a6,p1a7,p2a0,p2a1,p2a2,p2a3,p2a4,p2a5,p2a6,p2a7,p3a0,p3a1,p3a2,p3a3,p3a4,p3a5,p3a6,p3a7
rigid_group_q1 = p4a0,p4a1,p4a2,p4a3,p4a4,p4a5,p4a6,p4a7,p5a0,p5a1,p5a2,p5a3,p5a4,p5a5,p5a6,p5a7,p6a0,p6a1,p6a2,p6a3,p6a4,p6a5,p6a6,p6a7,p7a0,p7a1,p7a2,p7a3,p7a4,p7a5,p7a6,p7a7
rigid_group_q2 = p8a0,p8a1,p8a2,p8a3,p8a4,p8a5,p8a6,p8a7,p9a0,p9a1,p9a2,p9a3,p9a4,p9a5,p9a6,p9a7,p10a0,p10a1,p10a2,p10a3,p10a4,p10a5,p10a6,p10a7,p11a0,p11a1,p11a2,p11a3,p11a4,p11a5,p11a6,p11a7
rigid_group_q3 = p12a0,p12a1,p12a2,p12a3,p12a4,p12a5,p12a6,p12a7,p13a0,p13a1,p13a2,p13a3,p13a4,p13a5,p13a6,p13a7,p14a0,p14a1,p14a2,p14a3,p14a4,p14a5,p14a6,p14a7,p15a0,p15a1,p15a2,p15a3,p15a4,p15a5,p15a6,p15a7
rigid_group_collection_quads = q0,q1,q2,q3

p0a0/dim1 = 0
p0a0/dim2 = ss
p0a0/dim3 = fs
p0a0/min_fs = 0
p0a0/min_ss = 0
p0a0/max_fs = 127
p0a0/max_ss = 63
p0a0/fs = +0.003363999999999999x -0.999995y
p0a0/ss = +0.999995x +0.003363999999999999y
p0a0/corner_x = -528.694
p0a0/corner_y = 637.243
p0a0/coffset = 0.005007

p0a1/dim1 = 0
p0a1/dim2 = ss
p0a1/dim3 = fs
p0a1/min_fs = 0
p0a1/min_ss = 64
p0a1/max_fs = 127
p0a1/max_ss = 127
p0a1/fs = +0.003363999999999999x -0.999995y
p0a1/ss = +0.999995x +0.003363999999999999y
p0a1/corner_x = -462.694
p0a1/corner_y = 637.4649999999999
p0a1/coffset = 0.005007

p0a2/dim1 = 0
p0a2/dim2 = ss
p0a2/dim3 = fs
p0a2/min_fs = 0
p0a2/min_ss = 128
p0a2/max_fs = 127
p0a2/max_ss = 191
p0a2/fs = +0.003363999999999999x -0.999995y
p0a2/ss = +0.999995x +0.003363999999999999y
p0a2/corner_x = -396.6939999999999
p0a2/corner_y = 637.6869999999998
p0a2/coffset = 0.005007

p0a3/dim1 = 0
p0a3/dim2 = ss
p0a3/dim3 = fs
p0a3/min_fs = 0
p0a3/min_ss = 192
p0a3/max_fs = 127
p0a3/max_ss = 255
p0a3/fs = +0.003363999999999999x -0.999995y
p0a3/ss = +0.999995x +0.003363999999999999y
p0a3/corner_x = -330.695
p0a3/corner_y = 637.909
p0a3/coffset = 0.005007

p0a4/dim1 = 0
p0a4/dim2 = ss
p0a4/dim3 = fs
p0a4/min_fs = 0
p0a4/min_ss = 256
p0a4/max_fs = 127
p0a4/max_ss = 319
p0a4/fs = +0.003363999999999999x -0.999995y
p0a4/ss = +0.999995x +0.003363999999999999y
p0a4/corner_x = -264.695
p0a4/corner_y = 638.132
p0a4/coffset = 0.005007

p0a5/dim1 = 0
p0a5/dim2 = ss
p0a5/dim3 = fs
p0a5/min_fs = 0
p0a5/min_ss = 320
p0a5/max_fs = 127
p0a5/max_ss = 383
p0a5/fs = +0.003363999999999999x -0.999995y
p0a5/ss = +0.999995x +0.003363999999999999y
p0a5/corner_x = -198.69499999999994
p0a5/corner_y = 638.353
p0a5/coffset = 0.005007

p0a6/dim1 = 0
p0a6/dim2 = ss
p0a6/dim3 = fs
p0a6/min_fs = 0
p0a6/min_ss = 384
p0a6/max_fs = 127
p0a6/max_ss = 447
p0a6/fs = +0.003363999999999999x -0.999995y
p0a6/ss = +0.999995x +0.003363999999999999y
p0a6/corner_x = -132.696
p0a6/corner_y = 638.575
p0a6/coffset = 0.005007

p0a7/dim1 = 0
p0a7/dim2 = ss
p0a7/dim3 = fs
p0a7/min_fs = 0
p0a7/min_ss = 448
p0a7/max_fs = 127
p0a7/max_ss = 511
p0a7/fs = +0.003363999999999999x -0.999995y
p0a7/ss = +0.999995x +0.003363999999999999y
p0a7/corner_x = -66.696
p0a7/corner_y = 638.7969999999999
p0a7/coffset = 0.005007

p1a0/dim1 = 1
p1a0/dim2 = ss
p1a0/dim3 = fs
p1a0/min_fs = 0
p1a0/min_ss = 0
p1a0/max_fs = 127
p1a0/max_ss = 63
p1a0/fs = +0.006471x -0.99998y
p1a0/ss = +0.99998x +0.006471y
p1a0/corner_x = -528.1859999999999
p1a0/corner_y = 478.614
p1a0/coffset = -0.000282

p1a1/dim1 = 1
p1a1/dim2 = ss
p1a1/dim3 = fs
p1a1/min_fs = 0
p1a1/min_ss = 64
p1a1/max_fs = 127
p1a1/max_ss = 127
p1a1/fs = +0.006471x -0.99998y
p1a1/ss = +0.99998x +0.006471y
p1a1/corner_x = -462.187
p1a1/corner_y = 479.04
p1a1/coffset = -0.000282

p1a2/dim1 = 1
p1a2/dim2 = ss
p1a2/dim3 = fs
p1a2/min_fs = 0
p1a2/min_ss = 128
p1a2/max_fs = 127
p1a2/max_ss = 191
p1a2/fs = +0.006471x -0.99998y
p1a2/ss = +0.99998x +0.006471y
p1a2/corner_x = -396.18799999999993
p1a2/corner_y = 479.468
p1a2/coffset = -0.000282

p1a3/dim1 = 1
p1a3/dim2 = ss
p1a3/dim3 = fs
p1a3/min_fs = 0
p1a3/min_ss = 192
p1a3/max_fs = 127
p1a3/max_ss = 255
p1a3/fs = +0.006471x -0.99998y
p1a3/ss = +0.99998x +0.006471y
p1a3/corner_x = -330.19
p1a3/corner_y = 479.894
p1a3/coffset = -0.000282

p1a4/dim1 = 1
p1a4/dim2 = ss
p1a4/dim3 = fs
p1a4/min_fs = 0
p1a4/min_ss = 256
p1a4/max_fs = 127
p1a4/max_ss = 319
p1a4/fs = +0.006471x -0.99998y
p1a4/ss = +0.99998x +0.006471y
p1a4/corner_x = -264.191
p1a4/corner_y = 480.32199999999995
p1a4/coffset = -0.000282

p1a5/dim1 = 1
p1a5/dim2 = ss
p1a5/dim3 = fs
p1a5/min_fs = 0
p1a5/min_ss = 320
p1a5/max_fs = 127
p1a5/max_ss = 383
p1a5/fs = +0.006471x -0.99998y
p1a5/ss = +0.99998x +0.006471y
p1a5/corner_x = -198.192
p1a5/corner_y = 480.748
p1a5/coffset = -0.000282

p1a6/dim1 = 1
p1a6/dim2 = ss
p1a6/dim3 = fs
p1a6/min_fs = 0
p1a6/min_ss = 384
p1a6/max_fs = 127
p1a6/max_ss = 447
p1a6/fs = +0.006471x -0.99998y
p1a6/ss = +0.99998x +0.006471y
p1a6/corner_x = -132.194
p1a6/corner_y = 481.17599999999993
p1a6/coffset = -0.000282

p1a7/dim1 = 1
p1a7/dim2 = ss
p1a7/dim3 = fs
p1a7/min_fs = 0
p1a7/min_ss = 448
p1a7/max_fs = 127
p1a7/max_ss = 511
p1a7/fs = +0.006471x -0.99998y
p1a7/ss = +0.99998x +0.006471y
p1a7/corner_x = -66.1956
p1a7/corner_y = 481.602
p1a7/coffset = -0.000282

p2a0/dim1 = 2
p2a0/dim2 = ss
p2a0/dim3 = fs
p2a0/min_fs = 0
p2a0/min_ss = 0
p2a0/max_fs = 127
p2a0/max_ss = 63
p2a0/fs = +0.007873x -0.999969y
p2a0/ss = +0.999969x +0.007873y
p2a0/corner_x = -527.542
p2a0/corner_y = 322.1139999999999
p2a0/coffset = -0.000455

p2a1/dim1 = 2
p2a1/dim2 = ss
p2a1/dim3 = fs
p2a1/min_fs = 0
p2a1/min_ss = 64
p2a1/max_fs = 127
p2a1/max_ss = 127
p2a1/fs = +0.007873x -0.999969y
p2a1/ss = +0.999969x +0.007873y
p2a1/corner_x = -461.544
p2a1/corner_y = 322.634
p2a1/coffset = -0.000455

p2a2/dim1 = 2
p2a2/dim2 = ss
p2a2/dim3 = fs
p2a2/min_fs = 0
p2a2/min_ss = 128
p2a2/max_fs = 127
p2a2/max_ss = 191
p2a2/fs = +0.007873x -0.999969y
p2a2/ss = +0.999969x +0.007873y
p2a2/corner_x = -395.546
p2a2/corner_y = 323.153
p2a2/coffset = -0.000455

p2a3/dim1 = 2
p2a3/dim2 = ss
p2a3/dim3 = fs
p2a3/min_fs = 0
p2a3/min_ss = 192
p2a3/max_fs = 127
p2a3/max_ss = 255
p2a3/fs = +0.007873x -0.999969y
p2a3/ss = +0.999969x +0.007873y
p2a3/corner_x = -329.5489999999999
p2a3/corner_y = 323.672
p2a3/coffset = -0.000455

p2a4/dim1 = 2
p2a4/dim2 = ss
p2a4/dim3 = fs
p2a4/min_fs = 0
p2a4/min_ss = 256
p2a4/max_fs = 127
p2a4/max_ss = 319
p2a4/fs = +0.007873x -0.999969y
p2a4/ss = +0.999969x +0.007873y
p2a4/corner_x = -263.55
p2a4/corner_y = 324.193
p2a4/coffset = -0.000455

p2a5/dim1 = 2
p2a5/dim2 = ss
p2a5/dim3 = fs
p2a5/min_fs = 0
p2a5/min_ss = 320
p2a5/max_fs = 127
p2a5/max_ss = 383
p2a5/fs = +0.007873x -0.999969y
p2a5/ss = +0.999969x +0.007873y
p2a5/corner_x = -197.552
p2a5/corner_y = 324.711
p2a5/coffset = -0.000455

p2a6/dim1 = 2
p2a6/dim2 = ss
p2a6/dim3 = fs
p2a6/min_fs = 0
p2a6/min_ss = 384
p2a6/max_fs = 127
p2a6/max_ss = 447
p2a6/fs = +0.007873x -0.999969y
p2a6/ss = +0.999969x +0.007873y
p2a6/corner_x = -131.55399999999997
p2a6/corner_y = 325.232
p2a6/coffset = -0.000455

p2a7/dim1 = 2
p2a7/dim2 = ss
p2a7/dim3 = fs
p2a7/min_fs = 0
p2a7/min_ss = 448
p2a7/max_fs = 127
p2a7/max_ss = 511
p2a7/fs = +0.007873x -0.999969y
p2a7/ss = +0.999969x +0.007873y
p2a7/corner_x = -65.5568
p2a7/corner_y = 325.75099999999986
p2a7/coffset = -0.000455

p3a0/dim1 = 3
p3a0/dim2 = ss
p3a0/dim3 = fs
p3a0/min_fs = 0
p3a0/min_ss = 0
p3a0/max_fs = 127
p3a0/max_ss = 63
p3a0/fs = +0.006664x -0.9999779999999998y
p3a0/ss = +0.9999779999999998x +0.006664y
p3a0/corner_x = -527.228
p3a0/corner_y = 165.36099999999993
p3a0/coffset = -0.000523

p3a1/dim1 = 3
p3a1/dim2 = ss
p3a1/dim3 = fs
p3a1/min_fs = 0
p3a1/min_ss = 64
p3a1/max_fs = 127
p3a1/max_ss = 127
p3a1/fs = +0.006664x -0.9999779999999998y
p3a1/ss = +0.9999779999999998x +0.006664y
p3a1/corner_x = -461.229
p3a1/corner_y = 165.8
p3a1/coffset = -0.000523

p3a2/dim1 = 3
p3a2/dim2 = ss
p3a2/dim3 = fs
p3a2/min_fs = 0
p3a2/min_ss = 128
p3a2/max_fs = 127
p3a2/max_ss = 191
p3a2/fs = +0.006664x -0.9999779999999998y
p3a2/ss = +0.9999779999999998x +0.006664y
p3a2/corner_x = -395.231
p3a2/corner_y = 166.23999999999995
p3a2/coffset = -0.000523

p3a3/dim1 = 3
p3a3/dim2 = ss
p3a3/dim3 = fs
p3a3/min_fs = 0
p3a3/min_ss = 192
p3a3/max_fs = 127
p3a3/max_ss = 255
p3a3/fs = +0.006664x -0.9999779999999998y
p3a3/ss = +0.9999779999999998x +0.006664y
p3a3/corner_x = -329.23199999999997
p3a3/corner_y = 166.67899999999995
p3a3/coffset = -0.000523

p3a4/dim1 = 3
p3a4/dim2 = ss
p3a4/dim3 = fs
p3a4/min_fs = 0
p3a4/min_ss = 256
p3a4/max_fs = 127
p3a4/max_ss = 319
p3a4/fs = +0.006664x -0.9999779999999998y
p3a4/ss = +0.9999779999999998x +0.006664y
p3a4/corner_x = -263.23499999999996
p3a4/corner_y = 167.12
p3a4/coffset = -0.000523

p3a5/dim1 = 3
p3a5/dim2 = ss
p3a5/dim3 = fs
p3a5/min_fs = 0
p3a5/min_ss = 320
p3a5/max_fs = 127
p3a5/max_ss = 383
p3a5/fs = +0.006664x -0.9999779999999998y
p3a5/ss = +0.9999779999999998x +0.006664y
p3a5/corner_x = -197.235
p3a5/corner_y = 167.559
p3a5/coffset = -0.000523

p3a6/dim1 = 3
p3a6/dim2 = ss
p3a6/dim3 = fs
p3a6/min_fs = 0
p3a6/min_ss = 384
p3a6/max_fs = 127
p3a6/max_ss = 447
p3a6/fs = +0.006664x -0.9999779999999998y
p3a6/ss = +0.9999779999999998x +0.006664y
p3a6/corner_x = -131.237
p3a6/corner_y = 167.999
p3a6/coffset = -0.000523

p3a7/dim1 = 3
p3a7/dim2 = ss
p3a7/dim3 = fs
p3a7/min_fs = 0
p3a7/min_ss = 448
p3a7/max_fs = 127
p3a7/max_ss = 511
p3a7/fs = +0.006664x -0.9999779999999998y
p3a7/ss = +0.9999779999999998x +0.006664y
p3a7/corner_x = -65.2385
p3a7/corner_y = 168.43899999999994
p3a7/coffset = -0.000523

p4a0/dim1 = 4
p4a0/dim2 = ss
p4a0/dim3 = fs
p4a0/min_fs = 0
p4a0/min_ss = 0
p4a0/max_fs = 127
p4a0/max_ss = 63
p4a0/fs = +0.004449x -0.999991y
p4a0/ss = +0.999991x +0.004449y
p4a0/corner_x = -539.498
p4a0/corner_y = -21.107
p4a0/coffset = -0.000607

p4a1/dim1 = 4
p4a1/dim2 = ss
p4a1/dim3 = fs
p4a1/min_fs = 0
p4a1/min_ss = 64
p4a1/max_fs = 127
p4a1/max_ss = 127
p4a1/fs = +0.004449x -0.999991y
p4a1/ss = +0.999991x +0.004449y
p4a1/corner_x = -473.4979999999999
p4a1/corner_y = -20.8134
p4a1/coffset = -0.000607

p4a2/dim1 = 4
p4a2/dim2 = ss
p4a2/dim3 = fs
p4a2/min_fs = 0
p4a2/min_ss = 128
p4a2/max_fs = 127
p4a2/max_ss = 191
p4a2/fs = +0.004449x -0.999991y
p4a2/ss = +0.999991x +0.004449y
p4a2/corner_x = -407.5
p4a2/corner_y = -20.5197
p4a2/coffset = -0.000607

p4a3/dim1 = 4
p4a3/dim2 = ss
p4a3/dim3 = fs
p4a3/min_fs = 0
p4a3/min_ss = 192
p4a3/max_fs = 127
p4a3/max_ss = 255
p4a3/fs = +0.004449x -0.999991y
p4a3/ss = +0.999991x +0.004449y
p4a3/corner_x = -341.5
p4a3/corner_y = -20.22609999999999
p4a3/coffset = -0.000607

p4a4/dim1 = 4
p4a4/dim2 = ss
p4a4/dim3 = fs
p4a4/min_fs = 0
p4a4/min_ss = 256
p4a4/max_fs = 127
p4a4/max_ss = 319
p4a4/fs = +0.004449x -0.999991y
p4a4/ss = +0.999991x +0.004449y
p4a4/corner_x = -275.5
p4a4/corner_y = -19.9324
p4a4/coffset = -0.000607

p4a5/dim1 = 4
p4a5/dim2 = ss
p4a5/dim3 = fs
p4a5/min_fs = 0
p4a5/min_ss = 320
p4a5/max_fs = 127
p4a5/max_ss = 383
p4a5/fs = +0.004449x -0.999991y
p4a5/ss = +0.999991x +0.004449y
p4a5/corner_x = -209.5
p4a5/corner_y = -19.638799999999993
p4a5/coffset = -0.000607

p4a6/dim1 = 4
p4a6/dim2 = ss
p4a6/dim3 = fs
p4a6/min_fs = 0
p4a6/min_ss = 384
p4a6/max_fs = 127
p4a6/max_ss = 447
p4a6/fs = +0.004449x -0.999991y
p4a6/ss = +0.999991x +0.004449y
p4a6/corner_x = -143.501
p4a6/corner_y = -19.3451
p4a6/coffset = -0.000607

p4a7/dim1 = 4
p4a7/dim2 = ss
p4a7/dim3 = fs
p4a7/min_fs = 0
p4a7/min_ss = 448
p4a7/max_fs = 127
p4a7/max_ss = 511
p4a7/fs = +0.004449x -0.999991y
p4a7/ss = +0.999991x +0.004449y
p4a7/corner_x = -77.5025
p4a7/corner_y = -19.051499999999997
p4a7/coffset = -0.000607

p5a0/dim1 = 5
p5a0/dim2 = ss
p5a0/dim3 = fs
p5a0/min_fs = 0
p5a0/min_ss = 0
p5a0/max_fs = 127
p5a0/max_ss = 63
p5a0/fs = +0.004179x -0.999991y
p5a0/ss = +0.999991x +0.004179y
p5a0/corner_x = -539.582
p5a0/corner_y = -177.87399999999997
p5a0/coffset = -0.000563

p5a1/dim1 = 5
p5a1/dim2 = ss
p5a1/dim3 = fs
p5a1/min_fs = 0
p5a1/min_ss = 64
p5a1/max_fs = 127
p5a1/max_ss = 127
p5a1/fs = +0.004179x -0.999991y
p5a1/ss = +0.999991x +0.004179y
p5a1/corner_x = -473.583
p5a1/corner_y = -177.598
p5a1/coffset = -0.000563

p5a2/dim1 = 5
p5a2/dim2 = ss
p5a2/dim3 = fs
p5a2/min_fs = 0
p5a2/min_ss = 128
p5a2/max_fs = 127
p5a2/max_ss = 191
p5a2/fs = +0.004179x -0.999991y
p5a2/ss = +0.999991x +0.004179y
p5a2/corner_x = -407.583
p5a2/corner_y = -177.322
p5a2/coffset = -0.000563

p5a3/dim1 = 5
p5a3/dim2 = ss
p5a3/dim3 = fs
p5a3/min_fs = 0
p5a3/min_ss = 192
p5a3/max_fs = 127
p5a3/max_ss = 255
p5a3/fs = +0.004179x -0.999991y
p5a3/ss = +0.999991x +0.004179y
p5a3/corner_x = -341.5839999999999
p5a3/corner_y = -177.046
p5a3/coffset = -0.000563

p5a4/dim1 = 5
p5a4/dim2 = ss
p5a4/dim3 = fs
p5a4/min_fs = 0
p5a4/min_ss = 256
p5a4/max_fs = 127
p5a4/max_ss = 319
p5a4/fs = +0.004179x -0.999991y
p5a4/ss = +0.999991x +0.004179y
p5a4/corner_x = -275.584
p5a4/corner_y = -176.77
p5a4/coffset = -0.000563

p5a5/dim1 = 5
p5a5/dim2 = ss
p5a5/dim3 = fs
p5a5/min_fs = 0
p5a5/min_ss = 320
p5a5/max_fs = 127
p5a5/max_ss = 383
p5a5/fs = +0.004179x -0.999991y
p5a5/ss = +0.999991x +0.004179y
p5a5/corner_x = -209.585
p5a5/corner_y = -176.495
p5a5/coffset = -0.000563

p5a6/dim1 = 5
p5a6/dim2 = ss
p5a6/dim3 = fs
p5a6/min_fs = 0
p5a6/min_ss = 384
p5a6/max_fs = 127
p5a6/max_ss = 447
p5a6/fs = +0.004179x -0.999991y
p5a6/ss = +0.999991x +0.004179y
p5a6/corner_x = -143.586
p5a6/corner_y = -176.219
p5a6/coffset = -0.000563

p5a7/dim1 = 5
p5a7/dim2 = ss
p5a7/dim3 = fs
p5a7/min_fs = 0
p5a7/min_ss = 448
p5a7/max_fs = 127
p5a7/max_ss = 511
p5a7/fs = +0.004179x -0.999991y
p5a7/ss = +0.999991x +0.004179y
p5a7/corner_x = -77.58629999999998
p5a7/corner_y = -175.94299999999996
p5a7/coffset = -0.000563

p6a0/dim1 = 6
p6a0/dim2 = ss
p6a0/dim3 = fs
p6a0/min_fs = 0
p6a0/min_ss = 0
p6a0/max_fs = 127
p6a0/max_ss = 63
p6a0/fs = +0.002516x -0.999996y
p6a0/ss = +0.999996x +0.002516y
p6a0/corner_x = -539.185
p6a0/corner_y = -334.222
p6a0/coffset = -0.000497

p6a1/dim1 = 6
p6a1/dim2 = ss
p6a1/dim3 = fs
p6a1/min_fs = 0
p6a1/min_ss = 64
p6a1/max_fs = 127
p6a1/max_ss = 127
p6a1/fs = +0.002516x -0.999996y
p6a1/ss = +0.999996x +0.002516y
p6a1/corner_x = -473.1849999999999
p6a1/corner_y = -334.056
p6a1/coffset = -0.000497

p6a2/dim1 = 6
p6a2/dim2 = ss
p6a2/dim3 = fs
p6a2/min_fs = 0
p6a2/min_ss = 128
p6a2/max_fs = 127
p6a2/max_ss = 191
p6a2/fs = +0.002516x -0.999996y
p6a2/ss = +0.999996x +0.002516y
p6a2/corner_x = -407.18499999999995
p6a2/corner_y = -333.89
p6a2/coffset = -0.000497

p6a3/dim1 = 6
p6a3/dim2 = ss
p6a3/dim3 = fs
p6a3/min_fs = 0
p6a3/min_ss = 192
p6a3/max_fs = 127
p6a3/max_ss = 255
p6a3/fs = +0.002516x -0.999996y
p6a3/ss = +0.999996x +0.002516y
p6a3/corner_x = -341.185
p6a3/corner_y = -333.72499999999997
p6a3/coffset = -0.000497

p6a4/dim1 = 6
p6a4/dim2 = ss
p6a4/dim3 = fs
p6a4/min_fs = 0
p6a4/min_ss = 256
p6a4/max_fs = 127
p6a4/max_ss = 319
p6a4/fs = +0.002516x -0.999996y
p6a4/ss = +0.999996x +0.002516y
p6a4/corner_x = -275.186
p6a4/corner_y = -333.5589999999999
p6a4/coffset = -0.000497

p6a5/dim1 = 6
p6a5/dim2 = ss
p6a5/dim3 = fs
p6a5/min_fs = 0
p6a5/min_ss = 320
p6a5/max_fs = 127
p6a5/max_ss = 383
p6a5/fs = +0.002516x -0.999996y
p6a5/ss = +0.999996x +0.002516y
p6a5/corner_x = -209.18599999999998
p6a5/corner_y = -333.392
p6a5/coffset = -0.000497

p6a6/dim1 = 6
p6a6/dim2 = ss
p6a6/dim3 = fs
p6a6/min_fs = 0
p6a6/min_ss = 384
p6a6/max_fs = 127
p6a6/max_ss = 447
p6a6/fs = +0.002516x -0.999996y
p6a6/ss = +0.999996x +0.002516y
p6a6/corner_x = -143.186
p6a6/corner_y = -333.226
p6a6/coffset = -0.000497

p6a7/dim1 = 6
p6a7/dim2 = ss
p6a7/dim3 = fs
p6a7/min_fs = 0
p6a7/min_ss = 448
p6a7/max_fs = 127
p6a7/max_ss = 511
p6a7/fs = +0.002516x -0.999996y
p6a7/ss = +0.999996x +0.002516y
p6a7/corner_x = -77.1862
p6a7/corner_y = -333.06
p6a7/coffset = -0.000497

p7a0/dim1 = 7
p7a0/dim2 = ss
p7a0/dim3 = fs
p7a0/min_fs = 0
p7a0/min_ss = 0
p7a0/max_fs = 127
p7a0/max_ss = 63
p7a0/fs = +0.003771x -0.9999929999999999y
p7a0/ss = +0.9999929999999999x +0.003771y
p7a0/corner_x = -538.734
p7a0/corner_y = -490.62299999999993
p7a0/coffset = -0.000436

p7a1/dim1 = 7
p7a1/dim2 = ss
p7a1/dim3 = fs
p7a1/min_fs = 0
p7a1/min_ss = 64
p7a1/max_fs = 127
p7a1/max_ss = 127
p7a1/fs = +0.003771x -0.9999929999999999y
p7a1/ss = +0.9999929999999999x +0.003771y
p7a1/corner_x = -472.73299999999995
p7a1/corner_y = -490.373
p7a1/coffset = -0.000436

p7a2/dim1 = 7
p7a2/dim2 = ss
p7a2/dim3 = fs
p7a2/min_fs = 0
p7a2/min_ss = 128
p7a2/max_fs = 127
p7a2/max_ss = 191
p7a2/fs = +0.003771x -0.9999929999999999y
p7a2/ss = +0.9999929999999999x +0.003771y
p7a2/corner_x = -406.734
p7a2/corner_y = -490.125
p7a2/coffset = -0.000436

p7a3/dim1 = 7
p7a3/dim2 = ss
p7a3/dim3 = fs
p7a3/min_fs = 0
p7a3/min_ss = 192
p7a3/max_fs = 127
p7a3/max_ss = 255
p7a3/fs = +0.003771x -0.9999929999999999y
p7a3/ss = +0.9999929999999999x +0.003771y
p7a3/corner_x = -340.73399999999987
p7a3/corner_y = -489.875
p7a3/coffset = -0.000436

p7a4/dim1 = 7
p7a4/dim2 = ss
p7a4/dim3 = fs
p7a4/min_fs = 0
p7a4/min_ss = 256
p7a4/max_fs = 127
p7a4/max_ss = 319
p7a4/fs = +0.003771x -0.9999929999999999y
p7a4/ss = +0.9999929999999999x +0.003771y
p7a4/corner_x = -274.734
p7a4/corner_y = -489.628
p7a4/coffset = -0.000436

p7a5/dim1 = 7
p7a5/dim2 = ss
p7a5/dim3 = fs
p7a5/min_fs = 0
p7a5/min_ss = 320
p7a5/max_fs = 127
p7a5/max_ss = 383
p7a5/fs = +0.003771x -0.9999929999999999y
p7a5/ss = +0.9999929999999999x +0.003771y
p7a5/corner_x = -208.73499999999999
p7a5/corner_y = -489.3779999999999
p7a5/coffset = -0.000436

p7a6/dim1 = 7
p7a6/dim2 = ss
p7a6/dim3 = fs
p7a6/min_fs = 0
p7a6/min_ss = 384
p7a6/max_fs = 127
p7a6/max_ss = 447
p7a6/fs = +0.003771x -0.9999929999999999y
p7a6/ss = +0.9999929999999999x +0.003771y
p7a6/corner_x = -142.736
p7a6/corner_y = -489.129
p7a6/coffset = -0.000436

p7a7/dim1 = 7
p7a7/dim2 = ss
p7a7/dim3 = fs
p7a7/min_fs = 0
p7a7/min_ss = 448
p7a7/max_fs = 127
p7a7/max_ss = 511
p7a7/fs = +0.003771x -0.9999929999999999y
p7a7/ss = +0.9999929999999999x +0.003771y
p7a7/corner_x = -76.7367
p7a7/corner_y = -488.881
p7a7/coffset = -0.000436

p8a0/dim1 = 8
p8a0/dim2 = ss
p8a0/dim3 = fs
p8a0/min_fs = 0
p8a0/min_ss = 0
p8a0/max_fs = 127
p8a0/max_ss = 63
p8a0/fs = -0.0028669999999999998x +0.999996y
p8a0/ss = -0.999996x -0.0028669999999999998y
p8a0/corner_x = 539.069
p8a0/corner_y = -153.634
p8a0/coffset = -0.000685

p8a1/dim1 = 8
p8a1/dim2 = ss
p8a1/dim3 = fs
p8a1/min_fs = 0
p8a1/min_ss = 64
p8a1/max_fs = 127
p8a1/max_ss = 127
p8a1/fs = -0.0028669999999999998x +0.999996y
p8a1/ss = -0.999996x -0.0028669999999999998y
p8a1/corner_x = 473.06899999999996
p8a1/corner_y = -153.824
p8a1/coffset = -0.000685

p8a2/dim1 = 8
p8a2/dim2 = ss
p8a2/dim3 = fs
p8a2/min_fs = 0
p8a2/min_ss = 128
p8a2/max_fs = 127
p8a2/max_ss = 191
p8a2/fs = -0.0028669999999999998x +0.999996y
p8a2/ss = -0.999996x -0.0028669999999999998y
p8a2/corner_x = 407.069
p8a2/corner_y = -154.01299999999998
p8a2/coffset = -0.000685

p8a3/dim1 = 8
p8a3/dim2 = ss
p8a3/dim3 = fs
p8a3/min_fs = 0
p8a3/min_ss = 192
p8a3/max_fs = 127
p8a3/max_ss = 255
p8a3/fs = -0.0028669999999999998x +0.999996y
p8a3/ss = -0.999996x -0.0028669999999999998y
p8a3/corner_x = 341.069
p8a3/corner_y = -154.202
p8a3/coffset = -0.000685

p8a4/dim1 = 8
p8a4/dim2 = ss
p8a4/dim3 = fs
p8a4/min_fs = 0
p8a4/min_ss = 256
p8a4/max_fs = 127
p8a4/max_ss = 319
p8a4/fs = -0.0028669999999999998x +0.999996y
p8a4/ss = -0.999996x -0.0028669999999999998y
p8a4/corner_x = 275.07
p8a4/corner_y = -154.392
p8a4/coffset = -0.000685

p8a5/dim1 = 8
p8a5/dim2 = ss
p8a5/dim3 = fs
p8a5/min_fs = 0
p8a5/min_ss = 320
p8a5/max_fs = 127
p8a5/max_ss = 383
p8a5/fs = -0.0028669999999999998x +0.999996y
p8a5/ss = -0.999996x -0.0028669999999999998y
p8a5/corner_x = 209.06999999999994
p8a5/corner_y = -154.58
p8a5/coffset = -0.000685

p8a6/dim1 = 8
p8a6/dim2 = ss
p8a6/dim3 = fs
p8a6/min_fs = 0
p8a6/min_ss = 384
p8a6/max_fs = 127
p8a6/max_ss = 447
p8a6/fs = -0.0028669999999999998x +0.999996y
p8a6/ss = -0.999996x -0.0028669999999999998y
p8a6/corner_x = 143.06999999999996
p8a6/corner_y = -154.77
p8a6/coffset = -0.000685

p8a7/dim1 = 8
p8a7/dim2 = ss
p8a7/dim3 = fs
p8a7/min_fs = 0
p8a7/min_ss = 448
p8a7/max_fs = 127
p8a7/max_ss = 511
p8a7/fs = -0.0028669999999999998x +0.999996y
p8a7/ss = -0.999996x -0.0028669999999999998y
p8a7/corner_x = 77.07049999999998
p8a7/corner_y = -154.958
p8a7/coffset = -0.000685

p9a0/dim1 = 9
p9a0/dim2 = ss
p9a0/dim3 = fs
p9a0/min_fs = 0
p9a0/min_ss = 0
p9a0/max_fs = 127
p9a0/max_ss = 63
p9a0/fs = -0.003058x +0.999995y
p9a0/ss = -0.999995x -0.003058y
p9a0/corner_x = 539.892
p9a0/corner_y = -310.68799999999993
p9a0/coffset = -0.000651

p9a1/dim1 = 9
p9a1/dim2 = ss
p9a1/dim3 = fs
p9a1/min_fs = 0
p9a1/min_ss = 64
p9a1/max_fs = 127
p9a1/max_ss = 127
p9a1/fs = -0.003058x +0.999995y
p9a1/ss = -0.999995x -0.003058y
p9a1/corner_x = 473.89299999999986
p9a1/corner_y = -310.89
p9a1/coffset = -0.000651

p9a2/dim1 = 9
p9a2/dim2 = ss
p9a2/dim3 = fs
p9a2/min_fs = 0
p9a2/min_ss = 128
p9a2/max_fs = 127
p9a2/max_ss = 191
p9a2/fs = -0.003058x +0.999995y
p9a2/ss = -0.999995x -0.003058y
p9a2/corner_x = 407.893
p9a2/corner_y = -311.092
p9a2/coffset = -0.000651

p9a3/dim1 = 9
p9a3/dim2 = ss
p9a3/dim3 = fs
p9a3/min_fs = 0
p9a3/min_ss = 192
p9a3/max_fs = 127
p9a3/max_ss = 255
p9a3/fs = -0.003058x +0.999995y
p9a3/ss = -0.999995x -0.003058y
p9a3/corner_x = 341.893
p9a3/corner_y = -311.294
p9a3/coffset = -0.000651

p9a4/dim1 = 9
p9a4/dim2 = ss
p9a4/dim3 = fs
p9a4/min_fs = 0
p9a4/min_ss = 256
p9a4/max_fs = 127
p9a4/max_ss = 319
p9a4/fs = -0.003058x +0.999995y
p9a4/ss = -0.999995x -0.003058y
p9a4/corner_x = 275.893
p9a4/corner_y = -311.495
p9a4/coffset = -0.000651

p9a5/dim1 = 9
p9a5/dim2 = ss
p9a5/dim3 = fs
p9a5/min_fs = 0
p9a5/min_ss = 320
p9a5/max_fs = 127
p9a5/max_ss = 383
p9a5/fs = -0.003058x +0.999995y
p9a5/ss = -0.999995x -0.003058y
p9a5/corner_x = 209.894
p9a5/corner_y = -311.697
p9a5/coffset = -0.000651

p9a6/dim1 = 9
p9a6/dim2 = ss
p9a6/dim3 = fs
p9a6/min_fs = 0
p9a6/min_ss = 384
p9a6/max_fs = 127
p9a6/max_ss = 447
p9a6/fs = -0.003058x +0.999995y
p9a6/ss = -0.999995x -0.003058y
p9a6/corner_x = 143.894
p9a6/corner_y = -311.898
p9a6/coffset = -0.000651

p9a7/dim1 = 9
p9a7/dim2 = ss
p9a7/dim3 = fs
p9a7/min_fs = 0
p9a7/min_ss = 448
p9a7/max_fs = 127
p9a7/max_ss = 511
p9a7/fs = -0.003058x +0.999995y
p9a7/ss = -0.999995x -0.003058y
p9a7/corner_x = 77.8938
p9a7/corner_y = -312.1
p9a7/coffset = -0.000651

p10a0/dim1 = 10
p10a0/dim2 = ss
p10a0/dim3 = fs
p10a0/min_fs = 0
p10a0/min_ss = 0
p10a0/max_fs = 127
p10a0/max_ss = 63
p10a0/fs = -0.002637x +0.999996y
p10a0/ss = -0.999996x -0.002637y
p10a0/corner_x = 540.021
p10a0/corner_y = -467.2589999999999
p10a0/coffset = -0.000591

p10a1/dim1 = 10
p10a1/dim2 = ss
p10a1/dim3 = fs
p10a1/min_fs = 0
p10a1/min_ss = 64
p10a1/max_fs = 127
p10a1/max_ss = 127
p10a1/fs = -0.002637x +0.999996y
p10a1/ss = -0.999996x -0.002637y
p10a1/corner_x = 474.021
p10a1/corner_y = -467.433
p10a1/coffset = -0.000591

p10a2/dim1 = 10
p10a2/dim2 = ss
p10a2/dim3 = fs
p10a2/min_fs = 0
p10a2/min_ss = 128
p10a2/max_fs = 127
p10a2/max_ss = 191
p10a2/fs = -0.002637x +0.999996y
p10a2/ss = -0.999996x -0.002637y
p10a2/corner_x = 408.021
p10a2/corner_y = -467.6079999999999
p10a2/coffset = -0.000591

p10a3/dim1 = 10
p10a3/dim2 = ss
p10a3/dim3 = fs
p10a3/min_fs = 0
p10a3/min_ss = 192
p10a3/max_fs = 127
p10a3/max_ss = 255
p10a3/fs = -0.002637x +0.999996y
p10a3/ss = -0.999996x -0.002637y
p10a3/corner_x = 342.02099999999996
p10a3/corner_y = -467.781
p10a3/coffset = -0.000591

p10a4/dim1 = 10
p10a4/dim2 = ss
p10a4/dim3 = fs
p10a4/min_fs = 0
p10a4/min_ss = 256
p10a4/max_fs = 127
p10a4/max_ss = 319
p10a4/fs = -0.002637x +0.999996y
p10a4/ss = -0.999996x -0.002637y
p10a4/corner_x = 276.022
p10a4/corner_y = -467.9559999999999
p10a4/coffset = -0.000591

p10a5/dim1 = 10
p10a5/dim2 = ss
p10a5/dim3 = fs
p10a5/min_fs = 0
p10a5/min_ss = 320
p10a5/max_fs = 127
p10a5/max_ss = 383
p10a5/fs = -0.002637x +0.999996y
p10a5/ss = -0.999996x -0.002637y
p10a5/corner_x = 210.02199999999993
p10a5/corner_y = -468.13
p10a5/coffset = -0.000591

p10a6/dim1 = 10
p10a6/dim2 = ss
p10a6/dim3 = fs
p10a6/min_fs = 0
p10a6/min_ss = 384
p10a6/max_fs = 127
p10a6/max_ss = 447
p10a6/fs = -0.002637x +0.999996y
p10a6/ss = -0.999996x -0.002637y
p10a6/corner_x = 144.022
p10a6/corner_y = -468.3029999999999
p10a6/coffset = -0.000591

p10a7/dim1 = 10
p10a7/dim2 = ss
p10a7/dim3 = fs
p10a7/min_fs = 0
p10a7/min_ss = 448
p10a7/max_fs = 127
p10a7/max_ss = 511
p10a7/fs = -0.002637x +0.999996y
p10a7/ss = -0.999996x -0.002637y
p10a7/corner_x = 78.0218
p10a7/corner_y = -468.478
p10a7/coffset = -0.000591

p11a0/dim1 = 11
p11a0/dim2 = ss
p11a0/dim3 = fs
p11a0/min_fs = 0
p11a0/min_ss = 0
p11a0/max_fs = 127
p11a0/max_ss = 63
p11a0/fs = -0.0019629999999999995x +0.999999y
p11a0/ss = -0.999999x -0.0019629999999999995y
p11a0/corner_x = 541.22
p11a0/corner_y = -624.972
p11a0/coffset = -0.000682

p11a1/dim1 = 11
p11a1/dim2 = ss
p11a1/dim3 = fs
p11a1/min_fs = 0
p11a1/min_ss = 64
p11a1/max_fs = 127
p11a1/max_ss = 127
p11a1/fs = -0.0019629999999999995x +0.999999y
p11a1/ss = -0.999999x -0.0019629999999999995y
p11a1/corner_x = 475.2199999999999
p11a1/corner_y = -625.101
p11a1/coffset = -0.000682

p11a2/dim1 = 11
p11a2/dim2 = ss
p11a2/dim3 = fs
p11a2/min_fs = 0
p11a2/min_ss = 128
p11a2/max_fs = 127
p11a2/max_ss = 191
p11a2/fs = -0.0019629999999999995x +0.999999y
p11a2/ss = -0.999999x -0.0019629999999999995y
p11a2/corner_x = 409.221
p11a2/corner_y = -625.231
p11a2/coffset = -0.000682

p11a3/dim1 = 11
p11a3/dim2 = ss
p11a3/dim3 = fs
p11a3/min_fs = 0
p11a3/min_ss = 192
p11a3/max_fs = 127
p11a3/max_ss = 255
p11a3/fs = -0.0019629999999999995x +0.999999y
p11a3/ss = -0.999999x -0.0019629999999999995y
p11a3/corner_x = 343.221
p11a3/corner_y = -625.361
p11a3/coffset = -0.000682

p11a4/dim1 = 11
p11a4/dim2 = ss
p11a4/dim3 = fs
p11a4/min_fs = 0
p11a4/min_ss = 256
p11a4/max_fs = 127
p11a4/max_ss = 319
p11a4/fs = -0.0019629999999999995x +0.999999y
p11a4/ss = -0.999999x -0.0019629999999999995y
p11a4/corner_x = 277.221
p11a4/corner_y = -625.489
p11a4/coffset = -0.000682

p11a5/dim1 = 11
p11a5/dim2 = ss
p11a5/dim3 = fs
p11a5/min_fs = 0
p11a5/min_ss = 320
p11a5/max_fs = 127
p11a5/max_ss = 383
p11a5/fs = -0.0019629999999999995x +0.999999y
p11a5/ss = -0.999999x -0.0019629999999999995y
p11a5/corner_x = 211.221
p11a5/corner_y = -625.619
p11a5/coffset = -0.000682

p11a6/dim1 = 11
p11a6/dim2 = ss
p11a6/dim3 = fs
p11a6/min_fs = 0
p11a6/min_ss = 384
p11a6/max_fs = 127
p11a6/max_ss = 447
p11a6/fs = -0.0019629999999999995x +0.999999y
p11a6/ss = -0.999999x -0.0019629999999999995y
p11a6/corner_x = 145.22099999999998
p11a6/corner_y = -625.749
p11a6/coffset = -0.000682

p11a7/dim1 = 11
p11a7/dim2 = ss
p11a7/dim3 = fs
p11a7/min_fs = 0
p11a7/min_ss = 448
p11a7/max_fs = 127
p11a7/max_ss = 511
p11a7/fs = -0.0019629999999999995x +0.999999y
p11a7/ss = -0.999999x -0.0019629999999999995y
p11a7/corner_x = 79.2203
p11a7/corner_y = -625.8779999999999
p11a7/coffset = -0.000682

p12a0/dim1 = 12
p12a0/dim2 = ss
p12a0/dim3 = fs
p12a0/min_fs = 0
p12a0/min_ss = 0
p12a0/max_fs = 127
p12a0/max_ss = 63
p12a0/fs = -0.002172x +0.9999979999999998y
p12a0/ss = -0.9999979999999998x -0.002172y
p12a0/corner_x = 553.728
p12a0/corner_y = 489.163
p12a0/coffset = 0.002664

p12a1/dim1 = 12
p12a1/dim2 = ss
p12a1/dim3 = fs
p12a1/min_fs = 0
p12a1/min_ss = 64
p12a1/max_fs = 127
p12a1/max_ss = 127
p12a1/fs = -0.002172x +0.9999979999999998y
p12a1/ss = -0.9999979999999998x -0.002172y
p12a1/corner_x = 487.729
p12a1/corner_y = 489.02
p12a1/coffset = 0.002664

p12a2/dim1 = 12
p12a2/dim2 = ss
p12a2/dim3 = fs
p12a2/min_fs = 0
p12a2/min_ss = 128
p12a2/max_fs = 127
p12a2/max_ss = 191
p12a2/fs = -0.002172x +0.9999979999999998y
p12a2/ss = -0.9999979999999998x -0.002172y
p12a2/corner_x = 421.728
p12a2/corner_y = 488.87599999999986
p12a2/coffset = 0.002664

p12a3/dim1 = 12
p12a3/dim2 = ss
p12a3/dim3 = fs
p12a3/min_fs = 0
p12a3/min_ss = 192
p12a3/max_fs = 127
p12a3/max_ss = 255
p12a3/fs = -0.002172x +0.9999979999999998y
p12a3/ss = -0.9999979999999998x -0.002172y
p12a3/corner_x = 355.728
p12a3/corner_y = 488.733
p12a3/coffset = 0.002664

p12a4/dim1 = 12
p12a4/dim2 = ss
p12a4/dim3 = fs
p12a4/min_fs = 0
p12a4/min_ss = 256
p12a4/max_fs = 127
p12a4/max_ss = 319
p12a4/fs = -0.002172x +0.9999979999999998y
p12a4/ss = -0.9999979999999998x -0.002172y
p12a4/corner_x = 289.73
p12a4/corner_y = 488.59
p12a4/coffset = 0.002664

p12a5/dim1 = 12
p12a5/dim2 = ss
p12a5/dim3 = fs
p12a5/min_fs = 0
p12a5/min_ss = 320
p12a5/max_fs = 127
p12a5/max_ss = 383
p12a5/fs = -0.002172x +0.9999979999999998y
p12a5/ss = -0.9999979999999998x -0.002172y
p12a5/corner_x = 223.73
p12a5/corner_y = 488.4459999999999
p12a5/coffset = 0.002664

p12a6/dim1 = 12
p12a6/dim2 = ss
p12a6/dim3 = fs
p12a6/min_fs = 0
p12a6/min_ss = 384
p12a6/max_fs = 127
p12a6/max_ss = 447
p12a6/fs = -0.002172x +0.9999979999999998y
p12a6/ss = -0.9999979999999998x -0.002172y
p12a6/corner_x = 157.729
p12a6/corner_y = 488.303
p12a6/coffset = 0.002664

p12a7/dim1 = 12
p12a7/dim2 = ss
p12a7/dim3 = fs
p12a7/min_fs = 0
p12a7/min_ss = 448
p12a7/max_fs = 127
p12a7/max_ss = 511
p12a7/fs = -0.002172x +0.9999979999999998y
p12a7/ss = -0.9999979999999998x -0.002172y
p12a7/corner_x = 91.72909999999999
p12a7/corner_y = 488.16
p12a7/coffset = 0.002664

p13a0/dim1 = 13
p13a0/dim2 = ss
p13a0/dim3 = fs
p13a0/min_fs = 0
p13a0/min_ss = 0
p13a0/max_fs = 127
p13a0/max_ss = 63
p13a0/fs = -0.006381x +0.99998y
p13a0/ss = -0.99998x -0.006381y
p13a0/corner_x = 552.2679999999999
p13a0/corner_y = 341.5959999999999
p13a0/coffset = -0.000452

p13a1/dim1 = 13
p13a1/dim2 = ss
p13a1/dim3 = fs
p13a1/min_fs = 0
p13a1/min_ss = 64
p13a1/max_fs = 127
p13a1/max_ss = 127
p13a1/fs = -0.006381x +0.99998y
p13a1/ss = -0.99998x -0.006381y
p13a1/corner_x = 486.268
p13a1/corner_y = 341.174
p13a1/coffset = -0.000452

p13a2/dim1 = 13
p13a2/dim2 = ss
p13a2/dim3 = fs
p13a2/min_fs = 0
p13a2/min_ss = 128
p13a2/max_fs = 127
p13a2/max_ss = 191
p13a2/fs = -0.006381x +0.99998y
p13a2/ss = -0.99998x -0.006381y
p13a2/corner_x = 420.2699999999999
p13a2/corner_y = 340.75399999999996
p13a2/coffset = -0.000452

p13a3/dim1 = 13
p13a3/dim2 = ss
p13a3/dim3 = fs
p13a3/min_fs = 0
p13a3/min_ss = 192
p13a3/max_fs = 127
p13a3/max_ss = 255
p13a3/fs = -0.006381x +0.99998y
p13a3/ss = -0.99998x -0.006381y
p13a3/corner_x = 354.271
p13a3/corner_y = 340.3339999999999
p13a3/coffset = -0.000452

p13a4/dim1 = 13
p13a4/dim2 = ss
p13a4/dim3 = fs
p13a4/min_fs = 0
p13a4/min_ss = 256
p13a4/max_fs = 127
p13a4/max_ss = 319
p13a4/fs = -0.006381x +0.99998y
p13a4/ss = -0.99998x -0.006381y
p13a4/corner_x = 288.272
p13a4/corner_y = 339.912
p13a4/coffset = -0.000452

p13a5/dim1 = 13
p13a5/dim2 = ss
p13a5/dim3 = fs
p13a5/min_fs = 0
p13a5/min_ss = 320
p13a5/max_fs = 127
p13a5/max_ss = 383
p13a5/fs = -0.006381x +0.99998y
p13a5/ss = -0.99998x -0.006381y
p13a5/corner_x = 222.274
p13a5/corner_y = 339.49099999999993
p13a5/coffset = -0.000452

p13a6/dim1 = 13
p13a6/dim2 = ss
p13a6/dim3 = fs
p13a6/min_fs = 0
p13a6/min_ss = 384
p13a6/max_fs = 127
p13a6/max_ss = 447
p13a6/fs = -0.006381x +0.99998y
p13a6/ss = -0.99998x -0.006381y
p13a6/corner_x = 156.27499999999998
p13a6/corner_y = 339.069
p13a6/coffset = -0.000452

p13a7/dim1 = 13
p13a7/dim2 = ss
p13a7/dim3 = fs
p13a7/min_fs = 0
p13a7/min_ss = 448
p13a7/max_fs = 127
p13a7/max_ss = 511
p13a7/fs = -0.006381x +0.99998y
p13a7/ss = -0.99998x -0.006381y
p13a7/corner_x = 90.27629999999999
p13a7/corner_y = 338.649
p13a7/coffset = -0.000452

p14a0/dim1 = 14
p14a0/dim2 = ss
p14a0/dim3 = fs
p14a0/min_fs = 0
p14a0/min_ss = 0
p14a0/max_fs = 127
p14a0/max_ss = 63
p14a0/fs = -0.005358x +0.9999849999999999y
p14a0/ss = -0.9999849999999999x -0.005358y
p14a0/corner_x = 552.739
p14a0/corner_y = 184.59299999999993
p14a0/coffset = -0.000568

p14a1/dim1 = 14
p14a1/dim2 = ss
p14a1/dim3 = fs
p14a1/min_fs = 0
p14a1/min_ss = 64
p14a1/max_fs = 127
p14a1/max_ss = 127
p14a1/fs = -0.005358x +0.9999849999999999y
p14a1/ss = -0.9999849999999999x -0.005358y
p14a1/corner_x = 486.739
p14a1/corner_y = 184.23999999999998
p14a1/coffset = -0.000568

p14a2/dim1 = 14
p14a2/dim2 = ss
p14a2/dim3 = fs
p14a2/min_fs = 0
p14a2/min_ss = 128
p14a2/max_fs = 127
p14a2/max_ss = 191
p14a2/fs = -0.005358x +0.9999849999999999y
p14a2/ss = -0.9999849999999999x -0.005358y
p14a2/corner_x = 420.74099999999993
p14a2/corner_y = 183.887
p14a2/coffset = -0.000568

p14a3/dim1 = 14
p14a3/dim2 = ss
p14a3/dim3 = fs
p14a3/min_fs = 0
p14a3/min_ss = 192
p14a3/max_fs = 127
p14a3/max_ss = 255
p14a3/fs = -0.005358x +0.9999849999999999y
p14a3/ss = -0.9999849999999999x -0.005358y
p14a3/corner_x = 354.741
p14a3/corner_y = 183.533
p14a3/coffset = -0.000568

p14a4/dim1 = 14
p14a4/dim2 = ss
p14a4/dim3 = fs
p14a4/min_fs = 0
p14a4/min_ss = 256
p14a4/max_fs = 127
p14a4/max_ss = 319
p14a4/fs = -0.005358x +0.9999849999999999y
p14a4/ss = -0.9999849999999999x -0.005358y
p14a4/corner_x = 288.742
p14a4/corner_y = 183.179
p14a4/coffset = -0.000568

p14a5/dim1 = 14
p14a5/dim2 = ss
p14a5/dim3 = fs
p14a5/min_fs = 0
p14a5/min_ss = 320
p14a5/max_fs = 127
p14a5/max_ss = 383
p14a5/fs = -0.005358x +0.9999849999999999y
p14a5/ss = -0.9999849999999999x -0.005358y
p14a5/corner_x = 222.744
p14a5/corner_y = 182.825
p14a5/coffset = -0.000568

p14a6/dim1 = 14
p14a6/dim2 = ss
p14a6/dim3 = fs
p14a6/min_fs = 0
p14a6/min_ss = 384
p14a6/max_fs = 127
p14a6/max_ss = 447
p14a6/fs = -0.005358x +0.9999849999999999y
p14a6/ss = -0.9999849999999999x -0.005358y
p14a6/corner_x = 156.745
p14a6/corner_y = 182.472
p14a6/coffset = -0.000568

p14a7/dim1 = 14
p14a7/dim2 = ss
p14a7/dim3 = fs
p14a7/min_fs = 0
p14a7/min_ss = 448
p14a7/max_fs = 127
p14a7/max_ss = 511
p14a7/fs = -0.005358x +0.9999849999999999y
p14a7/ss = -0.9999849999999999x -0.005358y
p14a7/corner_x = 90.7451
p14a7/corner_y = 182.118
p14a7/coffset = -0.000568

p15a0/dim1 = 15
p15a0/dim2 = ss
p15a0/dim3 = fs
p15a0/min_fs = 0
p15a0/min_ss = 0
p15a0/max_fs = 127
p15a0/max_ss = 63
p15a0/fs = -0.004081x +0.999992y
p15a0/ss = -0.999992x -0.004081y
p15a0/corner_x = 553.248
p15a0/corner_y = 28.4926
p15a0/coffset = -0.000679

p15a1/dim1 = 15
p15a1/dim2 = ss
p15a1/dim3 = fs
p15a1/min_fs = 0
p15a1/min_ss = 64
p15a1/max_fs = 127
p15a1/max_ss = 127
p15a1/fs = -0.004081x +0.999992y
p15a1/ss = -0.999992x -0.004081y
p15a1/corner_x = 487.24999999999994
p15a1/corner_y = 28.223299999999995
p15a1/coffset = -0.000679

p15a2/dim1 = 15
p15a2/dim2 = ss
p15a2/dim3 = fs
p15a2/min_fs = 0
p15a2/min_ss = 128
p15a2/max_fs = 127
p15a2/max_ss = 191
p15a2/fs = -0.004081x +0.999992y
p15a2/ss = -0.999992x -0.004081y
p15a2/corner_x = 421.25
p15a2/corner_y = 27.953999999999994
p15a2/coffset = -0.000679

p15a3/dim1 = 15
p15a3/dim2 = ss
p15a3/dim3 = fs
p15a3/min_fs = 0
p15a3/min_ss = 192
p15a3/max_fs = 127
p15a3/max_ss = 255
p15a3/fs = -0.004081x +0.999992y
p15a3/ss = -0.999992x -0.004081y
p15a3/corner_x = 355.251
p15a3/corner_y = 27.684799999999992
p15a3/coffset = -0.000679

p15a4/dim1 = 15
p15a4/dim2 = ss
p15a4/dim3 = fs
p15a4/min_fs = 0
p15a4/min_ss = 256
p15a4/max_fs = 127
p15a4/max_ss = 319
p15a4/fs = -0.004081x +0.999992y
p15a4/ss = -0.999992x -0.004081y
p15a4/corner_x = 289.251
p15a4/corner_y = 27.4153
p15a4/coffset = -0.000679

p15a5/dim1 = 15
p15a5/dim2 = ss
p15a5/dim3 = fs
p15a5/min_fs = 0
p15a5/min_ss = 320
p15a5/max_fs = 127
p15a5/max_ss = 383
p15a5/fs = -0.004081x +0.999992y
p15a5/ss = -0.999992x -0.004081y
p15a5/corner_x = 223.25199999999995
p15a5/corner_y = 27.1461
p15a5/coffset = -0.000679

p15a6/dim1 = 15
p15a6/dim2 = ss
p15a6/dim3 = fs
p15a6/min_fs = 0
p15a6/min_ss = 384
p15a6/max_fs = 127
p15a6/max_ss = 447
p15a6/fs = -0.004081x +0.999992y
p15a6/ss = -0.999992x -0.004081y
p15a6/corner_x = 157.25099999999998
p15a6/corner_y = 26.8768
p15a6/coffset = -0.000679

p15a7/dim1 = 15
p15a7/dim2 = ss
p15a7/dim3 = fs
p15a7/min_fs = 0
p15a7/min_ss = 448
p15a7/max_fs = 127
p15a7/max_ss = 511
p15a7/fs = -0.004081x +0.999992y
p15a7/ss = -0.999992x -0.004081y
p15a7/corner_x = 91.2518
p15a7/corner_y = 26.6075
p15a7/coffset = -0.000679