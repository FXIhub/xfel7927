; AGIPD-1M geometry file written by EXtra-geom 1.12.0
; You may need to edit this file to add:
; - data and mask locations in the file
; - mask_good & mask_bad values to interpret the mask
; - adu_per_eV & photon_energy
; - clen (detector distance)
;
; See: http://www.desy.de/~twhite/crystfel/manual-crystfel_geometry.html

;XGEOM MOTORS=4,2
;XGEOM MOTOR_Q1=55.10413360595703,84.70966339111328
;XGEOM MOTOR_Q2=-9.022299766540527,-71.01100158691406
;XGEOM MOTOR_Q3=35.99867248535156,694.8634033203125
;XGEOM MOTOR_Q4=2.0068368911743164,-32.179752349853516

data = /entry_1/instrument_1/detector_1/data ;
dim0 = %
res = 5000.0 ; pixels per metre

; Beam energy in eV
; photon_energy = SET ME

; Camera length, aka detector distance
; clen = SET ME

; Analogue Digital Units per eV
; adu_per_eV = SET ME
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
p0a0/corner_x = -533.1835553588867
p0a0/corner_y = 643.717723815918
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
p0a1/corner_x = -467.18355535888674
p0a1/corner_y = 643.9397238159179
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
p0a2/corner_x = -401.1835553588867
p0a2/corner_y = 644.1617238159178
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
p0a3/corner_x = -335.18455535888677
p0a3/corner_y = 644.383723815918
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
p0a4/corner_x = -269.1845553588867
p0a4/corner_y = 644.6067238159179
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
p0a5/corner_x = -203.18455535888666
p0a5/corner_y = 644.8277238159179
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
p0a6/corner_x = -137.1855553588867
p0a6/corner_y = 645.049723815918
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
p0a7/corner_x = -71.18555535888672
p0a7/corner_y = 645.2717238159179
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
p1a0/corner_x = -532.6755553588866
p1a0/corner_y = 485.08872381591794
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
p1a1/corner_x = -466.67655535888673
p1a1/corner_y = 485.51472381591793
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
p1a2/corner_x = -400.6775553588867
p1a2/corner_y = 485.942723815918
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
p1a3/corner_x = -334.6795553588867
p1a3/corner_y = 486.368723815918
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
p1a4/corner_x = -268.6805553588867
p1a4/corner_y = 486.7967238159179
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
p1a5/corner_x = -202.68155535888673
p1a5/corner_y = 487.22272381591796
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
p1a6/corner_x = -136.6835553588867
p1a6/corner_y = 487.6507238159179
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
p1a7/corner_x = -70.68515535888672
p1a7/corner_y = 488.07672381591794
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
p2a0/corner_x = -532.0315553588867
p2a0/corner_y = 328.5887238159179
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
p2a1/corner_x = -466.0335553588867
p2a1/corner_y = 329.108723815918
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
p2a2/corner_x = -400.03555535888677
p2a2/corner_y = 329.627723815918
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
p2a3/corner_x = -334.0385553588867
p2a3/corner_y = 330.146723815918
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
p2a4/corner_x = -268.0395553588867
p2a4/corner_y = 330.66772381591795
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
p2a5/corner_x = -202.0415553588867
p2a5/corner_y = 331.185723815918
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
p2a6/corner_x = -136.0435553588867
p2a6/corner_y = 331.70672381591794
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
p2a7/corner_x = -70.04635535888671
p2a7/corner_y = 332.22572381591783
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
p3a0/corner_x = -531.7175553588867
p3a0/corner_y = 171.8357238159179
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
p3a1/corner_x = -465.71855535888676
p3a1/corner_y = 172.27472381591795
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
p3a2/corner_x = -399.7205553588867
p3a2/corner_y = 172.71472381591792
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
p3a3/corner_x = -333.72155535888675
p3a3/corner_y = 173.1537238159179
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
p3a4/corner_x = -267.7245553588867
p3a4/corner_y = 173.59472381591797
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
p3a5/corner_x = -201.7245553588867
p3a5/corner_y = 174.03372381591797
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
p3a6/corner_x = -135.7265553588867
p3a6/corner_y = 174.47372381591794
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
p3a7/corner_x = -69.72805535888672
p3a7/corner_y = 174.9137238159179
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
p4a0/corner_x = -546.0045002441407
p4a0/corner_y = 9.68666683959961
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
p4a1/corner_x = -480.0045002441405
p4a1/corner_y = 9.98026683959961
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
p4a2/corner_x = -414.0065002441406
p4a2/corner_y = 10.273966839599609
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
p4a3/corner_x = -348.0065002441406
p4a3/corner_y = 10.567566839599618
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
p4a4/corner_x = -282.0065002441406
p4a4/corner_y = 10.861266839599608
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
p4a5/corner_x = -216.00650024414062
p4a5/corner_y = 11.154866839599617
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
p4a6/corner_x = -150.00750024414063
p4a6/corner_y = 11.448566839599613
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
p4a7/corner_x = -84.00900024414062
p4a7/corner_y = 11.742166839599614
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
p5a0/corner_x = -546.0885002441406
p5a0/corner_y = -147.08033316040036
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
p5a1/corner_x = -480.08950024414065
p5a1/corner_y = -146.8043331604004
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
p5a2/corner_x = -414.08950024414065
p5a2/corner_y = -146.52833316040036
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
p5a3/corner_x = -348.0905002441405
p5a3/corner_y = -146.25233316040038
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
p5a4/corner_x = -282.0905002441406
p5a4/corner_y = -145.9763331604004
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
p5a5/corner_x = -216.09150024414063
p5a5/corner_y = -145.7013331604004
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
p5a6/corner_x = -150.0925002441406
p5a6/corner_y = -145.42533316040036
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
p5a7/corner_x = -84.0928002441406
p5a7/corner_y = -145.14933316040035
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
p6a0/corner_x = -545.6915002441406
p6a0/corner_y = -303.42833316040037
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
p6a1/corner_x = -479.69150024414046
p6a1/corner_y = -303.2623331604004
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
p6a2/corner_x = -413.69150024414057
p6a2/corner_y = -303.0963331604004
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
p6a3/corner_x = -347.6915002441406
p6a3/corner_y = -302.93133316040036
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
p6a4/corner_x = -281.6925002441406
p6a4/corner_y = -302.7653331604003
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
p6a5/corner_x = -215.6925002441406
p6a5/corner_y = -302.5983331604004
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
p6a6/corner_x = -149.69250024414063
p6a6/corner_y = -302.4323331604004
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
p6a7/corner_x = -83.69270024414062
p6a7/corner_y = -302.2663331604004
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
p7a0/corner_x = -545.2405002441407
p7a0/corner_y = -459.8293331604004
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
p7a1/corner_x = -479.2395002441406
p7a1/corner_y = -459.5793331604004
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
p7a2/corner_x = -413.2405002441406
p7a2/corner_y = -459.3313331604004
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
p7a3/corner_x = -347.2405002441405
p7a3/corner_y = -459.08133316040045
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
p7a4/corner_x = -281.2405002441406
p7a4/corner_y = -458.8343331604004
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
p7a5/corner_x = -215.2415002441406
p7a5/corner_y = -458.5843331604003
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
p7a6/corner_x = -149.24250024414061
p7a6/corner_y = -458.3353331604004
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
p7a7/corner_x = -83.24320024414062
p7a7/corner_y = -458.08733316040036
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
p8a0/corner_x = 525.5582578125
p8a0/corner_y = -198.6310245361328
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
p8a1/corner_x = 459.55825781249996
p8a1/corner_y = -198.8210245361328
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
p8a2/corner_x = 393.5582578125
p8a2/corner_y = -199.0100245361328
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
p8a3/corner_x = 327.5582578125
p8a3/corner_y = -199.1990245361328
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
p8a4/corner_x = 261.5592578125
p8a4/corner_y = -199.38902453613278
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
p8a5/corner_x = 195.55925781249994
p8a5/corner_y = -199.57702453613283
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
p8a6/corner_x = 129.55925781249996
p8a6/corner_y = -199.7670245361328
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
p8a7/corner_x = 63.55975781249998
p8a7/corner_y = -199.9550245361328
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
p9a0/corner_x = 526.3812578125
p9a0/corner_y = -355.68502453613274
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
p9a1/corner_x = 460.38225781249986
p9a1/corner_y = -355.88702453613274
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
p9a2/corner_x = 394.38225781250003
p9a2/corner_y = -356.0890245361328
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
p9a3/corner_x = 328.3822578125
p9a3/corner_y = -356.2910245361328
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
p9a4/corner_x = 262.3822578125
p9a4/corner_y = -356.4920245361328
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
p9a5/corner_x = 196.3832578125
p9a5/corner_y = -356.6940245361328
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
p9a6/corner_x = 130.3832578125
p9a6/corner_y = -356.89502453613284
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
p9a7/corner_x = 64.3830578125
p9a7/corner_y = -357.0970245361328
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
p10a0/corner_x = 526.5102578125
p10a0/corner_y = -512.2560245361326
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
p10a1/corner_x = 460.5102578125
p10a1/corner_y = -512.4300245361328
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
p10a2/corner_x = 394.5102578125
p10a2/corner_y = -512.6050245361326
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
p10a3/corner_x = 328.5102578125
p10a3/corner_y = -512.7780245361328
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
p10a4/corner_x = 262.5112578125
p10a4/corner_y = -512.9530245361327
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
p10a5/corner_x = 196.51125781249993
p10a5/corner_y = -513.1270245361328
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
p10a6/corner_x = 130.5112578125
p10a6/corner_y = -513.3000245361327
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
p10a7/corner_x = 64.5110578125
p10a7/corner_y = -513.4750245361328
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
p11a0/corner_x = 527.7092578125
p11a0/corner_y = -669.9690245361328
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
p11a1/corner_x = 461.70925781249997
p11a1/corner_y = -670.0980245361328
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
p11a2/corner_x = 395.7102578125
p11a2/corner_y = -670.2280245361328
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
p11a3/corner_x = 329.7102578125
p11a3/corner_y = -670.3580245361327
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
p11a4/corner_x = 263.71025781249995
p11a4/corner_y = -670.4860245361328
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
p11a5/corner_x = 197.71025781249998
p11a5/corner_y = -670.6160245361328
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
p11a6/corner_x = 131.71025781249998
p11a6/corner_y = -670.7460245361328
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
p11a7/corner_x = 65.7095578125
p11a7/corner_y = -670.8750245361327
p11a7/coffset = -0.000682

p12a0/dim1 = 12
p12a0/dim2 = ss
p12a0/dim3 = fs
p12a0/min_fs = 0
p12a0/min_ss = 0
p12a0/max_fs = 127
p12a0/max_ss = 63
p12a0/fs = -0.002172x +0.9999979999999997y
p12a0/ss = -0.9999979999999997x -0.002172y
p12a0/corner_x = 546.7303288726806
p12a0/corner_y = 427.28793801116944
p12a0/coffset = 0.002664

p12a1/dim1 = 12
p12a1/dim2 = ss
p12a1/dim3 = fs
p12a1/min_fs = 0
p12a1/min_ss = 64
p12a1/max_fs = 127
p12a1/max_ss = 127
p12a1/fs = -0.002172x +0.9999979999999997y
p12a1/ss = -0.9999979999999997x -0.002172y
p12a1/corner_x = 480.73132887268065
p12a1/corner_y = 427.1449380111694
p12a1/coffset = 0.002664

p12a2/dim1 = 12
p12a2/dim2 = ss
p12a2/dim3 = fs
p12a2/min_fs = 0
p12a2/min_ss = 128
p12a2/max_fs = 127
p12a2/max_ss = 191
p12a2/fs = -0.002172x +0.9999979999999997y
p12a2/ss = -0.9999979999999997x -0.002172y
p12a2/corner_x = 414.7303288726807
p12a2/corner_y = 427.0009380111693
p12a2/coffset = 0.002664

p12a3/dim1 = 12
p12a3/dim2 = ss
p12a3/dim3 = fs
p12a3/min_fs = 0
p12a3/min_ss = 192
p12a3/max_fs = 127
p12a3/max_ss = 255
p12a3/fs = -0.002172x +0.9999979999999997y
p12a3/ss = -0.9999979999999997x -0.002172y
p12a3/corner_x = 348.7303288726807
p12a3/corner_y = 426.85793801116944
p12a3/coffset = 0.002664

p12a4/dim1 = 12
p12a4/dim2 = ss
p12a4/dim3 = fs
p12a4/min_fs = 0
p12a4/min_ss = 256
p12a4/max_fs = 127
p12a4/max_ss = 319
p12a4/fs = -0.002172x +0.9999979999999997y
p12a4/ss = -0.9999979999999997x -0.002172y
p12a4/corner_x = 282.7323288726807
p12a4/corner_y = 426.7149380111694
p12a4/coffset = 0.002664

p12a5/dim1 = 12
p12a5/dim2 = ss
p12a5/dim3 = fs
p12a5/min_fs = 0
p12a5/min_ss = 320
p12a5/max_fs = 127
p12a5/max_ss = 383
p12a5/fs = -0.002172x +0.9999979999999997y
p12a5/ss = -0.9999979999999997x -0.002172y
p12a5/corner_x = 216.73232887268065
p12a5/corner_y = 426.5709380111693
p12a5/coffset = 0.002664

p12a6/dim1 = 12
p12a6/dim2 = ss
p12a6/dim3 = fs
p12a6/min_fs = 0
p12a6/min_ss = 384
p12a6/max_fs = 127
p12a6/max_ss = 447
p12a6/fs = -0.002172x +0.9999979999999997y
p12a6/ss = -0.9999979999999997x -0.002172y
p12a6/corner_x = 150.73132887268068
p12a6/corner_y = 426.4279380111694
p12a6/coffset = 0.002664

p12a7/dim1 = 12
p12a7/dim2 = ss
p12a7/dim3 = fs
p12a7/min_fs = 0
p12a7/min_ss = 448
p12a7/max_fs = 127
p12a7/max_ss = 511
p12a7/fs = -0.002172x +0.9999979999999997y
p12a7/ss = -0.9999979999999997x -0.002172y
p12a7/corner_x = 84.73142887268065
p12a7/corner_y = 426.28493801116946
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
p13a0/corner_x = 545.2703288726806
p13a0/corner_y = 279.7209380111693
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
p13a1/corner_x = 479.27032887268064
p13a1/corner_y = 279.2989380111694
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
p13a2/corner_x = 413.2723288726806
p13a2/corner_y = 278.8789380111694
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
p13a3/corner_x = 347.2733288726807
p13a3/corner_y = 278.4589380111693
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
p13a4/corner_x = 281.27432887268066
p13a4/corner_y = 278.0369380111694
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
p13a5/corner_x = 215.2763288726807
p13a5/corner_y = 277.61593801116936
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
p13a6/corner_x = 149.27732887268064
p13a6/corner_y = 277.19393801116945
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
p13a7/corner_x = 83.27862887268066
p13a7/corner_y = 276.77393801116943
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
p14a0/corner_x = 545.7413288726807
p14a0/corner_y = 122.71793801116938
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
p14a1/corner_x = 479.74132887268064
p14a1/corner_y = 122.36493801116943
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
p14a2/corner_x = 413.74332887268065
p14a2/corner_y = 122.01193801116943
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
p14a3/corner_x = 347.7433288726807
p14a3/corner_y = 121.65793801116942
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
p14a4/corner_x = 281.7443288726807
p14a4/corner_y = 121.30393801116945
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
p14a5/corner_x = 215.74632887268066
p14a5/corner_y = 120.94993801116944
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
p14a6/corner_x = 149.74732887268067
p14a6/corner_y = 120.59693801116944
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
p14a7/corner_x = 83.74742887268066
p14a7/corner_y = 120.24293801116943
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
p15a0/corner_x = 546.2503288726807
p15a0/corner_y = -33.38246198883056
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
p15a1/corner_x = 480.2523288726806
p15a1/corner_y = -33.651761988830565
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
p15a2/corner_x = 414.25232887268066
p15a2/corner_y = -33.921061988830566
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
p15a3/corner_x = 348.25332887268064
p15a3/corner_y = -34.19026198883057
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
p15a4/corner_x = 282.25332887268064
p15a4/corner_y = -34.459761988830564
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
p15a5/corner_x = 216.25432887268065
p15a5/corner_y = -34.72896198883056
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
p15a6/corner_x = 150.25332887268064
p15a6/corner_y = -34.99826198883056
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
p15a7/corner_x = 84.25412887268067
p15a7/corner_y = -35.26756198883056
p15a7/coffset = -0.000679
