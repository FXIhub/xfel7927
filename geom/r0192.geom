; AGIPD-1M geometry file written by EXtra-geom 1.12.0
; You may need to edit this file to add:
; - data and mask locations in the file
; - mask_good & mask_bad values to interpret the mask
; - adu_per_eV & photon_energy
; - clen (detector distance)
;
; See: http://www.desy.de/~twhite/crystfel/manual-crystfel_geometry.html

;XGEOM MOTORS=4,2
;XGEOM MOTOR_Q1=55.10418701171875,84.70984649658203
;XGEOM MOTOR_Q2=-9.022299766540527,-71.01100158691406
;XGEOM MOTOR_Q3=35.99871826171875,694.863525390625
;XGEOM MOTOR_Q4=2.006953239440918,-32.17916488647461

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
p0a0/corner_x = -533.1844708862304
p0a0/corner_y = 643.7174567871094
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
p0a1/corner_x = -467.1844708862305
p0a1/corner_y = 643.9394567871093
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
p0a2/corner_x = -401.1844708862304
p0a2/corner_y = 644.1614567871092
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
p0a3/corner_x = -335.18547088623046
p0a3/corner_y = 644.3834567871093
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
p0a4/corner_x = -269.18547088623046
p0a4/corner_y = 644.6064567871093
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
p0a5/corner_x = -203.1854708862304
p0a5/corner_y = 644.8274567871092
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
p0a6/corner_x = -137.18647088623047
p0a6/corner_y = 645.0494567871094
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
p0a7/corner_x = -71.18647088623045
p0a7/corner_y = 645.2714567871093
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
p1a0/corner_x = -532.6764708862304
p1a0/corner_y = 485.08845678710935
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
p1a1/corner_x = -466.6774708862305
p1a1/corner_y = 485.5144567871094
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
p1a2/corner_x = -400.6784708862304
p1a2/corner_y = 485.9424567871094
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
p1a3/corner_x = -334.68047088623047
p1a3/corner_y = 486.3684567871094
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
p1a4/corner_x = -268.68147088623044
p1a4/corner_y = 486.7964567871093
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
p1a5/corner_x = -202.68247088623048
p1a5/corner_y = 487.22245678710937
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
p1a6/corner_x = -136.68447088623046
p1a6/corner_y = 487.6504567871093
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
p1a7/corner_x = -70.68607088623047
p1a7/corner_y = 488.07645678710935
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
p2a0/corner_x = -532.0324708862305
p2a0/corner_y = 328.5884567871093
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
p2a1/corner_x = -466.03447088623045
p2a1/corner_y = 329.1084567871094
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
p2a2/corner_x = -400.03647088623046
p2a2/corner_y = 329.6274567871094
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
p2a3/corner_x = -334.0394708862304
p2a3/corner_y = 330.1464567871094
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
p2a4/corner_x = -268.0404708862305
p2a4/corner_y = 330.66745678710936
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
p2a5/corner_x = -202.04247088623046
p2a5/corner_y = 331.1854567871094
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
p2a6/corner_x = -136.04447088623044
p2a6/corner_y = 331.7064567871094
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
p2a7/corner_x = -70.04727088623046
p2a7/corner_y = 332.22545678710924
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
p3a0/corner_x = -531.7184708862304
p3a0/corner_y = 171.8354567871093
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
p3a1/corner_x = -465.71947088623045
p3a1/corner_y = 172.2744567871094
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
p3a2/corner_x = -399.72147088623046
p3a2/corner_y = 172.71445678710933
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
p3a3/corner_x = -333.72247088623044
p3a3/corner_y = 173.15345678710932
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
p3a4/corner_x = -267.7254708862304
p3a4/corner_y = 173.59445678710938
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
p3a5/corner_x = -201.72547088623048
p3a5/corner_y = 174.03345678710937
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
p3a6/corner_x = -135.72747088623046
p3a6/corner_y = 174.47345678710937
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
p3a7/corner_x = -69.72897088623047
p3a7/corner_y = 174.9134567871093
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
p8a0/corner_x = 525.5588681640625
p8a0/corner_y = -198.63079565429686
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
p8a1/corner_x = 459.55886816406246
p8a1/corner_y = -198.8207956542969
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
p8a2/corner_x = 393.5588681640625
p8a2/corner_y = -199.00979565429685
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
p8a3/corner_x = 327.5588681640625
p8a3/corner_y = -199.19879565429687
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
p8a4/corner_x = 261.5598681640625
p8a4/corner_y = -199.38879565429687
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
p8a5/corner_x = 195.55986816406244
p8a5/corner_y = -199.57679565429686
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
p8a6/corner_x = 129.55986816406246
p8a6/corner_y = -199.76679565429689
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
p8a7/corner_x = 63.560368164062474
p8a7/corner_y = -199.95479565429687
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
p9a0/corner_x = 526.3818681640626
p9a0/corner_y = -355.6847956542968
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
p9a1/corner_x = 460.38286816406236
p9a1/corner_y = -355.88679565429686
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
p9a2/corner_x = 394.3828681640625
p9a2/corner_y = -356.08879565429686
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
p9a3/corner_x = 328.3828681640625
p9a3/corner_y = -356.2907956542968
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
p9a4/corner_x = 262.3828681640625
p9a4/corner_y = -356.4917956542968
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
p9a5/corner_x = 196.3838681640625
p9a5/corner_y = -356.6937956542969
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
p9a6/corner_x = 130.3838681640625
p9a6/corner_y = -356.8947956542969
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
p9a7/corner_x = 64.3836681640625
p9a7/corner_y = -357.0967956542969
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
p10a0/corner_x = 526.5108681640625
p10a0/corner_y = -512.2557956542968
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
p10a1/corner_x = 460.5108681640625
p10a1/corner_y = -512.4297956542969
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
p10a2/corner_x = 394.5108681640625
p10a2/corner_y = -512.6047956542968
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
p10a3/corner_x = 328.51086816406246
p10a3/corner_y = -512.7777956542969
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
p10a4/corner_x = 262.5118681640625
p10a4/corner_y = -512.9527956542968
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
p10a5/corner_x = 196.51186816406243
p10a5/corner_y = -513.1267956542969
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
p10a6/corner_x = 130.5118681640625
p10a6/corner_y = -513.2997956542968
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
p10a7/corner_x = 64.51166816406248
p10a7/corner_y = -513.4747956542968
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
p11a0/corner_x = 527.7098681640625
p11a0/corner_y = -669.9687956542969
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
p11a1/corner_x = 461.7098681640624
p11a1/corner_y = -670.0977956542968
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
p11a2/corner_x = 395.7108681640625
p11a2/corner_y = -670.2277956542968
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
p11a3/corner_x = 329.7108681640625
p11a3/corner_y = -670.3577956542968
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
p11a4/corner_x = 263.7108681640625
p11a4/corner_y = -670.4857956542968
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
p11a5/corner_x = 197.7108681640625
p11a5/corner_y = -670.6157956542968
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
p11a6/corner_x = 131.71086816406248
p11a6/corner_y = -670.7457956542968
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
p11a7/corner_x = 65.7101681640625
p11a7/corner_y = -670.8747956542968
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
p12a0/corner_x = 546.733266189575
p12a0/corner_y = 427.2873562698364
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
p12a1/corner_x = 480.7342661895752
p12a1/corner_y = 427.1443562698364
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
p12a2/corner_x = 414.7332661895752
p12a2/corner_y = 427.0003562698363
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
p12a3/corner_x = 348.7332661895752
p12a3/corner_y = 426.8573562698364
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
p12a4/corner_x = 282.7352661895752
p12a4/corner_y = 426.71435626983634
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
p12a5/corner_x = 216.73526618957519
p12a5/corner_y = 426.5703562698363
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
p12a6/corner_x = 150.7342661895752
p12a6/corner_y = 426.42735626983637
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
p12a7/corner_x = 84.73436618957517
p12a7/corner_y = 426.28435626983645
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
p13a0/corner_x = 545.2732661895751
p13a0/corner_y = 279.7203562698363
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
p13a1/corner_x = 479.27326618957517
p13a1/corner_y = 279.2983562698364
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
p13a2/corner_x = 413.2752661895751
p13a2/corner_y = 278.8783562698364
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
p13a3/corner_x = 347.2762661895752
p13a3/corner_y = 278.4583562698363
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
p13a4/corner_x = 281.2772661895752
p13a4/corner_y = 278.0363562698364
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
p13a5/corner_x = 215.2792661895752
p13a5/corner_y = 277.61535626983635
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
p13a6/corner_x = 149.28026618957517
p13a6/corner_y = 277.19335626983644
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
p13a7/corner_x = 83.2815661895752
p13a7/corner_y = 276.7733562698364
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
p14a0/corner_x = 545.7442661895751
p14a0/corner_y = 122.71735626983637
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
p14a1/corner_x = 479.74426618957517
p14a1/corner_y = 122.36435626983642
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
p14a2/corner_x = 413.7462661895751
p14a2/corner_y = 122.01135626983644
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
p14a3/corner_x = 347.7462661895752
p14a3/corner_y = 121.65735626983641
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
p14a4/corner_x = 281.7472661895752
p14a4/corner_y = 121.30335626983644
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
p14a5/corner_x = 215.74926618957517
p14a5/corner_y = 120.94935626983643
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
p14a6/corner_x = 149.7502661895752
p14a6/corner_y = 120.59635626983643
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
p14a7/corner_x = 83.75036618957519
p14a7/corner_y = 120.24235626983642
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
p15a0/corner_x = 546.2532661895752
p15a0/corner_y = -33.38304373016358
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
p15a1/corner_x = 480.25526618957514
p15a1/corner_y = -33.65234373016358
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
p15a2/corner_x = 414.2552661895752
p15a2/corner_y = -33.92164373016358
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
p15a3/corner_x = 348.2562661895751
p15a3/corner_y = -34.190843730163586
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
p15a4/corner_x = 282.2562661895752
p15a4/corner_y = -34.46034373016358
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
p15a5/corner_x = 216.25726618957515
p15a5/corner_y = -34.72954373016358
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
p15a6/corner_x = 150.25626618957517
p15a6/corner_y = -34.99884373016358
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
p15a7/corner_x = 84.25706618957518
p15a7/corner_y = -35.26814373016357
p15a7/coffset = -0.000679
