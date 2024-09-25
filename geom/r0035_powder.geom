; AGIPD-1M geometry file written by EXtra-geom 1.11.0
; You may need to edit this file to add:
; - data and mask locations in the file
; - mask_good & mask_bad values to interpret the mask
; - adu_per_eV & photon_energy
; - clen (detector distance)
;
; See: http://www.desy.de/~twhite/crystfel/manual-crystfel_geometry.html

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
p0a0/fs = -0.001356x -0.9999939999999997y
p0a0/ss = +0.9999939999999997x -0.001356y
p0a0/corner_x = -518.7662501098633
p0a0/corner_y = 681.0337913269041
p0a0/coffset = 0.701089

p0a1/dim1 = 0
p0a1/dim2 = ss
p0a1/dim3 = fs
p0a1/min_fs = 0
p0a1/min_ss = 64
p0a1/max_fs = 127
p0a1/max_ss = 127
p0a1/fs = -0.001356x -0.9999939999999997y
p0a1/ss = +0.9999939999999997x -0.001356y
p0a1/corner_x = -452.76725010986326
p0a1/corner_y = 680.9237913269042
p0a1/coffset = 0.701089

p0a2/dim1 = 0
p0a2/dim2 = ss
p0a2/dim3 = fs
p0a2/min_fs = 0
p0a2/min_ss = 128
p0a2/max_fs = 127
p0a2/max_ss = 191
p0a2/fs = -0.001356x -0.9999939999999997y
p0a2/ss = +0.9999939999999997x -0.001356y
p0a2/corner_x = -386.77325010986317
p0a2/corner_y = 680.8127913269042
p0a2/coffset = 0.701089

p0a3/dim1 = 0
p0a3/dim2 = ss
p0a3/dim3 = fs
p0a3/min_fs = 0
p0a3/min_ss = 192
p0a3/max_fs = 127
p0a3/max_ss = 255
p0a3/fs = -0.001356x -0.9999939999999997y
p0a3/ss = +0.9999939999999997x -0.001356y
p0a3/corner_x = -320.7712501098632
p0a3/corner_y = 680.7027913269043
p0a3/coffset = 0.701089

p0a4/dim1 = 0
p0a4/dim2 = ss
p0a4/dim3 = fs
p0a4/min_fs = 0
p0a4/min_ss = 256
p0a4/max_fs = 127
p0a4/max_ss = 319
p0a4/fs = -0.001356x -0.9999939999999997y
p0a4/ss = +0.9999939999999997x -0.001356y
p0a4/corner_x = -254.7712501098632
p0a4/corner_y = 680.5877913269043
p0a4/coffset = 0.701089

p0a5/dim1 = 0
p0a5/dim2 = ss
p0a5/dim3 = fs
p0a5/min_fs = 0
p0a5/min_ss = 320
p0a5/max_fs = 127
p0a5/max_ss = 383
p0a5/fs = -0.001356x -0.9999939999999997y
p0a5/ss = +0.9999939999999997x -0.001356y
p0a5/corner_x = -188.7692501098633
p0a5/corner_y = 680.5887913269041
p0a5/coffset = 0.701089

p0a6/dim1 = 0
p0a6/dim2 = ss
p0a6/dim3 = fs
p0a6/min_fs = 0
p0a6/min_ss = 384
p0a6/max_fs = 127
p0a6/max_ss = 447
p0a6/fs = -0.001356x -0.9999939999999997y
p0a6/ss = +0.9999939999999997x -0.001356y
p0a6/corner_x = -122.76925010986325
p0a6/corner_y = 680.4987913269042
p0a6/coffset = 0.701089

p0a7/dim1 = 0
p0a7/dim2 = ss
p0a7/dim3 = fs
p0a7/min_fs = 0
p0a7/min_ss = 448
p0a7/max_fs = 127
p0a7/max_ss = 511
p0a7/fs = -0.001356x -0.9999939999999997y
p0a7/ss = +0.9999939999999997x -0.001356y
p0a7/corner_x = -56.76875010986328
p0a7/corner_y = 680.4087913269041
p0a7/coffset = 0.701089

p1a0/dim1 = 1
p1a0/dim2 = ss
p1a0/dim3 = fs
p1a0/min_fs = 0
p1a0/min_ss = 0
p1a0/max_fs = 127
p1a0/max_ss = 63
p1a0/fs = -0.000066x -1.y
p1a0/ss = +1.x -0.000066y
p1a0/corner_x = -519.5572501098632
p1a0/corner_y = 523.7587913269041
p1a0/coffset = 0.700911

p1a1/dim1 = 1
p1a1/dim2 = ss
p1a1/dim3 = fs
p1a1/min_fs = 0
p1a1/min_ss = 64
p1a1/max_fs = 127
p1a1/max_ss = 127
p1a1/fs = -0.000066x -1.y
p1a1/ss = +1.x -0.000066y
p1a1/corner_x = -453.5582501098633
p1a1/corner_y = 523.7657913269042
p1a1/coffset = 0.700911

p1a2/dim1 = 1
p1a2/dim2 = ss
p1a2/dim3 = fs
p1a2/min_fs = 0
p1a2/min_ss = 128
p1a2/max_fs = 127
p1a2/max_ss = 191
p1a2/fs = -0.000066x -1.y
p1a2/ss = +1.x -0.000066y
p1a2/corner_x = -387.5592501098632
p1a2/corner_y = 523.7737913269042
p1a2/coffset = 0.700911

p1a3/dim1 = 1
p1a3/dim2 = ss
p1a3/dim3 = fs
p1a3/min_fs = 0
p1a3/min_ss = 192
p1a3/max_fs = 127
p1a3/max_ss = 255
p1a3/fs = -0.000066x -1.y
p1a3/ss = +1.x -0.000066y
p1a3/corner_x = -321.5632501098633
p1a3/corner_y = 523.7797913269042
p1a3/coffset = 0.700911

p1a4/dim1 = 1
p1a4/dim2 = ss
p1a4/dim3 = fs
p1a4/min_fs = 0
p1a4/min_ss = 256
p1a4/max_fs = 127
p1a4/max_ss = 319
p1a4/fs = -0.000066x -1.y
p1a4/ss = +1.x -0.000066y
p1a4/corner_x = -255.56525010986317
p1a4/corner_y = 523.7877913269041
p1a4/coffset = 0.700911

p1a5/dim1 = 1
p1a5/dim2 = ss
p1a5/dim3 = fs
p1a5/min_fs = 0
p1a5/min_ss = 320
p1a5/max_fs = 127
p1a5/max_ss = 383
p1a5/fs = -0.000066x -1.y
p1a5/ss = +1.x -0.000066y
p1a5/corner_x = -189.56925010986325
p1a5/corner_y = 523.7947913269041
p1a5/coffset = 0.700911

p1a6/dim1 = 1
p1a6/dim2 = ss
p1a6/dim3 = fs
p1a6/min_fs = 0
p1a6/min_ss = 384
p1a6/max_fs = 127
p1a6/max_ss = 447
p1a6/fs = -0.000066x -1.y
p1a6/ss = +1.x -0.000066y
p1a6/corner_x = -123.54925010986327
p1a6/corner_y = 523.8037913269042
p1a6/coffset = 0.700911

p1a7/dim1 = 1
p1a7/dim2 = ss
p1a7/dim3 = fs
p1a7/min_fs = 0
p1a7/min_ss = 448
p1a7/max_fs = 127
p1a7/max_ss = 511
p1a7/fs = -0.000066x -1.y
p1a7/ss = +1.x -0.000066y
p1a7/corner_x = -57.55045010986326
p1a7/corner_y = 523.8107913269041
p1a7/coffset = 0.700911

p2a0/dim1 = 2
p2a0/dim2 = ss
p2a0/dim3 = fs
p2a0/min_fs = 0
p2a0/min_ss = 0
p2a0/max_fs = 127
p2a0/max_ss = 63
p2a0/fs = +0.0013719999999999997x -1.y
p2a0/ss = +1.x +0.0013719999999999997y
p2a0/corner_x = -519.0262501098632
p2a0/corner_y = 366.1547913269043
p2a0/coffset = 0.700662

p2a1/dim1 = 2
p2a1/dim2 = ss
p2a1/dim3 = fs
p2a1/min_fs = 0
p2a1/min_ss = 64
p2a1/max_fs = 127
p2a1/max_ss = 127
p2a1/fs = +0.0013719999999999997x -1.y
p2a1/ss = +1.x +0.0013719999999999997y
p2a1/corner_x = -453.02825010986317
p2a1/corner_y = 366.23679132690404
p2a1/coffset = 0.700662

p2a2/dim1 = 2
p2a2/dim2 = ss
p2a2/dim3 = fs
p2a2/min_fs = 0
p2a2/min_ss = 128
p2a2/max_fs = 127
p2a2/max_ss = 191
p2a2/fs = +0.0013719999999999997x -1.y
p2a2/ss = +1.x +0.0013719999999999997y
p2a2/corner_x = -387.0292501098632
p2a2/corner_y = 366.31979132690407
p2a2/coffset = 0.700662

p2a3/dim1 = 2
p2a3/dim2 = ss
p2a3/dim3 = fs
p2a3/min_fs = 0
p2a3/min_ss = 192
p2a3/max_fs = 127
p2a3/max_ss = 255
p2a3/fs = +0.0013719999999999997x -1.y
p2a3/ss = +1.x +0.0013719999999999997y
p2a3/corner_x = -321.02925010986326
p2a3/corner_y = 366.4037913269041
p2a3/coffset = 0.700662

p2a4/dim1 = 2
p2a4/dim2 = ss
p2a4/dim3 = fs
p2a4/min_fs = 0
p2a4/min_ss = 256
p2a4/max_fs = 127
p2a4/max_ss = 319
p2a4/fs = +0.0013719999999999997x -1.y
p2a4/ss = +1.x +0.0013719999999999997y
p2a4/corner_x = -255.03025010986323
p2a4/corner_y = 366.48479132690426
p2a4/coffset = 0.700662

p2a5/dim1 = 2
p2a5/dim2 = ss
p2a5/dim3 = fs
p2a5/min_fs = 0
p2a5/min_ss = 320
p2a5/max_fs = 127
p2a5/max_ss = 383
p2a5/fs = +0.0013719999999999997x -1.y
p2a5/ss = +1.x +0.0013719999999999997y
p2a5/corner_x = -189.03225010986327
p2a5/corner_y = 366.5677913269043
p2a5/coffset = 0.700662

p2a6/dim1 = 2
p2a6/dim2 = ss
p2a6/dim3 = fs
p2a6/min_fs = 0
p2a6/min_ss = 384
p2a6/max_fs = 127
p2a6/max_ss = 447
p2a6/fs = +0.0013719999999999997x -1.y
p2a6/ss = +1.x +0.0013719999999999997y
p2a6/corner_x = -123.01925010986324
p2a6/corner_y = 366.65079132690414
p2a6/coffset = 0.700662

p2a7/dim1 = 2
p2a7/dim2 = ss
p2a7/dim3 = fs
p2a7/min_fs = 0
p2a7/min_ss = 448
p2a7/max_fs = 127
p2a7/max_ss = 511
p2a7/fs = +0.0013719999999999997x -1.y
p2a7/ss = +1.x +0.0013719999999999997y
p2a7/corner_x = -57.015850109863266
p2a7/corner_y = 366.7317913269042
p2a7/coffset = 0.700662

p3a0/dim1 = 3
p3a0/dim2 = ss
p3a0/dim3 = fs
p3a0/min_fs = 0
p3a0/min_ss = 0
p3a0/max_fs = 127
p3a0/max_ss = 63
p3a0/fs = +0.003912x -0.9999939999999997y
p3a0/ss = +0.9999939999999997x +0.003912y
p3a0/corner_x = -519.4152501098632
p3a0/corner_y = 208.6847913269043
p3a0/coffset = 0.700897

p3a1/dim1 = 3
p3a1/dim2 = ss
p3a1/dim3 = fs
p3a1/min_fs = 0
p3a1/min_ss = 64
p3a1/max_fs = 127
p3a1/max_ss = 127
p3a1/fs = +0.003912x -0.9999939999999997y
p3a1/ss = +0.9999939999999997x +0.003912y
p3a1/corner_x = -453.4202501098632
p3a1/corner_y = 208.94679132690428
p3a1/coffset = 0.700897

p3a2/dim1 = 3
p3a2/dim2 = ss
p3a2/dim3 = fs
p3a2/min_fs = 0
p3a2/min_ss = 128
p3a2/max_fs = 127
p3a2/max_ss = 191
p3a2/fs = +0.003912x -0.9999939999999997y
p3a2/ss = +0.9999939999999997x +0.003912y
p3a2/corner_x = -387.4222501098633
p3a2/corner_y = 209.20579132690426
p3a2/coffset = 0.700897

p3a3/dim1 = 3
p3a3/dim2 = ss
p3a3/dim3 = fs
p3a3/min_fs = 0
p3a3/min_ss = 192
p3a3/max_fs = 127
p3a3/max_ss = 255
p3a3/fs = +0.003912x -0.9999939999999997y
p3a3/ss = +0.9999939999999997x +0.003912y
p3a3/corner_x = -321.42825010986314
p3a3/corner_y = 209.4647913269043
p3a3/coffset = 0.700897

p3a4/dim1 = 3
p3a4/dim2 = ss
p3a4/dim3 = fs
p3a4/min_fs = 0
p3a4/min_ss = 256
p3a4/max_fs = 127
p3a4/max_ss = 319
p3a4/fs = +0.003912x -0.9999939999999997y
p3a4/ss = +0.9999939999999997x +0.003912y
p3a4/corner_x = -255.42725010986317
p3a4/corner_y = 209.7247913269043
p3a4/coffset = 0.700897

p3a5/dim1 = 3
p3a5/dim2 = ss
p3a5/dim3 = fs
p3a5/min_fs = 0
p3a5/min_ss = 320
p3a5/max_fs = 127
p3a5/max_ss = 383
p3a5/fs = +0.003912x -0.9999939999999997y
p3a5/ss = +0.9999939999999997x +0.003912y
p3a5/corner_x = -189.41325010986327
p3a5/corner_y = 209.98579132690426
p3a5/coffset = 0.700897

p3a6/dim1 = 3
p3a6/dim2 = ss
p3a6/dim3 = fs
p3a6/min_fs = 0
p3a6/min_ss = 384
p3a6/max_fs = 127
p3a6/max_ss = 447
p3a6/fs = +0.003912x -0.9999939999999997y
p3a6/ss = +0.9999939999999997x +0.003912y
p3a6/corner_x = -123.41025010986324
p3a6/corner_y = 210.24779132690423
p3a6/coffset = 0.700897

p3a7/dim1 = 3
p3a7/dim2 = ss
p3a7/dim3 = fs
p3a7/min_fs = 0
p3a7/min_ss = 448
p3a7/max_fs = 127
p3a7/max_ss = 511
p3a7/fs = +0.003912x -0.9999939999999997y
p3a7/ss = +0.9999939999999997x +0.003912y
p3a7/corner_x = -57.412850109863264
p3a7/corner_y = 210.50779132690425
p3a7/coffset = 0.700897

p4a0/dim1 = 4
p4a0/dim2 = ss
p4a0/dim3 = fs
p4a0/min_fs = 0
p4a0/min_ss = 0
p4a0/max_fs = 127
p4a0/max_ss = 63
p4a0/fs = -0.001188x -1.000001y
p4a0/ss = +1.000001x -0.001188y
p4a0/corner_x = -550.4834198974609
p4a0/corner_y = 42.10144364318847
p4a0/coffset = 0.700727

p4a1/dim1 = 4
p4a1/dim2 = ss
p4a1/dim3 = fs
p4a1/min_fs = 0
p4a1/min_ss = 64
p4a1/max_fs = 127
p4a1/max_ss = 127
p4a1/fs = -0.001188x -1.000001y
p4a1/ss = +1.000001x -0.001188y
p4a1/corner_x = -484.48441989746084
p4a1/corner_y = 42.03124364318847
p4a1/coffset = 0.700727

p4a2/dim1 = 4
p4a2/dim2 = ss
p4a2/dim3 = fs
p4a2/min_fs = 0
p4a2/min_ss = 128
p4a2/max_fs = 127
p4a2/max_ss = 191
p4a2/fs = -0.001188x -1.000001y
p4a2/ss = +1.000001x -0.001188y
p4a2/corner_x = -418.4894198974609
p4a2/corner_y = 41.96104364318847
p4a2/coffset = 0.700727

p4a3/dim1 = 4
p4a3/dim2 = ss
p4a3/dim3 = fs
p4a3/min_fs = 0
p4a3/min_ss = 192
p4a3/max_fs = 127
p4a3/max_ss = 255
p4a3/fs = -0.001188x -1.000001y
p4a3/ss = +1.000001x -0.001188y
p4a3/corner_x = -352.4904198974609
p4a3/corner_y = 41.89084364318847
p4a3/coffset = 0.700727

p4a4/dim1 = 4
p4a4/dim2 = ss
p4a4/dim3 = fs
p4a4/min_fs = 0
p4a4/min_ss = 256
p4a4/max_fs = 127
p4a4/max_ss = 319
p4a4/fs = -0.001188x -1.000001y
p4a4/ss = +1.000001x -0.001188y
p4a4/corner_x = -286.4914198974609
p4a4/corner_y = 41.82104364318847
p4a4/coffset = 0.700727

p4a5/dim1 = 4
p4a5/dim2 = ss
p4a5/dim3 = fs
p4a5/min_fs = 0
p4a5/min_ss = 320
p4a5/max_fs = 127
p4a5/max_ss = 383
p4a5/fs = -0.001188x -1.000001y
p4a5/ss = +1.000001x -0.001188y
p4a5/corner_x = -220.49341989746094
p4a5/corner_y = 41.75094364318847
p4a5/coffset = 0.700727

p4a6/dim1 = 4
p4a6/dim2 = ss
p4a6/dim3 = fs
p4a6/min_fs = 0
p4a6/min_ss = 384
p4a6/max_fs = 127
p4a6/max_ss = 447
p4a6/fs = -0.001188x -1.000001y
p4a6/ss = +1.000001x -0.001188y
p4a6/corner_x = -154.4964198974609
p4a6/corner_y = 41.68094364318847
p4a6/coffset = 0.700727

p4a7/dim1 = 4
p4a7/dim2 = ss
p4a7/dim3 = fs
p4a7/min_fs = 0
p4a7/min_ss = 448
p4a7/max_fs = 127
p4a7/max_ss = 511
p4a7/fs = -0.001188x -1.000001y
p4a7/ss = +1.000001x -0.001188y
p4a7/corner_x = -88.50031989746091
p4a7/corner_y = 41.61074364318847
p4a7/coffset = 0.700727

p5a0/dim1 = 5
p5a0/dim2 = ss
p5a0/dim3 = fs
p5a0/min_fs = 0
p5a0/min_ss = 0
p5a0/max_fs = 127
p5a0/max_ss = 63
p5a0/fs = +0.0017189999999999996x -1.y
p5a0/ss = +1.x +0.0017189999999999996y
p5a0/corner_x = -548.8724198974609
p5a0/corner_y = -115.7463563568115
p5a0/coffset = 0.7008399999999999

p5a1/dim1 = 5
p5a1/dim2 = ss
p5a1/dim3 = fs
p5a1/min_fs = 0
p5a1/min_ss = 64
p5a1/max_fs = 127
p5a1/max_ss = 127
p5a1/fs = +0.0017189999999999996x -1.y
p5a1/ss = +1.x +0.0017189999999999996y
p5a1/corner_x = -482.8744198974608
p5a1/corner_y = -115.63735635681152
p5a1/coffset = 0.7008399999999999

p5a2/dim1 = 5
p5a2/dim2 = ss
p5a2/dim3 = fs
p5a2/min_fs = 0
p5a2/min_ss = 128
p5a2/max_fs = 127
p5a2/max_ss = 191
p5a2/fs = +0.0017189999999999996x -1.y
p5a2/ss = +1.x +0.0017189999999999996y
p5a2/corner_x = -416.87741989746075
p5a2/corner_y = -115.53335635681154
p5a2/coffset = 0.7008399999999999

p5a3/dim1 = 5
p5a3/dim2 = ss
p5a3/dim3 = fs
p5a3/min_fs = 0
p5a3/min_ss = 192
p5a3/max_fs = 127
p5a3/max_ss = 255
p5a3/fs = +0.0017189999999999996x -1.y
p5a3/ss = +1.x +0.0017189999999999996y
p5a3/corner_x = -350.8804198974609
p5a3/corner_y = -115.42835635681152
p5a3/coffset = 0.7008399999999999

p5a4/dim1 = 5
p5a4/dim2 = ss
p5a4/dim3 = fs
p5a4/min_fs = 0
p5a4/min_ss = 256
p5a4/max_fs = 127
p5a4/max_ss = 319
p5a4/fs = +0.0017189999999999996x -1.y
p5a4/ss = +1.x +0.0017189999999999996y
p5a4/corner_x = -284.8814198974609
p5a4/corner_y = -115.31935635681153
p5a4/coffset = 0.7008399999999999

p5a5/dim1 = 5
p5a5/dim2 = ss
p5a5/dim3 = fs
p5a5/min_fs = 0
p5a5/min_ss = 320
p5a5/max_fs = 127
p5a5/max_ss = 383
p5a5/fs = +0.0017189999999999996x -1.y
p5a5/ss = +1.x +0.0017189999999999996y
p5a5/corner_x = -218.8684198974609
p5a5/corner_y = -115.2143563568115
p5a5/coffset = 0.7008399999999999

p5a6/dim1 = 5
p5a6/dim2 = ss
p5a6/dim3 = fs
p5a6/min_fs = 0
p5a6/min_ss = 384
p5a6/max_fs = 127
p5a6/max_ss = 447
p5a6/fs = +0.0017189999999999996x -1.y
p5a6/ss = +1.x +0.0017189999999999996y
p5a6/corner_x = -152.8674198974609
p5a6/corner_y = -115.10835635681151
p5a6/coffset = 0.7008399999999999

p5a7/dim1 = 5
p5a7/dim2 = ss
p5a7/dim3 = fs
p5a7/min_fs = 0
p5a7/min_ss = 448
p5a7/max_fs = 127
p5a7/max_ss = 511
p5a7/fs = +0.0017189999999999996x -1.y
p5a7/ss = +1.x +0.0017189999999999996y
p5a7/corner_x = -86.86741989746089
p5a7/corner_y = -115.0003563568115
p5a7/coffset = 0.7008399999999999

p6a0/dim1 = 6
p6a0/dim2 = ss
p6a0/dim3 = fs
p6a0/min_fs = 0
p6a0/min_ss = 0
p6a0/max_fs = 127
p6a0/max_ss = 63
p6a0/fs = -0.001441x -1.y
p6a0/ss = +1.x -0.001441y
p6a0/corner_x = -549.6504198974609
p6a0/corner_y = -272.2323563568115
p6a0/coffset = 0.7009029999999999

p6a1/dim1 = 6
p6a1/dim2 = ss
p6a1/dim3 = fs
p6a1/min_fs = 0
p6a1/min_ss = 64
p6a1/max_fs = 127
p6a1/max_ss = 127
p6a1/fs = -0.001441x -1.y
p6a1/ss = +1.x -0.001441y
p6a1/corner_x = -483.6534198974609
p6a1/corner_y = -272.3373563568115
p6a1/coffset = 0.7009029999999999

p6a2/dim1 = 6
p6a2/dim2 = ss
p6a2/dim3 = fs
p6a2/min_fs = 0
p6a2/min_ss = 128
p6a2/max_fs = 127
p6a2/max_ss = 191
p6a2/fs = -0.001441x -1.y
p6a2/ss = +1.x -0.001441y
p6a2/corner_x = -417.65541989746094
p6a2/corner_y = -272.44235635681156
p6a2/coffset = 0.7009029999999999

p6a3/dim1 = 6
p6a3/dim2 = ss
p6a3/dim3 = fs
p6a3/min_fs = 0
p6a3/min_ss = 192
p6a3/max_fs = 127
p6a3/max_ss = 255
p6a3/fs = -0.001441x -1.y
p6a3/ss = +1.x -0.001441y
p6a3/corner_x = -351.65641989746086
p6a3/corner_y = -272.5453563568115
p6a3/coffset = 0.7009029999999999

p6a4/dim1 = 6
p6a4/dim2 = ss
p6a4/dim3 = fs
p6a4/min_fs = 0
p6a4/min_ss = 256
p6a4/max_fs = 127
p6a4/max_ss = 319
p6a4/fs = -0.001441x -1.y
p6a4/ss = +1.x -0.001441y
p6a4/corner_x = -285.65941989746085
p6a4/corner_y = -272.6483563568115
p6a4/coffset = 0.7009029999999999

p6a5/dim1 = 6
p6a5/dim2 = ss
p6a5/dim3 = fs
p6a5/min_fs = 0
p6a5/min_ss = 320
p6a5/max_fs = 127
p6a5/max_ss = 383
p6a5/fs = -0.001441x -1.y
p6a5/ss = +1.x -0.001441y
p6a5/corner_x = -219.65941989746094
p6a5/corner_y = -272.75135635681147
p6a5/coffset = 0.7009029999999999

p6a6/dim1 = 6
p6a6/dim2 = ss
p6a6/dim3 = fs
p6a6/min_fs = 0
p6a6/min_ss = 384
p6a6/max_fs = 127
p6a6/max_ss = 447
p6a6/fs = -0.001441x -1.y
p6a6/ss = +1.x -0.001441y
p6a6/corner_x = -153.6614198974609
p6a6/corner_y = -272.85435635681154
p6a6/coffset = 0.7009029999999999

p6a7/dim1 = 6
p6a7/dim2 = ss
p6a7/dim3 = fs
p6a7/min_fs = 0
p6a7/min_ss = 448
p6a7/max_fs = 127
p6a7/max_ss = 511
p6a7/fs = -0.001441x -1.y
p6a7/ss = +1.x -0.001441y
p6a7/corner_x = -87.6642198974609
p6a7/corner_y = -272.96035635681153
p6a7/coffset = 0.7009029999999999

p7a0/dim1 = 7
p7a0/dim2 = ss
p7a0/dim3 = fs
p7a0/min_fs = 0
p7a0/min_ss = 0
p7a0/max_fs = 127
p7a0/max_ss = 63
p7a0/fs = +0.001693x -0.999996y
p7a0/ss = +0.999996x +0.001693y
p7a0/corner_x = -550.1154198974608
p7a0/corner_y = -432.27035635681153
p7a0/coffset = 0.701398

p7a1/dim1 = 7
p7a1/dim2 = ss
p7a1/dim3 = fs
p7a1/min_fs = 0
p7a1/min_ss = 64
p7a1/max_fs = 127
p7a1/max_ss = 127
p7a1/fs = +0.001693x -0.999996y
p7a1/ss = +0.999996x +0.001693y
p7a1/corner_x = -484.12141989746095
p7a1/corner_y = -432.1593563568115
p7a1/coffset = 0.701398

p7a2/dim1 = 7
p7a2/dim2 = ss
p7a2/dim3 = fs
p7a2/min_fs = 0
p7a2/min_ss = 128
p7a2/max_fs = 127
p7a2/max_ss = 191
p7a2/fs = +0.001693x -0.999996y
p7a2/ss = +0.999996x +0.001693y
p7a2/corner_x = -418.1254198974609
p7a2/corner_y = -432.0533563568115
p7a2/coffset = 0.701398

p7a3/dim1 = 7
p7a3/dim2 = ss
p7a3/dim3 = fs
p7a3/min_fs = 0
p7a3/min_ss = 192
p7a3/max_fs = 127
p7a3/max_ss = 255
p7a3/fs = +0.001693x -0.999996y
p7a3/ss = +0.999996x +0.001693y
p7a3/corner_x = -352.1264198974608
p7a3/corner_y = -431.94235635681144
p7a3/coffset = 0.701398

p7a4/dim1 = 7
p7a4/dim2 = ss
p7a4/dim3 = fs
p7a4/min_fs = 0
p7a4/min_ss = 256
p7a4/max_fs = 127
p7a4/max_ss = 319
p7a4/fs = +0.001693x -0.999996y
p7a4/ss = +0.999996x +0.001693y
p7a4/corner_x = -286.1324198974609
p7a4/corner_y = -431.8333563568115
p7a4/coffset = 0.701398

p7a5/dim1 = 7
p7a5/dim2 = ss
p7a5/dim3 = fs
p7a5/min_fs = 0
p7a5/min_ss = 320
p7a5/max_fs = 127
p7a5/max_ss = 383
p7a5/fs = +0.001693x -0.999996y
p7a5/ss = +0.999996x +0.001693y
p7a5/corner_x = -220.13541989746088
p7a5/corner_y = -431.7253563568115
p7a5/coffset = 0.701398

p7a6/dim1 = 7
p7a6/dim2 = ss
p7a6/dim3 = fs
p7a6/min_fs = 0
p7a6/min_ss = 384
p7a6/max_fs = 127
p7a6/max_ss = 447
p7a6/fs = +0.001693x -0.999996y
p7a6/ss = +0.999996x +0.001693y
p7a6/corner_x = -154.1374198974609
p7a6/corner_y = -431.6153563568116
p7a6/coffset = 0.701398

p7a7/dim1 = 7
p7a7/dim2 = ss
p7a7/dim3 = fs
p7a7/min_fs = 0
p7a7/min_ss = 448
p7a7/max_fs = 127
p7a7/max_ss = 511
p7a7/fs = +0.001693x -0.999996y
p7a7/ss = +0.999996x +0.001693y
p7a7/corner_x = -88.13941989746091
p7a7/corner_y = -431.50535635681143
p7a7/coffset = 0.701398

p8a0/dim1 = 8
p8a0/dim2 = ss
p8a0/dim3 = fs
p8a0/min_fs = 0
p8a0/min_ss = 0
p8a0/max_fs = 127
p8a0/max_ss = 63
p8a0/fs = -0.00044199999999999996x +1.y
p8a0/ss = -1.x -0.00044199999999999996y
p8a0/corner_x = 526.0707128906248
p8a0/corner_y = -203.45838723144524
p8a0/coffset = 0.7005889999999999

p8a1/dim1 = 8
p8a1/dim2 = ss
p8a1/dim3 = fs
p8a1/min_fs = 0
p8a1/min_ss = 64
p8a1/max_fs = 127
p8a1/max_ss = 127
p8a1/fs = -0.00044199999999999996x +1.y
p8a1/ss = -1.x -0.00044199999999999996y
p8a1/corner_x = 460.0727128906247
p8a1/corner_y = -203.49138723144526
p8a1/coffset = 0.7005889999999999

p8a2/dim1 = 8
p8a2/dim2 = ss
p8a2/dim3 = fs
p8a2/min_fs = 0
p8a2/min_ss = 128
p8a2/max_fs = 127
p8a2/max_ss = 191
p8a2/fs = -0.00044199999999999996x +1.y
p8a2/ss = -1.x -0.00044199999999999996y
p8a2/corner_x = 394.0737128906248
p8a2/corner_y = -203.52538723144528
p8a2/coffset = 0.7005889999999999

p8a3/dim1 = 8
p8a3/dim2 = ss
p8a3/dim3 = fs
p8a3/min_fs = 0
p8a3/min_ss = 192
p8a3/max_fs = 127
p8a3/max_ss = 255
p8a3/fs = -0.00044199999999999996x +1.y
p8a3/ss = -1.x -0.00044199999999999996y
p8a3/corner_x = 328.0767128906248
p8a3/corner_y = -203.55838723144524
p8a3/coffset = 0.7005889999999999

p8a4/dim1 = 8
p8a4/dim2 = ss
p8a4/dim3 = fs
p8a4/min_fs = 0
p8a4/min_ss = 256
p8a4/max_fs = 127
p8a4/max_ss = 319
p8a4/fs = -0.00044199999999999996x +1.y
p8a4/ss = -1.x -0.00044199999999999996y
p8a4/corner_x = 262.07671289062483
p8a4/corner_y = -203.59638723144522
p8a4/coffset = 0.7005889999999999

p8a5/dim1 = 8
p8a5/dim2 = ss
p8a5/dim3 = fs
p8a5/min_fs = 0
p8a5/min_ss = 320
p8a5/max_fs = 127
p8a5/max_ss = 383
p8a5/fs = -0.00044199999999999996x +1.y
p8a5/ss = -1.x -0.00044199999999999996y
p8a5/corner_x = 196.07871289062473
p8a5/corner_y = -203.6283872314452
p8a5/coffset = 0.7005889999999999

p8a6/dim1 = 8
p8a6/dim2 = ss
p8a6/dim3 = fs
p8a6/min_fs = 0
p8a6/min_ss = 384
p8a6/max_fs = 127
p8a6/max_ss = 447
p8a6/fs = -0.00044199999999999996x +1.y
p8a6/ss = -1.x -0.00044199999999999996y
p8a6/corner_x = 130.08171289062483
p8a6/corner_y = -203.66338723144528
p8a6/coffset = 0.7005889999999999

p8a7/dim1 = 8
p8a7/dim2 = ss
p8a7/dim3 = fs
p8a7/min_fs = 0
p8a7/min_ss = 448
p8a7/max_fs = 127
p8a7/max_ss = 511
p8a7/fs = -0.00044199999999999996x +1.y
p8a7/ss = -1.x -0.00044199999999999996y
p8a7/corner_x = 64.08411289062482
p8a7/corner_y = -203.70138723144527
p8a7/coffset = 0.7005889999999999

p9a0/dim1 = 9
p9a0/dim2 = ss
p9a0/dim3 = fs
p9a0/min_fs = 0
p9a0/min_ss = 0
p9a0/max_fs = 127
p9a0/max_ss = 63
p9a0/fs = -0.0034019999999999996x +0.999995y
p9a0/ss = -0.999995x -0.0034019999999999996y
p9a0/corner_x = 525.4177128906248
p9a0/corner_y = -360.0603872314453
p9a0/coffset = 0.700742

p9a1/dim1 = 9
p9a1/dim2 = ss
p9a1/dim3 = fs
p9a1/min_fs = 0
p9a1/min_ss = 64
p9a1/max_fs = 127
p9a1/max_ss = 127
p9a1/fs = -0.0034019999999999996x +0.999995y
p9a1/ss = -0.999995x -0.0034019999999999996y
p9a1/corner_x = 459.4207128906247
p9a1/corner_y = -360.2893872314452
p9a1/coffset = 0.700742

p9a2/dim1 = 9
p9a2/dim2 = ss
p9a2/dim3 = fs
p9a2/min_fs = 0
p9a2/min_ss = 128
p9a2/max_fs = 127
p9a2/max_ss = 191
p9a2/fs = -0.0034019999999999996x +0.999995y
p9a2/ss = -0.999995x -0.0034019999999999996y
p9a2/corner_x = 393.42271289062484
p9a2/corner_y = -360.51638723144526
p9a2/coffset = 0.700742

p9a3/dim1 = 9
p9a3/dim2 = ss
p9a3/dim3 = fs
p9a3/min_fs = 0
p9a3/min_ss = 192
p9a3/max_fs = 127
p9a3/max_ss = 255
p9a3/fs = -0.0034019999999999996x +0.999995y
p9a3/ss = -0.999995x -0.0034019999999999996y
p9a3/corner_x = 327.4277128906248
p9a3/corner_y = -360.74338723144524
p9a3/coffset = 0.700742

p9a4/dim1 = 9
p9a4/dim2 = ss
p9a4/dim3 = fs
p9a4/min_fs = 0
p9a4/min_ss = 256
p9a4/max_fs = 127
p9a4/max_ss = 319
p9a4/fs = -0.0034019999999999996x +0.999995y
p9a4/ss = -0.999995x -0.0034019999999999996y
p9a4/corner_x = 261.4297128906247
p9a4/corner_y = -360.97238723144534
p9a4/coffset = 0.700742

p9a5/dim1 = 9
p9a5/dim2 = ss
p9a5/dim3 = fs
p9a5/min_fs = 0
p9a5/min_ss = 320
p9a5/max_fs = 127
p9a5/max_ss = 383
p9a5/fs = -0.0034019999999999996x +0.999995y
p9a5/ss = -0.999995x -0.0034019999999999996y
p9a5/corner_x = 195.43371289062478
p9a5/corner_y = -361.2003872314453
p9a5/coffset = 0.700742

p9a6/dim1 = 9
p9a6/dim2 = ss
p9a6/dim3 = fs
p9a6/min_fs = 0
p9a6/min_ss = 384
p9a6/max_fs = 127
p9a6/max_ss = 447
p9a6/fs = -0.0034019999999999996x +0.999995y
p9a6/ss = -0.999995x -0.0034019999999999996y
p9a6/corner_x = 129.41371289062482
p9a6/corner_y = -361.4273872314454
p9a6/coffset = 0.700742

p9a7/dim1 = 9
p9a7/dim2 = ss
p9a7/dim3 = fs
p9a7/min_fs = 0
p9a7/min_ss = 448
p9a7/max_fs = 127
p9a7/max_ss = 511
p9a7/fs = -0.0034019999999999996x +0.999995y
p9a7/ss = -0.999995x -0.0034019999999999996y
p9a7/corner_x = 63.4124128906248
p9a7/corner_y = -361.6553872314453
p9a7/coffset = 0.700742

p10a0/dim1 = 10
p10a0/dim2 = ss
p10a0/dim3 = fs
p10a0/min_fs = 0
p10a0/min_ss = 0
p10a0/max_fs = 127
p10a0/max_ss = 63
p10a0/fs = -0.001794x +0.999996y
p10a0/ss = -0.999996x -0.001794y
p10a0/corner_x = 525.3137128906249
p10a0/corner_y = -517.4093872314453
p10a0/coffset = 0.700942

p10a1/dim1 = 10
p10a1/dim2 = ss
p10a1/dim3 = fs
p10a1/min_fs = 0
p10a1/min_ss = 64
p10a1/max_fs = 127
p10a1/max_ss = 127
p10a1/fs = -0.001794x +0.999996y
p10a1/ss = -0.999996x -0.001794y
p10a1/corner_x = 459.3167128906248
p10a1/corner_y = -517.5353872314453
p10a1/coffset = 0.700942

p10a2/dim1 = 10
p10a2/dim2 = ss
p10a2/dim3 = fs
p10a2/min_fs = 0
p10a2/min_ss = 128
p10a2/max_fs = 127
p10a2/max_ss = 191
p10a2/fs = -0.001794x +0.999996y
p10a2/ss = -0.999996x -0.001794y
p10a2/corner_x = 393.3187128906248
p10a2/corner_y = -517.6663872314451
p10a2/coffset = 0.700942

p10a3/dim1 = 10
p10a3/dim2 = ss
p10a3/dim3 = fs
p10a3/min_fs = 0
p10a3/min_ss = 192
p10a3/max_fs = 127
p10a3/max_ss = 255
p10a3/fs = -0.001794x +0.999996y
p10a3/ss = -0.999996x -0.001794y
p10a3/corner_x = 327.3237128906248
p10a3/corner_y = -517.7913872314452
p10a3/coffset = 0.700942

p10a4/dim1 = 10
p10a4/dim2 = ss
p10a4/dim3 = fs
p10a4/min_fs = 0
p10a4/min_ss = 256
p10a4/max_fs = 127
p10a4/max_ss = 319
p10a4/fs = -0.001794x +0.999996y
p10a4/ss = -0.999996x -0.001794y
p10a4/corner_x = 261.32571289062474
p10a4/corner_y = -517.9173872314453
p10a4/coffset = 0.700942

p10a5/dim1 = 10
p10a5/dim2 = ss
p10a5/dim3 = fs
p10a5/min_fs = 0
p10a5/min_ss = 320
p10a5/max_fs = 127
p10a5/max_ss = 383
p10a5/fs = -0.001794x +0.999996y
p10a5/ss = -0.999996x -0.001794y
p10a5/corner_x = 195.32871289062476
p10a5/corner_y = -518.0473872314453
p10a5/coffset = 0.700942

p10a6/dim1 = 10
p10a6/dim2 = ss
p10a6/dim3 = fs
p10a6/min_fs = 0
p10a6/min_ss = 384
p10a6/max_fs = 127
p10a6/max_ss = 447
p10a6/fs = -0.001794x +0.999996y
p10a6/ss = -0.999996x -0.001794y
p10a6/corner_x = 129.33271289062483
p10a6/corner_y = -518.1753872314453
p10a6/coffset = 0.700942

p10a7/dim1 = 10
p10a7/dim2 = ss
p10a7/dim3 = fs
p10a7/min_fs = 0
p10a7/min_ss = 448
p10a7/max_fs = 127
p10a7/max_ss = 511
p10a7/fs = -0.001794x +0.999996y
p10a7/ss = -0.999996x -0.001794y
p10a7/corner_x = 63.33311289062479
p10a7/corner_y = -518.2983872314452
p10a7/coffset = 0.700942

p11a0/dim1 = 11
p11a0/dim2 = ss
p11a0/dim3 = fs
p11a0/min_fs = 0
p11a0/min_ss = 0
p11a0/max_fs = 127
p11a0/max_ss = 63
p11a0/fs = -0.002832x +0.999996y
p11a0/ss = -0.999996x -0.002832y
p11a0/corner_x = 525.9947128906246
p11a0/corner_y = -674.3903872314452
p11a0/coffset = 0.701051

p11a1/dim1 = 11
p11a1/dim2 = ss
p11a1/dim3 = fs
p11a1/min_fs = 0
p11a1/min_ss = 64
p11a1/max_fs = 127
p11a1/max_ss = 127
p11a1/fs = -0.002832x +0.999996y
p11a1/ss = -0.999996x -0.002832y
p11a1/corner_x = 459.9957128906247
p11a1/corner_y = -674.5813872314451
p11a1/coffset = 0.701051

p11a2/dim1 = 11
p11a2/dim2 = ss
p11a2/dim3 = fs
p11a2/min_fs = 0
p11a2/min_ss = 128
p11a2/max_fs = 127
p11a2/max_ss = 191
p11a2/fs = -0.002832x +0.999996y
p11a2/ss = -0.999996x -0.002832y
p11a2/corner_x = 393.99871289062474
p11a2/corner_y = -674.7773872314452
p11a2/coffset = 0.701051

p11a3/dim1 = 11
p11a3/dim2 = ss
p11a3/dim3 = fs
p11a3/min_fs = 0
p11a3/min_ss = 192
p11a3/max_fs = 127
p11a3/max_ss = 255
p11a3/fs = -0.002832x +0.999996y
p11a3/ss = -0.999996x -0.002832y
p11a3/corner_x = 328.0017128906248
p11a3/corner_y = -674.9693872314452
p11a3/coffset = 0.701051

p11a4/dim1 = 11
p11a4/dim2 = ss
p11a4/dim3 = fs
p11a4/min_fs = 0
p11a4/min_ss = 256
p11a4/max_fs = 127
p11a4/max_ss = 319
p11a4/fs = -0.002832x +0.999996y
p11a4/ss = -0.999996x -0.002832y
p11a4/corner_x = 262.00771289062476
p11a4/corner_y = -675.1623872314451
p11a4/coffset = 0.701051

p11a5/dim1 = 11
p11a5/dim2 = ss
p11a5/dim3 = fs
p11a5/min_fs = 0
p11a5/min_ss = 320
p11a5/max_fs = 127
p11a5/max_ss = 383
p11a5/fs = -0.002832x +0.999996y
p11a5/ss = -0.999996x -0.002832y
p11a5/corner_x = 196.00971289062483
p11a5/corner_y = -675.3563872314451
p11a5/coffset = 0.701051

p11a6/dim1 = 11
p11a6/dim2 = ss
p11a6/dim3 = fs
p11a6/min_fs = 0
p11a6/min_ss = 384
p11a6/max_fs = 127
p11a6/max_ss = 447
p11a6/fs = -0.002832x +0.999996y
p11a6/ss = -0.999996x -0.002832y
p11a6/corner_x = 130.01471289062482
p11a6/corner_y = -675.5463872314451
p11a6/coffset = 0.701051

p11a7/dim1 = 11
p11a7/dim2 = ss
p11a7/dim3 = fs
p11a7/min_fs = 0
p11a7/min_ss = 448
p11a7/max_fs = 127
p11a7/max_ss = 511
p11a7/fs = -0.002832x +0.999996y
p11a7/ss = -0.999996x -0.002832y
p11a7/corner_x = 64.0130128906248
p11a7/corner_y = -675.7403872314453
p11a7/coffset = 0.701051

p12a0/dim1 = 12
p12a0/dim2 = ss
p12a0/dim3 = fs
p12a0/min_fs = 0
p12a0/min_ss = 0
p12a0/max_fs = 127
p12a0/max_ss = 63
p12a0/fs = -0.004863x +0.9999889999999999y
p12a0/ss = -0.9999889999999999x -0.004863y
p12a0/corner_x = 549.4152438659668
p12a0/corner_y = 434.9466091262818
p12a0/coffset = 0.700575

p12a1/dim1 = 12
p12a1/dim2 = ss
p12a1/dim3 = fs
p12a1/min_fs = 0
p12a1/min_ss = 64
p12a1/max_fs = 127
p12a1/max_ss = 127
p12a1/fs = -0.004863x +0.9999889999999999y
p12a1/ss = -0.9999889999999999x -0.004863y
p12a1/corner_x = 483.41624386596675
p12a1/corner_y = 434.5826091262818
p12a1/coffset = 0.700575

p12a2/dim1 = 12
p12a2/dim2 = ss
p12a2/dim3 = fs
p12a2/min_fs = 0
p12a2/min_ss = 128
p12a2/max_fs = 127
p12a2/max_ss = 191
p12a2/fs = -0.004863x +0.9999889999999999y
p12a2/ss = -0.9999889999999999x -0.004863y
p12a2/corner_x = 417.41924386596673
p12a2/corner_y = 434.21360912628165
p12a2/coffset = 0.700575

p12a3/dim1 = 12
p12a3/dim2 = ss
p12a3/dim3 = fs
p12a3/min_fs = 0
p12a3/min_ss = 192
p12a3/max_fs = 127
p12a3/max_ss = 255
p12a3/fs = -0.004863x +0.9999889999999999y
p12a3/ss = -0.9999889999999999x -0.004863y
p12a3/corner_x = 351.41624386596663
p12a3/corner_y = 433.98360912628164
p12a3/coffset = 0.700575

p12a4/dim1 = 12
p12a4/dim2 = ss
p12a4/dim3 = fs
p12a4/min_fs = 0
p12a4/min_ss = 256
p12a4/max_fs = 127
p12a4/max_ss = 319
p12a4/fs = -0.004863x +0.9999889999999999y
p12a4/ss = -0.9999889999999999x -0.004863y
p12a4/corner_x = 285.4162438659668
p12a4/corner_y = 433.66560912628177
p12a4/coffset = 0.700575

p12a5/dim1 = 12
p12a5/dim2 = ss
p12a5/dim3 = fs
p12a5/min_fs = 0
p12a5/min_ss = 320
p12a5/max_fs = 127
p12a5/max_ss = 383
p12a5/fs = -0.004863x +0.9999889999999999y
p12a5/ss = -0.9999889999999999x -0.004863y
p12a5/corner_x = 219.4192438659668
p12a5/corner_y = 433.34260912628173
p12a5/coffset = 0.700575

p12a6/dim1 = 12
p12a6/dim2 = ss
p12a6/dim3 = fs
p12a6/min_fs = 0
p12a6/min_ss = 384
p12a6/max_fs = 127
p12a6/max_ss = 447
p12a6/fs = -0.004863x +0.9999889999999999y
p12a6/ss = -0.9999889999999999x -0.004863y
p12a6/corner_x = 153.41924386596676
p12a6/corner_y = 433.0236091262817
p12a6/coffset = 0.700575

p12a7/dim1 = 12
p12a7/dim2 = ss
p12a7/dim3 = fs
p12a7/min_fs = 0
p12a7/min_ss = 448
p12a7/max_fs = 127
p12a7/max_ss = 511
p12a7/fs = -0.004863x +0.9999889999999999y
p12a7/ss = -0.9999889999999999x -0.004863y
p12a7/corner_x = 87.41834386596676
p12a7/corner_y = 432.70060912628173
p12a7/coffset = 0.700575

p13a0/dim1 = 13
p13a0/dim2 = ss
p13a0/dim3 = fs
p13a0/min_fs = 0
p13a0/min_ss = 0
p13a0/max_fs = 127
p13a0/max_ss = 63
p13a0/fs = -0.001143x +0.999999y
p13a0/ss = -0.999999x -0.001143y
p13a0/corner_x = 548.0732438659666
p13a0/corner_y = 276.3266091262817
p13a0/coffset = 0.700561

p13a1/dim1 = 13
p13a1/dim2 = ss
p13a1/dim3 = fs
p13a1/min_fs = 0
p13a1/min_ss = 64
p13a1/max_fs = 127
p13a1/max_ss = 127
p13a1/fs = -0.001143x +0.999999y
p13a1/ss = -0.999999x -0.001143y
p13a1/corner_x = 482.0742438659666
p13a1/corner_y = 276.24260912628176
p13a1/coffset = 0.700561

p13a2/dim1 = 13
p13a2/dim2 = ss
p13a2/dim3 = fs
p13a2/min_fs = 0
p13a2/min_ss = 128
p13a2/max_fs = 127
p13a2/max_ss = 191
p13a2/fs = -0.001143x +0.999999y
p13a2/ss = -0.999999x -0.001143y
p13a2/corner_x = 416.07924386596665
p13a2/corner_y = 276.1546091262818
p13a2/coffset = 0.700561

p13a3/dim1 = 13
p13a3/dim2 = ss
p13a3/dim3 = fs
p13a3/min_fs = 0
p13a3/min_ss = 192
p13a3/max_fs = 127
p13a3/max_ss = 255
p13a3/fs = -0.001143x +0.999999y
p13a3/ss = -0.999999x -0.001143y
p13a3/corner_x = 350.08024386596674
p13a3/corner_y = 276.0666091262817
p13a3/coffset = 0.700561

p13a4/dim1 = 13
p13a4/dim2 = ss
p13a4/dim3 = fs
p13a4/min_fs = 0
p13a4/min_ss = 256
p13a4/max_fs = 127
p13a4/max_ss = 319
p13a4/fs = -0.001143x +0.999999y
p13a4/ss = -0.999999x -0.001143y
p13a4/corner_x = 284.0822438659668
p13a4/corner_y = 275.98060912628165
p13a4/coffset = 0.700561

p13a5/dim1 = 13
p13a5/dim2 = ss
p13a5/dim3 = fs
p13a5/min_fs = 0
p13a5/min_ss = 320
p13a5/max_fs = 127
p13a5/max_ss = 383
p13a5/fs = -0.001143x +0.999999y
p13a5/ss = -0.999999x -0.001143y
p13a5/corner_x = 218.06724386596673
p13a5/corner_y = 275.89160912628165
p13a5/coffset = 0.700561

p13a6/dim1 = 13
p13a6/dim2 = ss
p13a6/dim3 = fs
p13a6/min_fs = 0
p13a6/min_ss = 384
p13a6/max_fs = 127
p13a6/max_ss = 447
p13a6/fs = -0.001143x +0.999999y
p13a6/ss = -0.999999x -0.001143y
p13a6/corner_x = 152.0682438659668
p13a6/corner_y = 275.80560912628175
p13a6/coffset = 0.700561

p13a7/dim1 = 13
p13a7/dim2 = ss
p13a7/dim3 = fs
p13a7/min_fs = 0
p13a7/min_ss = 448
p13a7/max_fs = 127
p13a7/max_ss = 511
p13a7/fs = -0.001143x +0.999999y
p13a7/ss = -0.999999x -0.001143y
p13a7/corner_x = 86.06744386596675
p13a7/corner_y = 275.7186091262817
p13a7/coffset = 0.700561

p14a0/dim1 = 14
p14a0/dim2 = ss
p14a0/dim3 = fs
p14a0/min_fs = 0
p14a0/min_ss = 0
p14a0/max_fs = 127
p14a0/max_ss = 63
p14a0/fs = -0.004884x +0.999988y
p14a0/ss = -0.999988x -0.004884y
p14a0/corner_x = 550.3832438659665
p14a0/corner_y = 120.71360912628174
p14a0/coffset = 0.700526

p14a1/dim1 = 14
p14a1/dim2 = ss
p14a1/dim3 = fs
p14a1/min_fs = 0
p14a1/min_ss = 64
p14a1/max_fs = 127
p14a1/max_ss = 127
p14a1/fs = -0.004884x +0.999988y
p14a1/ss = -0.999988x -0.004884y
p14a1/corner_x = 484.3872438659667
p14a1/corner_y = 120.38960912628175
p14a1/coffset = 0.700526

p14a2/dim1 = 14
p14a2/dim2 = ss
p14a2/dim3 = fs
p14a2/min_fs = 0
p14a2/min_ss = 128
p14a2/max_fs = 127
p14a2/max_ss = 191
p14a2/fs = -0.004884x +0.999988y
p14a2/ss = -0.999988x -0.004884y
p14a2/corner_x = 418.3902438659666
p14a2/corner_y = 120.06460912628174
p14a2/coffset = 0.700526

p14a3/dim1 = 14
p14a3/dim2 = ss
p14a3/dim3 = fs
p14a3/min_fs = 0
p14a3/min_ss = 192
p14a3/max_fs = 127
p14a3/max_ss = 255
p14a3/fs = -0.004884x +0.999988y
p14a3/ss = -0.999988x -0.004884y
p14a3/corner_x = 352.3952438659666
p14a3/corner_y = 119.73860912628172
p14a3/coffset = 0.700526

p14a4/dim1 = 14
p14a4/dim2 = ss
p14a4/dim3 = fs
p14a4/min_fs = 0
p14a4/min_ss = 256
p14a4/max_fs = 127
p14a4/max_ss = 319
p14a4/fs = -0.004884x +0.999988y
p14a4/ss = -0.999988x -0.004884y
p14a4/corner_x = 286.39624386596677
p14a4/corner_y = 119.41460912628172
p14a4/coffset = 0.700526

p14a5/dim1 = 14
p14a5/dim2 = ss
p14a5/dim3 = fs
p14a5/min_fs = 0
p14a5/min_ss = 320
p14a5/max_fs = 127
p14a5/max_ss = 383
p14a5/fs = -0.004884x +0.999988y
p14a5/ss = -0.999988x -0.004884y
p14a5/corner_x = 220.38424386596674
p14a5/corner_y = 119.0866091262817
p14a5/coffset = 0.700526

p14a6/dim1 = 14
p14a6/dim2 = ss
p14a6/dim3 = fs
p14a6/min_fs = 0
p14a6/min_ss = 384
p14a6/max_fs = 127
p14a6/max_ss = 447
p14a6/fs = -0.004884x +0.999988y
p14a6/ss = -0.999988x -0.004884y
p14a6/corner_x = 154.38324386596673
p14a6/corner_y = 118.76160912628174
p14a6/coffset = 0.700526

p14a7/dim1 = 14
p14a7/dim2 = ss
p14a7/dim3 = fs
p14a7/min_fs = 0
p14a7/min_ss = 448
p14a7/max_fs = 127
p14a7/max_ss = 511
p14a7/fs = -0.004884x +0.999988y
p14a7/ss = -0.999988x -0.004884y
p14a7/corner_x = 88.38144386596679
p14a7/corner_y = 118.43660912628171
p14a7/coffset = 0.700526

p15a0/dim1 = 15
p15a0/dim2 = ss
p15a0/dim3 = fs
p15a0/min_fs = 0
p15a0/min_ss = 0
p15a0/max_fs = 127
p15a0/max_ss = 63
p15a0/fs = -0.001961x +1.000001y
p15a0/ss = -1.000001x -0.001961y
p15a0/corner_x = 550.3812438659667
p15a0/corner_y = -35.92479087371826
p15a0/coffset = 0.7004579999999999

p15a1/dim1 = 15
p15a1/dim2 = ss
p15a1/dim3 = fs
p15a1/min_fs = 0
p15a1/min_ss = 64
p15a1/max_fs = 127
p15a1/max_ss = 127
p15a1/fs = -0.001961x +1.000001y
p15a1/ss = -1.000001x -0.001961y
p15a1/corner_x = 484.38324386596673
p15a1/corner_y = -36.060290873718266
p15a1/coffset = 0.7004579999999999

p15a2/dim1 = 15
p15a2/dim2 = ss
p15a2/dim3 = fs
p15a2/min_fs = 0
p15a2/min_ss = 128
p15a2/max_fs = 127
p15a2/max_ss = 191
p15a2/fs = -0.001961x +1.000001y
p15a2/ss = -1.000001x -0.001961y
p15a2/corner_x = 418.3852438659667
p15a2/corner_y = -36.19619087371826
p15a2/coffset = 0.7004579999999999

p15a3/dim1 = 15
p15a3/dim2 = ss
p15a3/dim3 = fs
p15a3/min_fs = 0
p15a3/min_ss = 192
p15a3/max_fs = 127
p15a3/max_ss = 255
p15a3/fs = -0.001961x +1.000001y
p15a3/ss = -1.000001x -0.001961y
p15a3/corner_x = 352.3882438659667
p15a3/corner_y = -36.33149087371825
p15a3/coffset = 0.7004579999999999

p15a4/dim1 = 15
p15a4/dim2 = ss
p15a4/dim3 = fs
p15a4/min_fs = 0
p15a4/min_ss = 256
p15a4/max_fs = 127
p15a4/max_ss = 319
p15a4/fs = -0.001961x +1.000001y
p15a4/ss = -1.000001x -0.001961y
p15a4/corner_x = 286.3882438659668
p15a4/corner_y = -36.467190873718266
p15a4/coffset = 0.7004579999999999

p15a5/dim1 = 15
p15a5/dim2 = ss
p15a5/dim3 = fs
p15a5/min_fs = 0
p15a5/min_ss = 320
p15a5/max_fs = 127
p15a5/max_ss = 383
p15a5/fs = -0.001961x +1.000001y
p15a5/ss = -1.000001x -0.001961y
p15a5/corner_x = 220.36924386596672
p15a5/corner_y = -36.60269087371826
p15a5/coffset = 0.7004579999999999

p15a6/dim1 = 15
p15a6/dim2 = ss
p15a6/dim3 = fs
p15a6/min_fs = 0
p15a6/min_ss = 384
p15a6/max_fs = 127
p15a6/max_ss = 447
p15a6/fs = -0.001961x +1.000001y
p15a6/ss = -1.000001x -0.001961y
p15a6/corner_x = 154.37124386596673
p15a6/corner_y = -36.73879087371825
p15a6/coffset = 0.7004579999999999

p15a7/dim1 = 15
p15a7/dim2 = ss
p15a7/dim3 = fs
p15a7/min_fs = 0
p15a7/min_ss = 448
p15a7/max_fs = 127
p15a7/max_ss = 511
p15a7/fs = -0.001961x +1.000001y
p15a7/ss = -1.000001x -0.001961y
p15a7/corner_x = 88.37004386596676
p15a7/corner_y = -36.874090873718266
p15a7/coffset = 0.7004579999999999