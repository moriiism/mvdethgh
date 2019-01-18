#!/bin/sh

cd /home/morii/work/urakawa/
mkdir 19011801
cd 19011801

# get psf

mkdir setup
mkdir log

cat << EOF > setup/data.list
# 2d-fits  time(sec)
../data/tomoe_2010wc9/2010wc9_1.fits    0.25
../data/tomoe_2010wc9/2010wc9_2.fits    0.75
../data/tomoe_2010wc9/2010wc9_3.fits    1.25
../data/tomoe_2010wc9/2010wc9_4.fits    1.75
../data/tomoe_2010wc9/2010wc9_5.fits    2.25
../data/tomoe_2010wc9/2010wc9_6.fits    2.75
../data/tomoe_2010wc9/2010wc9_7.fits    3.25
../data/tomoe_2010wc9/2010wc9_8.fits    3.75
../data/tomoe_2010wc9/2010wc9_9.fits    4.25
../data/tomoe_2010wc9/2010wc9_10.fits   4.75
../data/tomoe_2010wc9/2010wc9_11.fits   5.25
../data/tomoe_2010wc9/2010wc9_12.fits   5.75
EOF

data_list=setup/data.list
subimg_dat="none"
val_smooth=5.0
nbin_kernel_half=10
nbin_psf_half=10
outdir=psf
outfile_head=psf

/home/morii/work/github/moriiism/mvdethgh/mkpsf/mkpsf \
$data_list \
$subimg_dat \
$val_smooth \
$nbin_kernel_half \
$nbin_psf_half \
$outdir \
$outfile_head > log/mkpsf.log 2>&1

# pre-process

data_list=setup/data.list
subimg_dat="none"
psf_dat=psf/psf_psf.dat
outdir=preproc
outfile_head=preproc

/home/morii/work/github/moriiism/mvdethgh/preproc/preproc \
$data_list \
$subimg_dat \
$psf_dat \
$outdir \
$outfile_head > log/preproc.log 2>&1

# flat

cat << EOF > setup/mvobj_conv.list
# 2d-fits  time(sec)
preproc/preproc_conv_00.fits  0.25
preproc/preproc_conv_01.fits  0.75
preproc/preproc_conv_02.fits  1.25
preproc/preproc_conv_03.fits  1.75
preproc/preproc_conv_04.fits  2.25
preproc/preproc_conv_05.fits  2.75
preproc/preproc_conv_06.fits  3.25
preproc/preproc_conv_07.fits  3.75
preproc/preproc_conv_08.fits  4.25
preproc/preproc_conv_09.fits  4.75
preproc/preproc_conv_10.fits  5.25
preproc/preproc_conv_11.fits  5.75
EOF

# time hist 
cat << EOF > setup/time.dat
# nbin lo up
12  0.0  6.0
EOF

# velocity range
# theta = atan(pixel/sec)
cat << EOF > setup/vel.dat
# velocity(pixel/sec)
# nbin lo up
6  0.75  3.75
EOF

# resolution for internal parameters
# number of bins for (rho, phi, psi)
cat << EOF > setup/res.dat
# nrho, nphi, npsi
1000  200  200
EOF

data_list=setup/mvobj_conv.list
subimg_dat="none"
time_dat=setup/time.dat
vel_dat=setup/vel.dat
res_dat=setup/res.dat
outdir=det
outfile_head=flat

/home/morii/work/github/moriiism/mvdethgh/mkflat/mkflat \
$data_list \
$subimg_dat \
$time_dat \
$vel_dat \
$res_dat \
$outdir \
$outfile_head > log/mkflat.log 2>&1

# det

data_list=setup/mvobj_conv.list
subimg_dat="none"
time_dat=setup/time.dat
vel_dat=setup/vel.dat
res_dat=setup/res.dat
sig=5
ndet=100
outdir=det
outfile_head=det

/home/morii/work/github/moriiism/mvdethgh/mvdethgh/mvdethgh4d \
$data_list \
$subimg_dat \
$time_dat \
$vel_dat \
$res_dat \
$sig \
$ndet \
$outdir \
$outfile_head > log/mvdethgh4d.log 2>&1
