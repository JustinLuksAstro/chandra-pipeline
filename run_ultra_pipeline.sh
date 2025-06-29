#!/bin/bash

# Usage: ./run_ultra_pipeline.sh <OBS_ID>

obsid=$1

if [ -z "$obsid" ]; then
  echo "Usage: $0 <OBS_ID>"
  exit 1
fi

# Download and reprocess
download_chandra_obsid $obsid
cd $obsid
chandra_repro . repro/
cd repro

# Launch DS9 so user can inspect source and get coordinates
ds9 acisf${obsid}_repro_evt2.fits &

# Wait for user to input coordinates
echo "Please inspect the file in DS9 and enter the coordinates when ready."
read -p "Enter pixel X center: " xcenter
read -p "Enter pixel Y center: " ycenter
read -p "Enter source radius (in pixels): " srcrad
read -p "Enter background X center: " x_bkg
read -p "Enter background Y center: " y_bkg
read -p "Enter background radius (in pixels): " bkgrad
read -p "Enter RA (HH:MM:SS): " ra_hms
read -p "Enter Dec (\u00b1DD:MM:SS): " dec_dms

# Convert RA/Dec to decimal
ra_dec=$(python3 -c "from astropy.coordinates import SkyCoord; c = SkyCoord('$ra_hms $dec_dms', unit=('hourangle', 'deg')); print(f'{c.ra.deg} {c.dec.deg}')")
ra=$(echo $ra_dec | awk '{print $1}')
dec=$(echo $ra_dec | awk '{print $2}')

# Run specextract
specextract \
  infile="acisf${obsid}_repro_evt2.fits[sky=circle(${xcenter},${ycenter},${srcrad})]" \
  outroot=${obsid} \
  bkgfile="acisf${obsid}_repro_evt2.fits[sky=circle(${x_bkg},${y_bkg},${bkgrad})]" \
  bkgresp=no weight=no correctpsf=yes grouptype=NONE verbose=0 mode=h

# Launch Sherpa
echo "Now launching Sherpa. Please run your fitting commands manually and save the spectrum as ${obsid}.dat"
sherpa

# Run simulate_psf
pset simulate_psf infile=acisf${obsid}_repro_evt2.fits
pset simulate_psf outroot=${obsid}_marxsim
pset simulate_psf ra=$ra
pset simulate_psf dec=$dec
pset simulate_psf spectrum=${obsid}.dat
pset simulate_psf blur=0.19
pset simulate_psf readout_streak=yes
pset simulate_psf pileup=no
pset simulate_psf ideal=no
pset simulate_psf extended=no

simulate_psf

# Deconvolution (1.0 binning)
dmcopy "acisf${obsid}_repro_evt2.fits[bin x=$(echo "$xcenter-6" | bc):$(echo "$xcenter+6" | bc):1,y=$(echo "$ycenter-6" | bc):$(echo "$ycenter+6" | bc):1]" obs_image.fits
arestore obs_image.fits ${obsid}_marxsim_psf_image.fits deconvolved_image.fits method=lucy numiter=100

# Deconvolution (0.25 binning)
dmcopy "${obsid}_marxsim_projrays.fits[bin x=$(echo "$xcenter-6" | bc):$(echo "$xcenter+6" | bc):0.25,y=$(echo "$ycenter-6" | bc):$(echo "$ycenter+6" | bc):0.25]" ${obsid}_marxsim_projrays_25.fits
dmcopy "acisf${obsid}_repro_evt2.fits[bin x=$(echo "$xcenter-6" | bc):$(echo "$xcenter+6" | bc):0.25,y=$(echo "$ycenter-6" | bc):$(echo "$ycenter+6" | bc):0.25]" obs_image_25.fits
arestore obs_image_25.fits ${obsid}_marxsim_projrays_25.fits deconvolved_image_bin.fits method=lucy numiter=100

# PSF Hook Feature
roll=$(dmkeypar acisf${obsid}_repro_evt2.fits ROLL_NOM echo+)
make_psf_asymmetry_region acisf${obsid}_repro_evt2.fits hook.reg $xcenter $ycenter format=ds9
ds9 acisf${obsid}_repro_evt2.fits -region hook.reg &

