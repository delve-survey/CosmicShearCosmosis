#! /bin/sh

# This script 
# Step 1. replaces the data vector in the simulated 2pt files with the treecorr measurements
# Step 2. blinds the data vector
# Step 3. removes the unblinded data vector
#
# To unblind, we would just comment out Step 2 and 3

sim_dv='/project/chihway/chihway/CosmicShearCosmosis/datavectors/baseline_Jan12th24.fits'
data_dir='/project/chihway/chihway/shearcat/shear_catalog/datavec/'
unblinded_dv='/project/chihway/chihway/CosmicShearCosmosis/datavectors/data_04052024.fits'
blinded_dv='/project/chihway/chihway/CosmicShearCosmosis/datavectors/data_04052024_BLINDED.fits'

# it is a little strange that it seems like we need to be in a large node to do this
#cd /project/chihway/chihway/CosmicShearCosmosis/
#source setup_cosmosis.sh

# Step 1
echo 'step 1'
python MakeBlindedDV.py $sim_dv $data_dir $unblinded_dv

# Step 2
echo 'step 2'
cd ../blinding
python -m blind_2pt_cosmosis -i blind_delve_fiducial.ini -u $unblinded_dv

# Step 3
echo 'step 3'
rm $unblinded_dv

