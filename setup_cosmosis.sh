conda activate /project/chihway/envs/cosmosis3


if [ "$USER" == "dhayaa" ]
then
    export REPO_DIR=/home/dhayaa/DECADE/CosmicShearCosmosis/
    export OUTDIR=/project/chihway/dhayaa/DECADE/cosmosis/
    export COSMOSIS_DIR=/home/dhayaa/DECADE/cosmosis-standard-library/
fi

if [ "$USER" == "chihway" ]
then
    export REPO_DIR=/project/chihway/chihway/CosmicShearCosmosis/
    export OUTDIR=/project/chihway/chihway/CosmicShearCosmosis/analysis/
    export COSMOSIS_DIR=/home/chihway/cosmosis-standard-library/
fi

if [ "$USER" == "nchicoine" ]
then
    export REPO_DIR=/home/nchicoine/delve/CosmicShearCosmosis/
    export OUTDIR=/project/chihway/nchicoine/DECADE/cosmosis/
    export COSMOSIS_DIR=/home/nchicoine/delve/cosmosis-standard-library/
fi
