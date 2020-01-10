#!/bin/bash
./make_webs.sh

OUTDIR="../data/c_webs/" ## ouput directory

mkdir -p $OUTDIR

for S in 400 ## total number of species 
do 

for SB in 50 100 150 200 ## number of basal species 
do 

SC=$(($S-$SB)) ## number of consumer species 
SEED=$(($SB*$SC)) ## set a seed to generate different webs (mandatory, if multiple webs with the same set of species are generated) 

export OUTDIR
export SEED
export SB
export SC

./web $OUTDIR $SEED $SB $SC

done
done  
