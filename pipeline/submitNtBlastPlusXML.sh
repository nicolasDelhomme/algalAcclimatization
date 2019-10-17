#!/bin/bash -l

## error and verbose
set -ex

# modules
module load bioinfo-tools blast

## default args
proj=u2019005
mail="nicolas.delhomme@slu.se"
in=../data/trinity/Trinity.fasta
out=../data/blastn
inx=../NCBI/nt

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
sbatch -A $proj --mail-user $mail -e $out/nt.err -o $out/nt.out \
../UPSCb-common/pipeline/runBlastPlusXML.sh blastn $in $inx $out 
