#!/bin/bash -l

## error and verbose
set -ex

# modules
module load bioinfo-tools blast

## default args
proj=u2019005
mail="nicolas.delhomme@slu.se"
in=$(realpath ../data/blastn/tmp/Trinity.fasta)
out=../data/blastn
inx=../NCBI/nt

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
sbatch -A $proj -t 2-00:00:00 --mail-user $mail -e $out/nt.err -o $out/nt.out --array=1-400 \
../UPSCb-common/pipeline/runBlastPlusXMLasArray.sh blastn $in $inx $out 
