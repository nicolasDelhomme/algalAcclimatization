#!/bin/bash -l

## be verbose and print
set -ex

proj=u2019005
mail=amit.bajhaiya@umu.se

## define a function
usage () {
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
}

## source functions
source $UPSCb/src/bash/functions.sh

## process the argument
in=~/algal-acclimatization/trinity/Trinity.fasta
out=~/algal-acclimatization/TransDecoder

## check vars
if [ -z $UPSCb ]; then
    echo "The UPSCb var needs to be set."
    usage
fi

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## execute
sbatch --mail-user=$mail -A $proj -J transdecoder \
-e $out/transdecoder.err -o $out/transdecoder.out \
$UPSCb/pipeline/runTrinityTransDecoder.sh $in $out



