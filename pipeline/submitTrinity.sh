#!/bin/bash -l

## be verbose (-x) and stop on error (-e)
set -ex

## define a function
usage () {
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
    exit 1
}

## variables
in=/mnt/picea/projects/algae/cfunk/algal-acclimatization/raw
out=/mnt/picea/projects/algae/cfunk/algal-acclimatization/trinity
mail=amit.bajhaiya@umu.se
account=u2019005
MEM=300G
CPU=32

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
sbatch -A $account --mem=$MEM --mail-user=$mail -e $out/trinity.err -o $out/trinity.out \
-p nolimit -n $CPU -t 10-00:00:00 \
$UPSCb/pipeline/runTrinity.sh -r -p $CPU -m $MEM $out \
$(find $in -name "*_1.fq.gz" | sort | paste -s --delimiters=,) \
$(find $in -name "*_2.fq.gz" | sort | paste -s --delimiters=,)
