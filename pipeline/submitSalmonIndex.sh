#!/bin/bash

# we need to be safe
# if there  is an error, we want to stop the progress
# we set -e (stop on error) for the whole of the script
set -e

# define variables
in=~/algal-acclimatization/trinity/Trinity.fasta
out=~/algal-acclimatization/indices/salmon
account=u2019005
mail=amit.bajhaiya@umu.se

# create output directory
if [ ! -d $out ]; then
  mkdir -p $out
fi

# submit to the queueing system using sbatch
sbatch -A $account -e $out/index.err -o $out/index.out -J salmon-index \
--mail-user=$mail $UPSCb/pipeline/runSalmonIndex.sh $in $out/Trinity.inx
