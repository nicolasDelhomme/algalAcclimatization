#!/bin/bash

# we need to be safe
# if there  is an error, we want to stop the progress
# we set -e (stop on error) for the whole of the script
set -e


# define variables
in=~/algal-acclimatization/raw
inx=~/algal-acclimatization/indices/salmon/Trinity.inx
out=~/algal-acclimatization/Salmon
account=u2019005
mail=amit.bajhaiya@umu.se

# create output directory
if [ ! -d $out ]; then
  mkdir -p $out
fi

# list all files in the in directory and run FastQC
for f in $(find $in -name "*_1.fq.gz"); do
  
  fnam=$(basename $f)
  prefix=${fnam/_1.fq.gz/}
  
  # submit to the queueing system using sbatch
  sbatch -A $account -e $out/${prefix}.err -o $out/${prefix}.out -J $prefix \
  --mail-user=$mail $UPSCb/pipeline/runSalmon.sh $inx  $f $in/${prefix}_2.fq.gz $out
  # for debugging, you can also use "echo command" to print the "command"
  # echo sbatch runFastQC.sh $out $f

done
