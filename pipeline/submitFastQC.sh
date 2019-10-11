#!/bin/bash

# we need to be safe
# if there  is an error, we want to stop the progress
# we set -e (stop on error) for the whole of the script
set -e

# we need the tools
# use the "module" tool manager
# use `module avail` to list all tools
module load bioinfo-tools FastQC

# define variables
in=~/algal-acclimatization/raw
out=~/algal-acclimatization/fastqc/raw
account=u2019005
mail=amit.bajhaiya@umu.se

# create output directory
if [ ! -d $out ]; then
  mkdir -p $out
fi

# list all files in the in directory and run FastQC
for f in $(find $in -name "*.fq.gz"); do
  
  fnam=$(basename $f)
  prefix=${fnam/.fq.gz/}
  
  # submit to the queueing system using sbatch
  sbatch -A $account -e $out/${prefix}.err -o $out/${prefix}.out -J $prefix \
  --mail-user=$mail $UPSCb/pipeline/runFastQC.sh $out $f
  
  # for debugging, you can also use "echo command" to print the "command"
  # echo sbatch runFastQC.sh $out $f
done
