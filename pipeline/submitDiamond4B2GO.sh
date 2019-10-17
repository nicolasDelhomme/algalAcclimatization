#!/bin/bash -l

set -ex

proj=u2019005
in=../data/TransDecoder
out=../data/diamond
inx=../uniref/indices/diamond/uniref90.dmnd
mapping=../taxonomy/prot.accession2taxid.gz
nodes=../taxonomy/nodes.dmp
names=../taxonomy/names.dmp

if [ ! -d $out ]; then
	mkdir -p $out
fi

# run
sbatch -n 32 -o $out/diamond4b2g.out \
-e $out/diamond4b2g.err ../UPSCb-common/pipeline/runDiamond4B2G0.sh \
-i $mapping -n $nodes -t $names \
-p 32 blastp $in/Trinity.fasta.transdecoder.pep $inx $out
