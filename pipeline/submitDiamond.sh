#!/bin/bash -l

set -ex

proj=u2019005
in=/mnt/picea/projects/algae/cfunk/algal-acclimatization/TransDecoder
out=/mnt/picea/projects/algae/cfunk/algal-acclimatization/diamond
inx=/mnt/picea/storage/reference/UniRef90/201908/indices/diamond/uniref90.dmnd
mapping=/mnt/picea/storage/reference/Taxonomy/20190825/prot.accession2taxid.gz
nodes=/mnt/picea/storage/reference/Taxonomy/20190825/nodes.dmp
names=/mnt/picea/storage/reference/Taxonomy/20190825/names.dmp

if [ ! -d $out ]; then
	mkdir -p $out
fi

# run
sbatch -n 64 -t 2:00:00 -o $out/diamond.out -e $out/diamond.err $UPSCb/pipeline/runDiamond.sh \
-m -s -i $mapping -n $nodes -t $names \
-p 64 blastp $in/Trinity.fasta.transdecoder.pep $inx $out
