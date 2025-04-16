#!/bin/bash
# XX.fas --> XX.fq --> XX.fq.gz --> XX_S00_L001_R1_001.fastq.gz

CMDNAME=`basename $0`

# parameter check
if [ "$#" -lt 1 ]; then
  echo "Usage: bash $CMDNAME list"
  echo '  $1  Strain List'
  echo '  $2  Extension (.fasta by default )'
  exit 1
fi

LIST=$1
if [ $2 ]; then
  EXT=$2
else
  EXT=.fasta
fi             
source activate art
for strain in `cat $LIST`; do
  strain=${strain}
  art_illumina -ss MSv3 -i ${strain}${EXT} -p -l 250 -f 100 -m 600 -s 10 --noALN -o ${strain}_
  gzip ${strain}_1.fq
  gzip ${strain}_2.fq
  mv ${strain}_1.fq.gz ${strain}_S00_L001_R1_001.fastq.gz
  mv ${strain}_2.fq.gz ${strain}_S00_L001_R2_001.fastq.gz
done
conda deactivate
