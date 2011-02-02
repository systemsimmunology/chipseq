#!/bin/bash
#
# For a chipseq BED file, computes overlap with annotation file with pybed.py
# then find the feature nearest to TSS with GenomicFeatureChIPSeqOverlap.py
# and finally, includes additional annotations.
#
# Output:$bedID-geneoverlap.tsv
#
WRONGARGS=1
if [ $# != 4 ]
then
  echo "Usage: `basename $0` <ChIPSeq BED file> <feature BED file> <bridged feature BED file> <Feature Prefix>" >&2
  exit $WRONGARGS
fi

PY3="/Library/Frameworks/Python.framework/Versions/3.1/bin/python3.1"

#featurefile="../annotation/GenomeAnnotationsM37.bed"
#featureprefix="ENSMUST"

#featurefile="../annotation/refGene.mouse.bed"
#featureprefix="NM_"

chipseqbedfile=$1
featurefile=$2
featurefileBridged=$3
featureprefix=$4

startdir=$PWD
bedID=`basename $1 .bed`

## Go to bedtools and run pybed.py
pushd ~/chipseq/bedtools >& /dev/null
$PY3 pybed.py nostrand -v 0 -f $1 $featurefileBridged > temp.bed
cp -p temp.bed $startdir/. 
popd >& /dev/null

## Derive parsed overlap file
ofile=$bedID.neartss.tsv

echo "Creating:" $ofile
~/chipseq/utils/TransformPyBedFile.py temp.bed $featurefile $featureprefix $bedID > simplified.bed
~/chipseq/utils/AddScore.py simplified.bed $chipseqbedfile > simpwithscore
~/chipseq/utils/ChIPSeqTSSCentered.py simpwithscore $featurefile > $ofile

##rm -f temp.bed

##
## Include mappings to Genes (Entrez genes or ENSMUSGs)
## 

if [ $featureprefix = "NM_" ]; then
    awk '{print $1}' $ofile > refseqs
    ~/bin/reMap.py refseqs  1 ~/data/ncbi/gene2refseqSimplified_NM_mouse 2 1 > enids
    ~/bin/reMap.py enids 1 ~/data/ncbi/gene_info_simplified_mouse 1 2 > gnames
    awk -F "\t" '{OFS="\t"; print $2,$3,$4,$5,$6,$7,$8}' $ofile > t1
    paste refseqs enids gnames t1 > t2
    rm -f refseqs enids gnames
    newheader="Genome Feature\tEntrez ID\tSymbol\tDistance to Midpoint\tSegment Length\tDistance to Segment Boundary\tSegment Score "
#    newheader="Genome Feature\tEntrez ID\tSymbol\tFractional Overlap\tLength of Overlap\tLength of Genome Feature\tFeature Chromosome\tFeature Start\tFeature End\tFeature Strand"
    tail +2 t2 > t3
    echo -e $newheader > t4
    finalfile=$bedID.tss.geneannot.tsv
    cat t4 t3 > $finalfile
    rm -f t1 t2 t3 t4
fi
