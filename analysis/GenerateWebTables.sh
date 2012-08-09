#!/bin/bash

##
head -1 FeatMatRefSeq.tsv > header 

##Full FeatureMatrix, filtered
##awk '{FS="\t" ; if( ($27!="NA") && ($27!="nPnRnI") ) print }' FeatMatRefSeq.tsv > FeatMatRefSeq.2.tsv
##awk '{FS="\t" ; if($11 > 4) print }' FeatMatRefSeq.2.tsv > FeatMatRefSeq.3.tsv

tags=(nPRI PnRI PRnI PnRnI nPRnI nPnRI)

for tag in "${tags[@]}"; do
    echo $tag

## Create $tag web folder and populate 
    mkdir $CR9/PublicDatasets/EpiGenomeLandscape/$tag
    pushd $CR9/PublicDatasets/EpiGenomeLandscape/$tag >& /dev/null
    cp -p ../PRIdev2/annotations_grid.json . 
    sed "s/PRI/$tag/g" annotations_grid.json > tempfile 
    mv -f tempfile annotations_grid.json 
    rm -f tempfile
    ln -s /titan/cancerregulome9/workspaces/users/vthorsson/EpiGenomeLandscape/KinPlots images ## require unix form of filesystem link
    popd >& /dev/null
## Create Feature matrix $tag and copy to web folder 
    awk -v t=$tag '{FS="\t" ; if ($27==t) print }' FeatMatRefSeq.tsv > tempfile ## Pass shell variables to awk http://www.cyberciti.biz/faq/linux-unix-appleosx-bsd-bash-passing-variables-to-awk/
    cat header tempfile > FeatMatRefSeq.$tag.tsv 
    rm -f tempfile
    cp -pf FeatMatRefSeq.$tag.tsv  $CR9/PublicDatasets/EpiGenomeLandscape/$tag/.
    pushd $CR9/PublicDatasets/EpiGenomeLandscape/$tag >& /dev/null
    ln -s FeatMatRefSeq.$tag.tsv annotations.tsv
## Update description with gene correct count
    n=`wc -l annotations.tsv | awk '{print $1}'`
    let "n -= 1" ## no of genes
    sed "s/118/$n/g" annotations_grid.json > tempfile
    mv -f tempfile annotations_grid.json
    rm -f tempfile
    popd >& /dev/null

done

rm -f header