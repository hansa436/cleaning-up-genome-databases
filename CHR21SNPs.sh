#!/bin/bash
# create a sliding windoe for coodorinates for all on chr21
#STEP is the amount moves per round
#WS= window size (bp)
#CHROM= chromosome 21
#START= Start position
#END= End position

STEP=20
WS=75
CHROM=21
START=5010000
END=48129895

#loop
rm OutputCHR21SNPs.csv
#set -xv
printf '%s\n' 'Start, SNPS, PHASTCONS' >> OutputCHR21SNPs.csv
while [ $START -lt $END ]
        do
        PHASTCONS=$(bigWigSummary -type=mean hg38.phastCons100way.bw chr21 $START $(( $START + $WS )) 1 2>&1)
        #if [[ $( echo $PHASTCONS | tr -cd "e" | wc -c ) == "1" ]] ; then PHASTCONS=0 ; else : ; fi
        [[ "$PHASTCONS" == *"e"* ]] && PHASTCONS=0
        #bc is needed to fix floating numbers
                if (( $( echo "$PHASTCONS > 0.4999999" | bc -l ) ))
                then 
		SNPS=$(tabix -h gnomad.genomes.v3.1.sites.chr21.vcf.bgz chr21:$START-$(( $START + $WS )) | grep -v "#\|insertion\|deletion" | cut -f 2 | sort | wc -l )
                echo $START, $SNPS, $PHASTCONS >> OutputCHR21SNPs.csv
                else
                :
        fi
        START=$(( $START + $STEP ))
done
echo "DONE!"

