#!/bin/bash

## IN ORDER TO RUN, NEED 5 ARGUMENTS (hg18 or hg19, path/to/directory/with/list/, name of the file with the list, regionID (e.g. 5kb, 10kb, etc.), path/to/directory/with/biologicalBED_files/
COUNT=1

Genome=$1
if [ -z "$1" ]; then
    echo "specify a genome, hg18 or hg19"
    exit
fi

WD=$2
if [ -z "$2" ]; then
    echo "Need the path to the list of Risk SNP file (e.g. /home/houtana/path/)"
    exit
fi
OUTPUT=$2results/
OUTPUT2=$OUTPUT$Genome/

cd $WD
SNPlist=$3
if [ -z "$3" ]; then
    echo "Need the file name with the list of Risk SNP"
    exit
fi

OUTPUT3=$OUTPUT2$4/

mkdir -p $OUTPUT
mkdir -p $OUTPUT2
mkdir -p $OUTPUT3
echo -e "Using $Genome as the reference\n"

biobed=$5

makeHaploview () {
    echo "makeHaploview is running"
    echo "current options are region: $region riskname: $riskname ethnic_group: $ethnic chromosome: $chrom"
    tabix -h /media/bigboy/shared_data/public/SNP/1000_genomes_${ethnic}_20100804.genotypes.vcf.gz $region > ALL.$riskname.$ethnic.$region.vcf
	vcftools --vcf ALL.$riskname.$ethnic.$region.vcf --plink --out plinkformat
    awk '{print "chr" $1 "\t" $4 "\t" $4+1 "\t" $2}' plinkformat.map > plinkformat.ALL.$riskname.$ethnic.$region.bed
    for file in `dir $biobed/*.bed`; do
        echo "the for loop was entered, we are working on $file"
		DIR=$(basename ${file%.bed})
        mkdir -p $DIR
        isbed=$(awk '{ if (NF != 4) print "no" }' $file | wc -l)
        if [ $isbed -ne "0" ] ; then
            echo >&2 \
            "bedfile $file, does not have exactly four columns"
            exit 1
        else
            echo "the bedfile $file is working"
        fi
        python /home/galaxy/galaxy-dist/tools/new_operations/gops_intersect.py plinkformat.ALL.$riskname.$ethnic.$region.bed $file $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename $file)
		cut -f1,2 $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename $file) > $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		grep $riskname ALL.$riskname.$ethnic.$region.vcf | cut -f1,2 >> $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		sed s/chr// $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos > $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).poss
		sed -e '1ichromosome \t position' $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).poss > $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		vcftools --vcf ALL.$riskname.$ethnic.$region.vcf --plink --out plinkformat --positions $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		cut -f2,4 plinkformat.map > $DIR/plinkformat.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers
		mv plinkformat.ped $DIR/plinkformat.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped
		java -jar /usr/local/bin/Haploview.jar -nogui -log haploview_results/ALL.$riskname.$ethnic.$region.$(basename ${file%.bed})_hap.log -out haploview_results/ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}) -pedfile $DIR/plinkformat.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped -info $DIR/plinkformat.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers -skipcheck -dprime -memory 1500
    done
}

parMakeHaploview () {
    makeHaploview &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge "$COUNT" ]; then
        wait
        NPROC=0
    fi
}

while read region riskname ethnic chrom; do
    cd $OUTPUT3
    mkdir -p $chrom.$riskname.$ethnic.$region/
    cd $chrom.$riskname.$ethnic.$region/
    mkdir -p haploview_results
    if [ "$COUNT" -gt "1" ]
    then
        parMakeHaploview
    else
        makeHaploview
    fi
	ls | grep -E -v '.bed|.poss|.markers|.ped' | xargs rm
done < "$WD/$SNPlist"


echo -e "Finished!, check $OUTPUT2 \n"
