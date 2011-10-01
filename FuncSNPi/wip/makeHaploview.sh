#!/bin/bash 
COUNT=1
GENOME=
WD=
SIZE=
BIOBED=
REGIONS=
SNP=
while getopts hc:g:r:b: opt
do
    case "$opt" in
        h)
            echo >&2 \
                "COMMAND LINE OPTIONS FOR FuncSNPi:
                -c : max cpus availible (ex. -c 8)
                -g : genome build to use (ex. -g hg19)
                -r : location of risk SNP file (ex. -r ~/riskSNP.txt)
                -b : location of of biological peaks file, in bed format (ex. -b ~/biopeaks.bed)

                riskSNP file must be in the following tab deliminated format
                <region>                <risk SNP id>   <ethnic group>
                10:122837334-123837335  rs2981579       EUR             <- example"
            exit 1;;
        c)
            echo >&2 \
                "$OPTARG cpus are in use"
            COUNT=$OPTARG;;
        g)
            echo >&2 \
                "genome $OPTARG is in use"
            GENOME=$OPTARG;;
        r)
            echo >&2 \
                "$(readlink -f $OPTARG) is the location of the risk SNP file"
            SNP=$(readlink -f $OPTARG);;
        b)
           echo >&2 \
               "$(readlink -f $OPTARG) is the location of the biological peaks file"
            BIOBED=$(readlink -f $OPTARG);;
        :)
            echo >&2 \
                "Option -$OPTARG requires an argument"
            exit 1;;
        \?)
            echo >&2 \
                "COMMAND LINE OPTIONS FOR FuncSNPi:
                -c : max cpus availible (ex. -c 8)
                -g : genome build to use (ex. -g hg19)
                -r : location of risk SNP file (ex. -r ~/riskSNP.txt)
                -b : location of of biological peaks file, in bed format (ex. -b ~/biopeaks.bed)

                riskSNP file must be in the following tab deliminated format
                <region>                <risk SNP id>   <ethnic group>
                10:122837334-123837335  rs2981579       EUR             <- example"
            exit 1;;
    esac
done
shift `expr $OPTIND - 1`

WD=`dirname $SNP`

if [ ${BIOBED: -4} == ".bed" ]; then
    BIOBED=$(dirname $BIOBED)
fi

SIZE=$( echo "-1 * ( $(head -1 $SNP | awk '{print $1}' | cut -d: -f2) ) - 1" | bc )
OUTPUT=$WD/results/$GENOME/$SIZE/
cd $WD

mkdir -p $OUTPUT


makeHaploview () {
for input in ALL ASN AFR EUR; do
    echo "makeHaploview is running"
    echo "current options are region: $region riskname: $riskname ethnic_group: $ethnic chromosome: $chrom"
    echo "INPUT: $input"
    tabix -h /media/bigboy/shared_data/public/SNP/1000_genomes_${input}_20100804.genotypes.vcf.gz $region > ${input}.$riskname.$ethnic.$region.vcf
	vcftools --vcf ${input}.$riskname.$ethnic.$region.vcf --plink --out plinkformat
    awk '{print "chr" $1 "\t" $4 "\t" $4+1 "\t" $2}' plinkformat.map > plinkformat.${input}.$riskname.$ethnic.$region.bed
    for file in `dir $BIOBED/*.bed`; do
        DIR=$(basename ${file%.bed})
        mkdir -p $DIR
        brokenbed=$(awk '{ if (NF != 4) print $1 "\t" $2 "\t" $3 "\t" "error"}' $file | grep error)
        isbed=$(awk '{ if (NF != 4) print "no" }' $file | wc -l)
        if [ $isbed -ne "0" ] ; then
            echo >&2 \
            "bedfile $file, does not have exactly four columns, lines with error are listed in bed_error.log"
            echo $brokenbed > bed_error.log
            exit 1
        else
            echo "the bedfile $file is working"
        fi
        python /home/galaxy/galaxy-dist/tools/new_operations/gops_intersect.py plinkformat.${input}.$riskname.$ethnic.$region.bed $file $DIR/intersect.${input}.$riskname.$ethnic.$region.$(basename $file)
		cut -f1,2 $DIR/intersect.${input}.$riskname.$ethnic.$region.$(basename $file) > $DIR/intersect.${input}.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		grep $riskname ${input}.$riskname.$ethnic.$region.vcf | cut -f1,2 >> $DIR/intersect.${input}.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		sed -i s/chr// $DIR/intersect.${input}.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		sed -i '1ichromosome \t position' $DIR/intersect.${input}.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		vcftools --vcf ${input}.$riskname.$ethnic.$region.vcf --plink --out plinkformat --positions $DIR/intersect.${input}.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		cut -f2,4 plinkformat.map > $DIR/plinkformat.${input}.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers
		mv plinkformat.ped $DIR/plinkformat.${input}.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped
		java -jar /usr/local/bin/Haploview.jar -nogui -log haploview_results/${input}.$riskname.$ethnic.$region.$(basename ${file%.bed})_hap.log -out haploview_results/${input}.$riskname.$ethnic.$region.$(basename ${file%.bed}) -pedfile $DIR/plinkformat.${input}.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped -info $DIR/plinkformat.${input}.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers -skipcheck -dprime -memory 1500
        
        ##
        ##
        #R script can go here. Read in LD files for biobed file N, subset by r.2 value and then plot only those that passed.
        ##
        #need for or while loop to go through the set that passed the r.2 cut off
        #java -jar /usr/local/bin/Haploview.jar -nogui -out haploview_plots/ -pedfile
        
        echo " ${input}.$riskname.$ethnic.$region complete!"
    done
    echo "still here"
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

while read region riskname ethnic chrom
do
    echo $chrom.$riskname.$ethnic.$region
    cd $OUTPUT
    mkdir -p $chrom.$riskname.$ethnic.$region/
    cd $chrom.$riskname.$ethnic.$region/
    mkdir -p haploview_results
    if [ "$COUNT" -gt "1" ]
    then
        parMakeHaploview
    else
    makeHaploview
    echo $chrom.$riskname.$ethnic.$region
    fi
   # ls | grep -E -v '.bed|.pos|.markers|.ped|.vcf' | xargs rm
    ls | grep -E -v '.bed|.pos|.markers|.ped' | xargs rm
done < "$SNP"

echo -e "Finished!"
