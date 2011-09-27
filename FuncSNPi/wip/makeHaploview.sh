#!/bin/bash

## IN ORDER TO RUN, NEED 5 ARGUMENTS (hg18 or hg19, path/to/directory/with/list/, name of the file with the list, regionID (e.g. 5kb, 10kb, etc.), path/to/directory/with/biologicalBED_files/

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

#cd $OUTPUT3

#for i in `cat $WD/$SNPlist`; do
while read region riskname ethnic chrom; do
cd $OUTPUT3
#	region="$($i | cut -f1)"
#	riskname="$($i | cut -f2)"
#	ethnic="$($i | cut -f3)"
#name = $chrom.$riskname.$ethnic.$region
mkdir -p $chrom.$riskname.$ethnic.$region/
cd $chrom.$riskname.$ethnic.$region/
mkdir -p haploview_results
##ALL
	tabix -h /media/bigboy/shared_data/public/SNP/1000_genomes_ALL_20100804.genotypes.vcf.gz $region > ALL.$riskname.$ethnic.$region.vcf
	vcftools --vcf ALL.$riskname.$ethnic.$region.vcf --plink --out plinkformat
	awk '{print "chr" $1 "\t" $4 "\t" $4+1 "\t" $2}' plinkformat.map > plinkformat.ALL.$riskname.$ethnic.$region.bed
	for file in `dir $5/*.bed`; do
		DIR=$(basename ${file%.bed})
		mkdir -p $DIR
		python /home/galaxy/galaxy-dist/tools/new_operations/gops_intersect.py plinkformat.ALL.$riskname.$ethnic.$region.bed $file $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename $file)
		cut -f1,2 $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename $file) > $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		#echo $riskname >> $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).snps
		grep $riskname ALL.$riskname.$ethnic.$region.vcf | cut -f1,2 >> $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		sed s/chr// $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos > $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).poss
		sed -e '1ichromosome \t position' $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).poss > $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		vcftools --vcf ALL.$riskname.$ethnic.$region.vcf --plink --out plinkformat --positions $DIR/intersect.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		cut -f2,4 plinkformat.map > $DIR/plinkformat.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers
		mv plinkformat.ped $DIR/plinkformat.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped
		java -jar /usr/local/bin/Haploview.jar -nogui -log haploview_results/ALL.$riskname.$ethnic.$region.$(basename ${file%.bed})_hap.log -out haploview_results/ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}) -pedfile $DIR/plinkformat.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped -info $DIR/plinkformat.ALL.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers -svg -ldcolorscheme RSQ -ldvalues RSQ -skipcheck -dprime -memory 1500
	done
##AFR
	tabix -h /media/bigboy/shared_data/public/SNP/1000_genomes_AFR_20100804.genotypes.vcf.gz $region > AFR.$riskname.$ethnic.$region.vcf
	vcftools --vcf AFR.$riskname.$ethnic.$region.vcf --plink --out plinkformat
	awk '{print "chr" $1 "\t" $4 "\t" $4+1 "\t" $2}' plinkformat.map > plinkformat.AFR.$riskname.$ethnic.$region.bed
	for file in `dir $5/*.bed`; do
		DIR=$(basename ${file%.bed})
		mkdir -p $DIR
		python /home/galaxy/galaxy-dist/tools/new_operations/gops_intersect.py plinkformat.AFR.$riskname.$ethnic.$region.bed $file $DIR/intersect.AFR.$riskname.$ethnic.$region.$(basename $file)
		cut -f1,2 $DIR/intersect.AFR.$riskname.$ethnic.$region.$(basename $file) > $DIR/intersect.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		#echo $riskname >> $DIR/intersect.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).snps
		grep $riskname AFR.$riskname.$ethnic.$region.vcf | cut -f1,2 >> $DIR/intersect.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		sed s/chr// $DIR/intersect.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos > $DIR/intersect.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).poss
		sed -e '1ichromosome \t position' $DIR/intersect.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).poss > $DIR/intersect.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		vcftools --vcf AFR.$riskname.$ethnic.$region.vcf --plink --out plinkformat --positions $DIR/intersect.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		cut -f2,4 plinkformat.map > $DIR/plinkformat.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers
		mv plinkformat.ped $DIR/plinkformat.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped
		java -jar /usr/local/bin/Haploview.jar -nogui -log haploview_results/AFR.$riskname.$ethnic.$region.$(basename ${file%.bed})_hap.log -out haploview_results/AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}) -pedfile $DIR/plinkformat.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped -info $DIR/plinkformat.AFR.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers -svg -ldcolorscheme RSQ -ldvalues RSQ -skipcheck -dprime -memory 1500
	done
##EUR
	tabix -h /media/bigboy/shared_data/public/SNP/1000_genomes_EUR_20100804.genotypes.vcf.gz $region > EUR.$riskname.$ethnic.$region.vcf
	vcftools --vcf EUR.$riskname.$ethnic.$region.vcf --plink --out plinkformat
	awk '{print "chr" $1 "\t" $4 "\t" $4+1 "\t" $2}' plinkformat.map > plinkformat.EUR.$riskname.$ethnic.$region.bed
	for file in `dir $5/*.bed`; do
		DIR=$(basename ${file%.bed})
		mkdir -p $DIR
		python /home/galaxy/galaxy-dist/tools/new_operations/gops_intersect.py plinkformat.EUR.$riskname.$ethnic.$region.bed $file $DIR/intersect.EUR.$riskname.$ethnic.$region.$(basename $file)
		cut -f1,2 $DIR/intersect.EUR.$riskname.$ethnic.$region.$(basename $file) > $DIR/intersect.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		#echo $riskname >> $DIR/intersect.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).snps
		grep $riskname EUR.$riskname.$ethnic.$region.vcf | cut -f1,2 >> $DIR/intersect.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		sed s/chr// $DIR/intersect.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos > $DIR/intersect.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).poss
		sed -e '1ichromosome \t position' $DIR/intersect.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).poss > $DIR/intersect.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		vcftools --vcf EUR.$riskname.$ethnic.$region.vcf --plink --out plinkformat --positions $DIR/intersect.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		cut -f2,4 plinkformat.map > $DIR/plinkformat.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers
		mv plinkformat.ped $DIR/plinkformat.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped
		java -jar /usr/local/bin/Haploview.jar -nogui -log haploview_results/EUR.$riskname.$ethnic.$region.$(basename ${file%.bed})_hap.log -out haploview_results/EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}) -pedfile $DIR/plinkformat.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped -info $DIR/plinkformat.EUR.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers -svg -ldcolorscheme RSQ -ldvalues RSQ -skipcheck -dprime -memory 1500
	done
##ASN
	tabix -h /media/bigboy/shared_data/public/SNP/1000_genomes_ASN_20100804.genotypes.vcf.gz $region > ASN.$riskname.$ethnic.$region.vcf
	vcftools --vcf ASN.$riskname.$ethnic.$region.vcf --plink --out plinkformat
	awk '{print "chr" $1 "\t" $4 "\t" $4+1 "\t" $2}' plinkformat.map > plinkformat.ASN.$riskname.$ethnic.$region.bed
	for file in `dir $5/*.bed`; do
		DIR=$(basename ${file%.bed})
		mkdir -p $DIR
		python /home/galaxy/galaxy-dist/tools/new_operations/gops_intersect.py plinkformat.ASN.$riskname.$ethnic.$region.bed $file $DIR/intersect.ASN.$riskname.$ethnic.$region.$(basename $file)
		cut -f1,2 $DIR/intersect.ASN.$riskname.$ethnic.$region.$(basename $file) > $DIR/intersect.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		#echo $riskname >> $DIR/intersect.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).snps
		grep $riskname ASN.$riskname.$ethnic.$region.vcf | cut -f1,2 >> $DIR/intersect.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		sed s/chr// $DIR/intersect.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos > $DIR/intersect.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).poss
		sed -e '1ichromosome \t position' $DIR/intersect.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).poss > $DIR/intersect.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		vcftools --vcf ASN.$riskname.$ethnic.$region.vcf --plink --out plinkformat --positions $DIR/intersect.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).pos
		cut -f2,4 plinkformat.map > $DIR/plinkformat.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers
		mv plinkformat.ped $DIR/plinkformat.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped
		java -jar /usr/local/bin/Haploview.jar -nogui -log haploview_results/ASN.$riskname.$ethnic.$region.$(basename ${file%.bed})_hap.log -out haploview_results/ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}) -pedfile $DIR/plinkformat.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).ped -info $DIR/plinkformat.ASN.$riskname.$ethnic.$region.$(basename ${file%.bed}).markers -svg -ldcolorscheme RSQ -ldvalues RSQ -skipcheck -dprime -memory 1500
	done
	ls | grep -E -v '.bed|.poss|.markers|.ped' | xargs rm



done <"$WD/$SNPlist"


echo -e "Finished!, check $OUTPUT2 \n"
