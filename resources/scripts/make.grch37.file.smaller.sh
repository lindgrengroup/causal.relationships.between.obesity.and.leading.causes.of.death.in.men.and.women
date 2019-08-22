$bin/bash
$!-cwd

#The pathway is to the file as detailed in the Readme file
echo "chr" "pos" "snp" > ../grch37.all.rsid.txt
zcat ../All_20150605.vcf.gz | awk '{ if (NR > 56) print $1 " " $2 " " $3 }' >> ../grch37.all.rsid.txt
