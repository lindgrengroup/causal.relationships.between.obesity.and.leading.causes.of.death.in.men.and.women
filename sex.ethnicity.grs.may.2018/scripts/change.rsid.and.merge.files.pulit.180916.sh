#!/bin/bash
#$-cwd

################################################################
######### Script for updating rsIDs and merge files #############
################# Jenny Censin, 2018-09-16 ####################

cp -r ../temporary.bed.pulit.180911 ../SAVE.EXTRA

for chr in {1..22}
do
awk 'BEGIN {FS = "[[:space:]]+"} \
{a[$5]; a[$6]; asorti(a,b)
$2=$2"_"b[1]"_"b[2]; delete a}1' \
../temporary.bed.pulit.180911/grs.chr${chr}.bim > \
../temporary.bed.pulit.180911/first${chr}

rm ../temporary.bed.pulit.180911/grs.chr${chr}.bim
mv ../temporary.bed.pulit.180911/first${chr} ../temporary.bed.pulit.180911/grs.chr${chr}.bim
done

#Merge the files
$plink --merge-list ../temporary.bed.pulit.180911/files_to_merge.txt \
--threads 3 \
--memory 45000 \
--make-bed \
--out ../pulit.180916

#Not really needed since that SNP is sorted out
#Change the rsID of "9:140079779_T_C_C_T" to "rs144926207_C_T" to match
#the GRS files
awk 'BEGIN {FS = "[ \t]+"; OFS = "\t"} \
{if ($2 == "9:140079779_T_C_C_T") $2 = "rs144926207_C_T"; print $0 }' \
../pulit.180916.bim > ../pulit.180916.TEMPORARY.bim
cat ../pulit.180916.TEMPORARY.bim > ../pulit.180916.bim
