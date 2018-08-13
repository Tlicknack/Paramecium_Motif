#!/bin/bash
#Iterate through all Paramecium .gff's to create new table with just genes

parameciumgff=*.gff

for species in $parameciumgff
do
bn=`basename $species -full.gff`
tableFile=$bn-gene.tab
echo "$species -> $bn -> $tableFile"
grep -P "\tgene\t" $species | cut -f 1,4,5,7,9 | cut -d ';' -f 1 > $tableFile
sed -i 's/ID=//' $tableFile
#       cat $species | grep "   gene    " | cut -f 1,4,5,6,9 > {$species}_table
