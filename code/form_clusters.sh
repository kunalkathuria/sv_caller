set -eux

${12} sort $2 ${2}.sorted
mv $2.sorted.bam $2
${12} index $2
time python ReadDiscordants.py $1 $2 $3 $4 $5 $6 $7 # 9 minutes on 50 x homozygous

echo Sorting
sort -T ../results/text/ -k 2,2 -k 3,3n ../results/text/All_Discords_P.txt > ../results/text/All_Discords_P_S.txt
echo DoneSorting

#Line below is for random insertions only
#sort -n -k 1,1 ../results/text/All_Discords_I.txt > ../results/text/All_Discords_I_S.txt

time (python FormClusters.py ../results/text/bam_stats.txt $8 $9 ${10} ${11}) # 50m
time (sort -T ../results/text/ -n -k 2,2 ../results/text/VariantMapInp_P.txt > ../results/text/VariantMapInp.txt) #1m

time python WriteClusterMap.py $8 # 1m
cp ../results/text/VariantMap.txt ../results/text/VariantMap_O.txt
