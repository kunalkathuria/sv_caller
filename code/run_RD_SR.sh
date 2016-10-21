cp ../results/text/ClassifiedVariantMap.txt ../results/text/VariantMap.txt
python add_SR.py $1 $2 $3 $4 # last 2 are "slop" and "refresh SR margin" used in code
cp ../results/text/VariantMap_SR.txt ../results/text/VariantMap.txt
cp ../results/text/All_Variants_SR.txt ../results/text/All_Variants.txt
cat ../results/text/All_Variants.txt | grep -v INV | grep -v known > ../results/text/inDels.txt # all INDELS
sort -k6,6 ../results/text/inDels.txt > ../results/text/inDels_S.txt
python add_RD.py .70 ../results/text/RDSegments.txt
cp ../results/text/VariantMap_RD.txt ../results/text/VariantMap.txt
cp ../results/text/All_Variants_RD.txt ../results/text/All_Variants.txt
python SetCover.py $5
python WriteBed.py
cp ../results/text/VariantMap_O.txt ../results/text/VariantMap.txt
cp ../results/text/All_Variants_O.txt ../results/text/All_Variants.txt
