#for i in {1..21} X; do 
i=21
 perl  Generate_all_variants.pl test/GencodeV19_splice.chr${i}.gz  test/GencodeV19_chr${i}_variants.txt > test/21.log
#done;
