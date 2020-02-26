cat "$1" | awk -F "\t" '{print $1"\t"$2-100000"\t"$3+100000}' | sed 's/-[0-9]*/0/g' | bedtools sort | bedtools merge
