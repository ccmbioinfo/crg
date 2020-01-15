report=$1
cohort=$2

cat ${report} | awk '{printf $1"\t"$2"\t"$3"\t"$4"\n"}' > ${cohort}".tsv"
