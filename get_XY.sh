#!/bin/bash

x=`samtools idxstats $1 | egrep  "X|chrX"`;
y=`samtools idxstats $1 | egrep  "Y|chrY"`;
xcov=`echo $x | awk '{ printf("%0.5f", $3/$2); }'`;
ycov=`echo $y | awk '{ printf("%0.5f", $3/$2); }'`;

rat=$(echo "scale=4; ${xcov}/${ycov}" | bc)
if (( $(echo "$rat > 5.0" | bc -l) )); then
    sex=female
else
    sex=male
fi
echo "$sex"