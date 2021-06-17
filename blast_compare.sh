S1=$1
S2=$2
SN1=${S1%%.*}
SN2=${S2%%.*}

makeblastdb -in ${S1} -out db/${SN1} -dbtype nucl
blastn -query ${S2} -db db/${SN1} -outfmt 6 -out ${SN1}_vs_${SN2}.tsv
