#!/bin/bash

if [ $# -ne 3 ];then
    echo "Usage: HMMAlignment.sh <src_file> <trg_file.es> <dest_dir>"
    exit 0
fi

prefix=$3

mkdir -p $prefix
text_1=$1 
text_2=$2

name_1=`basename ${text_1}`
name_2=`basename ${text_2}`

echo "Obtaining target and source vocabulary files (.vcb) and sentence pair file (.snt)."

${GIZA}/plain2snt ${text_1} ${text_2} -snt1  ${prefix}/${name_1}_${name_2}.snt -snt2  ${prefix}/${name_2}_${name_1}.snt -vcb1  ${prefix}/${name_1}.vcb -vcb2  ${prefix}/${name_2}.vcb
echo "Done"
echo "Creating coocurrence file"
${GIZA}/snt2cooc   ${prefix}/${name_1}_${name_2}.cooc  ${prefix}/${name_1}.vcb  ${prefix}/${name_2}.vcb  ${prefix}/${name_1}_${name_2}.snt
${GIZA}/snt2cooc   ${prefix}/${name_2}_${name_1}.cooc  ${prefix}/${name_2}.vcb  ${prefix}/${name_1}.vcb  ${prefix}/${name_2}_${name_1}.snt

echo "Aligning with GIZA++"
echo "One way (${name_1} -> ${name_2})"

${GIZA}/mgiza -s ${prefix}/${name_1}.vcb -t  ${prefix}/${name_2}.vcb -c  ${prefix}/${name_1}_${name_2}.snt -coocurrencefile  ${prefix}/${name_1}_${name_2}.cooc -o ${prefix}/${name_1}_${name_2} -m1 19 -m2 0 -m3 0 -m4 0 -mh 0 -m5 0 -model1dumpfrequency 1 -nodumps 0 -emprobforempty 0.0 -probsmooth 0.0

echo "Getting alignments"

python - ${prefix}/${name_1}.vcb ${prefix}/${name_2}.vcb ${prefix}/*.t1.19 <<EOF > ${prefix}/alignments
import sys
src = {}
for line in open(sys.argv[1]):
 src[int(line.split()[0])] = line.split()[1]
src[0] = 'NULL'
trg = {}
for line in open(sys.argv[2]):
 trg[int(line.split()[0])] = line.split()[1]
trg[0] = 'NULL'
for line in open(sys.argv[3]):
 print src[int(line.split()[0])], trg[int(line.split()[1])], line.split()[2]
EOF

rm ${prefix}/${name_1}* ${prefix}/${name_2}*

echo "Done"
