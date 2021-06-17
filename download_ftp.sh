#!/bin/bash
#Usage: bash download_ftp.sh [link] [dir]

echo "wget -nd -np --recursive -e robots=off --reject "index.html" --no-host-directories --cut-dirs=6 $1 -P $2"
wget -nd -np --recursive -e robots=off --reject "index.html" --no-host-directories --cut-dirs=6 $1 -P $2
cd $2
for i in *.gz
do
  gunzip $i
done
cd ..
ls $2
