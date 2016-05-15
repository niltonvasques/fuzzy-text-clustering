#!/bin/bash
FOLDER=$1
NAME=$2

grep ":real\." $FOLDER/discover.names | grep -o '"[a-z]\+"' | sed -s 's/"//g' > terms 
cat $FOLDER/discover.data  | awk 'BEGIN {FS=","} {for(x=0; x<=NF;x++){ if($x + 0 == $x) printf("%s ",$x); } print "" }' >  frequencys

cat <(wc -l terms | awk '{print $1}') <(wc -l frequencys | awk '{print $1}') terms frequencys > $NAME.in
rm frequencys
rm terms
