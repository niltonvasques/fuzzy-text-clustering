#!/bin/bash

grep ":real\." discover.names | grep -o '"[a-z]\+"' | sed -s 's/"//g' > terms 
cat discover.data  | awk 'BEGIN {FS=","} {for(x=0; x<=NF;x++){ if($x + 0 == $x) printf("%s ",$x); } print "" }' >  frequencys

cat <(wc -l terms | awk '{print $1}') <(wc -l frequencys | awk '{print $1}') terms frequencys > in
rm frequencys
rm terms
