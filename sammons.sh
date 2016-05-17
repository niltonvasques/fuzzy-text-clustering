#!/bin/bash
for i in $(ls -d results/*/); do 
  METHOD=$(echo $i | cut -f1 -d'-' | cut -f2 -d'/')
  echo $METHOD; BASE=$(echo $i | cut -f2 -d'-' | cut -f1 -d'/')
  echo $BASE
  Rscript sammons.r $BASE.frequencys $METHOD $i &
done

# Wait for all parallel jobs to finish
while [ 1 ]; do fg 2> /dev/null; [ $? == 1 ] && break; done
