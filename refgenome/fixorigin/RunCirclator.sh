#!/bin/bash

FILES=/home/chrisgaby/tmp/SEntericaChrom/refChrom/*.fa

for n in $FILES 
  do
  circlator fixstart $n /home/chrisgaby/tmp/SEntericaChrom/FixedOrigen/"$(basename ${n%%.fa})_OriginSet"
done
