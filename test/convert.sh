#!/bin/bash
nc=`grep -v "#" $1 |head -n 1 |awk '{print NF}' `
nl=`grep -v "#" $1 |wc -l |awk '{print $1}' `
./file2bin.x <(echo $nl" "$nc | cat - <(grep -v "#" $1)) ${1}_new.bin 
