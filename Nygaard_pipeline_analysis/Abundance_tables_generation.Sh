#!/bin/bash

for file in ./*.blast.tab; do
cut -f2 $file | sort | uniq -c > ASV_TABLE_$(basename $file)
done
