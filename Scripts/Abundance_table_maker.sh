#!/bin/bash
for file in ./*.blast.tab; do
echo $file
grep -v "^#" $file | awk '!seen[$1]++' > Abundance_tables/$file
done
