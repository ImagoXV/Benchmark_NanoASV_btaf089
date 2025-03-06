#!/bin/bash

for file in ./*ASV_*; do
sed -i 's/^[[:space:]]*//' $file
done
