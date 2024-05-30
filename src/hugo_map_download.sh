#!/bin/bash

# check if the data directory exists, and create it if not
if [ ! -d ./datas ]; then
    echo "Creating ./data"
    mkdir ./data
fi

# change into networks directory
cd ./data

wget https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt

