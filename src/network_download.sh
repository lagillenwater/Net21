#!/bin/bash

# check if the networks directory exists, and create it if not
if [ ! -d ./networks ]; then
    echo "Creating ./networks"
    mkdir ./networks
fi

# change into networks directory
cd ./networks

## download the full blood network from HumanBase
echo "Downloading the full blood network from HumanBase"
wget https://s3-us-west-2.amazonaws.com/humanbase/networks/blood_top.gz
echo "Download Complete"

### need to unzip network
echo "Unzipping network"
gunzip blood_top.gz
echo "Unzipping complete"
