#!/bin/bash

#####################################################################
# NEXTFLOW DOCUMENTATION
#####################################################################

WEBSITE_DIR=$(pwd)

# Make the empty target directories
mkdir -p output/docs/latest/
mkdir -p output/docs/edge/
mkdir -p output/docs/master/

# Fetch the Nextflow repo
git clone https://github.com/nextflow-io/nextflow.git
cd nextflow/docs/

# Find the latest stable and edge releases
STABLE_TAG=$(curl -s https://api.github.com/repos/nextflow-io/nextflow/releases | jq -r ". [].tag_name" | grep -v edge | head -n 1)
EDGE_TAG=$(curl -s https://api.github.com/repos/nextflow-io/nextflow/releases | jq -r ". [].tag_name" | grep edge | head -n 1)
echo "Latest stable release: $STABLE_TAG"
echo "Latest edge release: $EDGE_TAG"

# Sanity check: Assert string lengths
if [ ${#STABLE_TAG} -le 4 ]; then echo "Version string too short" ; exit 1; fi
if [ ${#EDGE_TAG} -le 4 ]; then echo "Version string too short" ; exit 1; fi

# Build edge docs
git checkout $EDGE_TAG
pip install -r requirements.txt
make clean html
mv _build/html/* $WEBSITE_DIR/output/docs/edge/

# Build stable docs
git checkout $STABLE_TAG
pip install -r requirements.txt
make clean html
mv _build/html/* $WEBSITE_DIR/output/docs/latest/

# Build current docs on master
git checkout master
pip install -r requirements.txt
make clean html
mv _build/html/* $WEBSITE_DIR/output/docs/master/
