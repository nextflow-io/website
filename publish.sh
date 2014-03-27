#!/bin/bash
s3cmd sync * s3://www.nextflow.io \
 -rr \
 --acl-public \
 --no-mime-magic \
 --access_key=$NXF_AWS_ACCESS \
 --secret_key=NXF_AWS_SECRET \
 --exclude 'publish.sh' \
 "$@"
