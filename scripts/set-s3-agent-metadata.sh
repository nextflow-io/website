#!/bin/bash
# Set Content-Type metadata on agent-discovery objects after S3 sync.
# Invoked by publish.sh after the main sync completes.

set -euo pipefail

BUCKET="${NXF_S3_BUCKET:-www2.nextflow.io}"

declare -A CONTENT_TYPES=(
  [".well-known/api-catalog"]="application/linkset+json"
  [".well-known/oauth-authorization-server"]="application/json"
  [".well-known/oauth-protected-resource"]="application/json"
  [".well-known/jwks.json"]="application/json"
  [".well-known/mcp/server-card.json"]="application/json"
  [".well-known/agent-skills/index.json"]="application/json"
  [".well-known/agent-skills/nextflow-docs/SKILL.md"]="text/markdown; charset=utf-8"
  ["auth.md"]="text/markdown; charset=utf-8"
  ["index.md"]="text/markdown; charset=utf-8"
  ["openapi.json"]="application/json"
)

for key in "${!CONTENT_TYPES[@]}"; do
  content_type="${CONTENT_TYPES[$key]}"
  if aws s3api head-object --bucket "$BUCKET" --key "$key" >/dev/null 2>&1; then
    aws s3 cp "s3://${BUCKET}/${key}" "s3://${BUCKET}/${key}" \
      --metadata-directive REPLACE \
      --content-type "$content_type" \
      --acl public-read
    echo "Set Content-Type for ${key} -> ${content_type}"
  else
    echo "Skipping ${key} (not found in bucket)"
  fi
done

# Link header on homepage objects
LINK_HEADER='</.well-known/api-catalog>; rel="api-catalog", </docs.seqera.io/nextflow/>; rel="service-doc", </openapi.json>; rel="service-desc"'

for key in index.html; do
  if aws s3api head-object --bucket "$BUCKET" --key "$key" >/dev/null 2>&1; then
    aws s3 cp "s3://${BUCKET}/${key}" "s3://${BUCKET}/${key}" \
      --metadata-directive REPLACE \
      --content-type "text/html; charset=utf-8" \
      --metadata "Link=${LINK_HEADER}" \
      --acl public-read
    echo "Set Link metadata for ${key}"
  fi
done
