publish:
	aws s3 sync output/ s3://www2.nextflow.io \
	 --cache-control max-age=2592000 \
	 --metadata-directive REPLACE \
	 --storage-class STANDARD \
	 --acl public-read \
	 --exclude "$0" \
	 "$@"

	aws s3 sync --delete output/docs s3://www2.nextflow.io/docs \
	 --cache-control max-age=2592000 \
	 --metadata-directive REPLACE \
	 --storage-class STANDARD \
	 --acl public-read \
	 --exclude "$0" \
	 "$@"

invalidate:
	aws cloudfront create-invalidation --distribution-id E3RPV5P71OW0UF --paths '/*'

make build:
	./jbake
