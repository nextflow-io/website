# CloudFront agent-readiness configuration

Production hosting uses S3 + CloudFront. Apply these after deploying static files.

## Link response headers (RFC 8288)

1. Create a response headers policy from `response-headers-policy.json`, or attach equivalent `Link` headers to the default cache behavior for `/` and `/index.html`.
2. Alternatively, rely on `scripts/set-s3-agent-metadata.sh` (runs during `make publish`) which sets S3 object metadata on `index.html`.

## Markdown content negotiation

1. Publish the CloudFront Function in `markdown-negotiation.js`.
2. Associate it with **viewer-request** on the default behavior.
3. Ensure `index.md` is deployed to S3 (copied from `public/index.md` via the Astro build).

The function rewrites `/` and `/index.html` to `/index.md` when `Accept: text/markdown` is present.

## Netlify previews

PR previews use `public/_headers` for Link headers and `netlify/edge-functions/markdown-negotiate.js` for markdown negotiation.
