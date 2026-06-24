# DNS for AI Discovery (DNS-AID)

Publish the records in `nextflow.io.zone` under the `_agents` namespace for `nextflow.io`.

## Requirements

1. Add `_index._agents.nextflow.io` and `_a2a._agents.nextflow.io` SVCB records (see zone file).
2. Enable DNSSEC on the public zone and publish DS records at your domain registrar.
3. Verify with DNS-over-HTTPS:

```bash
curl -s 'https://cloudflare-dns.com/dns-query?name=_index._agents.nextflow.io&type=SVCB' \
  -H 'accept: application/dns-json'
```

## Production note

DNS-AID records are managed outside this repository (Route 53, Cloudflare DNS, or your registrar). Apply the zone file through your DNS operations workflow after review.
