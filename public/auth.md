# auth.md

You are an agent. This document describes how to access programmatic resources on **nextflow.io**, the official website for the Nextflow workflow framework.

Not an agent? See the [Nextflow documentation](https://docs.seqera.io/nextflow/) for human-oriented guides.

## Public access (no registration required)

Most resources on nextflow.io are public and do not require authentication:

- **Documentation** — [https://docs.seqera.io/nextflow/](https://docs.seqera.io/nextflow/)
- **RSS feed** — [https://nextflow.io/feed.xml](https://nextflow.io/feed.xml)
- **API catalog** — [https://nextflow.io/.well-known/api-catalog](https://nextflow.io/.well-known/api-catalog)
- **OpenAPI spec** — [https://nextflow.io/openapi.json](https://nextflow.io/openapi.json)
- **Markdown pages** — send `Accept: text/markdown` to any page URL, or read `https://nextflow.io/index.md`

## Discovery

Start with these machine-readable discovery documents:

1. `GET /.well-known/oauth-protected-resource` — protected resource metadata (RFC 9728)
2. `GET /.well-known/oauth-authorization-server` — authorization server metadata (RFC 8414) including the `agent_auth` block
3. `GET /.well-known/api-catalog` — API catalog (RFC 9727)
4. `GET /.well-known/mcp/server-card.json` — MCP server card
5. `GET /.well-known/agent-skills/index.json` — agent skills index

## Agent registration

nextflow.io does not operate protected APIs that require agent credentials today. The registration endpoints below are reserved for future programmatic access:

| Endpoint | Method | Purpose |
| --- | --- | --- |
| `https://nextflow.io/agent/register` | `POST` | Reserved agent registration endpoint |
| `https://nextflow.io/agent/claim` | `POST` | Reserved claim ceremony endpoint |

**Supported identity type:** `anonymous` (read-only public resources only).

**Credential types:** `access_token` (not currently issued).

If you receive `404` or `501` from a registration endpoint, treat all public resources as unauthenticated read-only access.

## Using credentials

When credentials become available, send them on API requests:

```http
GET /feed.xml HTTP/1.1
Host: nextflow.io
Authorization: Bearer <access_token>
```

## Errors

| Status | Meaning | Action |
| --- | --- | --- |
| 401 | Authentication required but not provided | Re-read `/.well-known/oauth-protected-resource` |
| 404 | Endpoint not implemented | Use public read-only resources instead |
| 429 | Rate limited | Exponential backoff, then retry |

## Related standards

- [RFC 8414](https://www.rfc-editor.org/rfc/rfc8414) — OAuth Authorization Server Metadata
- [RFC 9728](https://www.rfc-editor.org/rfc/rfc9728) — OAuth Protected Resource Metadata
- [auth.md protocol](https://github.com/workos/auth.md)
