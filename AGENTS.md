# AGENTS.md

## Cursor Cloud specific instructions

This repo is the **nextflow.io** website, a static site built with [Astro](https://astro.build). Node 20+ is used (CI uses 20, the devcontainer uses 22). Dependencies are installed by the startup update script (`npm install`), so you normally do not need to install anything manually.

### Running the dev server

- Start it with `npm run dev` (alias for `astro dev`). It serves at `http://localhost:4321/` and live-reloads on file changes.
- Pages live in `src/pages/` (`.astro`/`.md`) and content collections in `src/content/` (blog, podcast). Each file maps to a route by filename.

### Non-obvious gotchas

- **No local blog/docs.** The dev server does NOT build the Nextflow documentation. The site's "Blog" navigation links point to the external `seqera.io` blog, so there is no local `/blog/` route (it 404s) — this is expected. Test rendering with real local pages instead, e.g. `/`, `/our_ambassadors`, `/examples`, `/about-us`.
- **Full build is heavy.** `npm run build` runs `astro check && astro build && bash build_docs.sh`. `build_docs.sh` clones the entire `nextflow-io/nextflow` repo and builds the Sphinx docs (requires network, Python, `jq`, and `pip` deps). For quick local verification of the website itself, run `npx astro build` (skips docs) and/or `npx astro check` (type check) instead. Build output goes to `output/` (gitignored).

### Lint / typecheck

- Formatting is enforced by Prettier via pre-commit (see `.pre-commit-config.yaml`). Check with `npx prettier --check "**/*.{astro,mdx,md,yml,yaml,ts,js,mjs}"`.
- Type checking: `npx astro check`.
