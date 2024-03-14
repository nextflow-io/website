<picture>
  <source media="(prefers-color-scheme: dark)" srcset="public/img/nextflow_darkbg.svg">
  <source media="(prefers-color-scheme: light)" srcset="public/img/nextflow.svg">
  <img alt="Nextflow Logo" src="public/img/nextflow.svg">
</picture>

# Nextflow web site

Source code for [https://nextflow.io](https://nextflow.io). Website is generated using [Astro](https://astro.build).

## Commands

To build the website you need Node.js and `npm` installed. See the [npm docs](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm) for instructions.

<!-- TODO: We can replace the Makefile with `npm run` scripts in package.json -->

There is a `Makefile` in the repo with a few standard commands for working with the website:

| Command           | Action                               |
| :---------------- | :----------------------------------- |
| `make dev`        | Run a local server for development   |
| `make build`      | Build the static site for deployment |
| `make publish`    | Publish the website to s3            |
| `make invalidate` | Invalidate the Cloudfront cache      |

The `dev` server builds pages on request and live-updates as you make changes. It includes additional debugging tools. Note that some type-checking on happens during `make build`.

There are also a bunch of `npm` commands for more fine-grained control:

<details>

| Command                   | Action                                           |
| :------------------------ | :----------------------------------------------- |
| `npm install`             | Installs dependencies                            |
| `npm run dev`             | Starts local dev server at `localhost:4321`      |
| `npm run build`           | Build your production site to `./output/`        |
| `npm run preview`         | Preview your build locally, before deploying     |
| `npm run astro ...`       | Run CLI commands like `astro add`, `astro check` |
| `npm run astro -- --help` | Get help using the Astro CLI                     |

</details>

Note that `make build` / `npm run build` fetch the Nextflow repo and build the docs automatically.
When running `make dev` / `npm run dev` this doesn't happen, so the docs will not be visible.

## Code formatting

Code syntax is enforced using [Prettier](https://prettier.io/). You can use a plugin in your code editor to run this automatically on save and it also runs as a CI check on push and pull-requests.

To make it easier to catch formatting errors before pushing, the repo is set up to use [pre-commit](https://pre-commit.com/). If installed, this will run Prettier on any edited files before allowing a commit.

To use, set up as follows:

```bash
pip install pre-commit  # Install the pre-commit tool itself
pre-commit install  # Set up in the Nextflow website repo
```

After that, checks will automatically run every time you do `git commit` in the repo. Note that formatting changes will be applied automatically but the commit will be aborted - you must do `git add .` and `git commit` again.

## Publishing

The Nextflow website is hosted on AWS s3, behind a Cloudfront cache. There are three components to the website:

1. Front end (this repository)
2. Docs (`/docs`, built from the main Nextflow repo using Sphinx)
3. Nextflow binary downloads (`https://get.nextflow.io`, built from the main Nextflow repo)

All three parts of the website are published separately. Care must be taken not to clobber anything by overwriting.

## Project Structure

The repository has the following rough structure:

```text
/
├── public/
└── src/
    ├── content/
    │   ├── blog/
    │   └── podcast/
    ├── components/
    ├── layouts/
    └── pages/
```

Blog posts and podcast episodes go into `content` as markdown files. Other top-level pages should be named as `.astro` or `.md` files in the `src/pages/` directory. Astro finds all of this content and builds a page. Each page is exposed as a route based on its file name.

Any static assets, like images, can be placed in the `public/` directory.

The `src/components/` and `src/layouts/` folders contain the templating files to build the site.

## Documentation

For more info, check the [Astro docs](https://docs.astro.build)
and the old [v3 Bootstrap docs](https://getbootstrap.com/docs/3.4/).
