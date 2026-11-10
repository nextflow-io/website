import { defineConfig } from "astro/config";
import sitemap from "@astrojs/sitemap";
import react from "@astrojs/react";
import tailwindcss from "@tailwindcss/vite";
import { createRenderer } from "astro-expressive-code";
import { pluginLineNumbers } from "@expressive-code/plugin-line-numbers";
import codeblocks from "starlight-codeblocks/astro";
import { nextflow } from "starlight-codeblocks/adapters/nextflow";

// Local modules `include`d by the example pipeline pages, linked to their source repos
const exampleModules = {
  RNASEQ: "rnaseq-nf/blob/master/modules/rnaseq.nf",
  Sample: "rnaseq-nf/blob/master/modules/rnaseq.nf",
  AlignedSample: "rnaseq-nf/blob/master/modules/rnaseq.nf",
  MULTIQC: "rnaseq-nf/tree/master/modules/multiqc",
  FETCH_DATASET: "hyperopt/tree/master/modules/fetch_dataset",
  SPLIT_TRAIN_TEST: "hyperopt/tree/master/modules/split_train_test",
  VISUALIZE: "hyperopt/tree/master/modules/visualize",
  TRAIN: "hyperopt/tree/master/modules/train",
  EVALUATE: "hyperopt/tree/master/modules/evaluate",
};

// bootstrap.css sets the root font size to 62.5%, so scale code block rem sizes back up
const scaleRem = (css) => css.replace(/(\d*\.?\d+)rem\b/g, "calc($1rem * 1.6)");

// https://astro.build/config
export default defineConfig({
  site: "https://nextflow.io/",
  outDir: "./output",
  // Astro 7's default ('jsx') drops spaces between inline elements in .astro templates
  compressHTML: true,
  vite: {
    plugins: [tailwindcss()],
  },
  integrations: [
    codeblocks({
      expressiveCode: {
        themes: ["light-plus", "dark-plus"],
        useDarkModeMediaQuery: false,
        themeCssSelector: (theme) => `[data-theme='${theme.type}']`,
        plugins: [pluginLineNumbers()],
        customCreateRenderer: async (options) => {
          const renderer = await createRenderer(options);
          renderer.baseStyles = scaleRem(renderer.baseStyles);
          renderer.themeStyles = scaleRem(renderer.themeStyles);
          return renderer;
        },
        defaultProps: {
          wrap: true,
          overridesByLang: { "bash,sh,shell,console": { showLineNumbers: false } },
        },
      },
      // The bundled vscode-icons rules only match a file named exactly `.nf` (ewels/starlight-codeblocks#29)
      fileIcons: { files: { ".nf": "file-type-nextflow" } },
      apiLinks: {
        adapters: [
          nextflow({
            modules: ({ name }) => {
              const path = exampleModules[name];
              if (!path) return undefined;
              const [repo] = path.split("/");
              return { href: `https://github.com/nextflow-io/${path}`, source: `nextflow-io/${repo}`, icon: "github" };
            },
          }),
        ],
      },
    }),
    react(),
    sitemap({
      // Just copying what was manually added in the old jbake site for now
      customPages: [
        "https://www.nextflow.io/docs/latest/index.html",
        "https://www.nextflow.io/docs/latest/getstarted.html",
        "https://www.nextflow.io/docs/latest/basic.html",
        "https://www.nextflow.io/docs/latest/script.html",
        "https://www.nextflow.io/docs/latest/process.html",
        "https://www.nextflow.io/docs/latest/channel.html",
        "https://www.nextflow.io/docs/latest/operator.html",
        "https://www.nextflow.io/docs/latest/executor.html",
        "https://www.nextflow.io/docs/latest/config.html",
        "https://www.nextflow.io/docs/latest/amazons3.html",
        "https://www.nextflow.io/docs/latest/docker.html",
        "https://www.nextflow.io/docs/latest/dnanexus.html",
        "https://www.nextflow.io/docs/latest/ignite.html",
        "https://www.nextflow.io/docs/latest/tracing.html",
        "https://www.nextflow.io/docs/latest/sharing.html",
        "https://www.nextflow.io/docs/latest/metadata.html",
      ],
    }),
  ],
  build: {
    // Same as old site - eg. /ambassadors.html instead of /ambassadors/
    format: "file",
  },
});
