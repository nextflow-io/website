import { defineConfig } from "astro/config";

// https://astro.build/config
export default defineConfig({
  site: "https://nextflow.io/",
  markdown: {
    shikiConfig: {
      // Choose from Shiki's built-in themes (or add your own)
      // https://github.com/shikijs/shiki/blob/main/docs/themes.md
      theme: "min-light",
    },
  },
});
