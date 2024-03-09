import { defineConfig } from "astro/config";

// https://astro.build/config
export default defineConfig({
  site: "https://nextflow.io/",
  markdown: {
    shikiConfig: {
      // Choose from Shiki's built-in themes (or add your own)
      // https://shiki.style/themes
      theme: "light-plus",
    },
  },
});
