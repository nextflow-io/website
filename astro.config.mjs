import { defineConfig } from "astro/config";

import sitemap from "@astrojs/sitemap";

const netlifyPrimeURL = import.meta.env.DEPLOY_PRIME_URL;
const siteURL = netlifyPrimeURL ? netlifyPrimeURL : "https://nextflow.io";

// https://astro.build/config
export default defineConfig({
  site: siteURL + "/",
  outDir: "./output",
  markdown: {
    shikiConfig: {
      // Choose from Shiki's built-in themes (or add your own)
      // https://shiki.style/themes
      theme: "light-plus",
    },
  },
  integrations: [
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
