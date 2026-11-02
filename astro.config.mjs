import admonitionsPlugin from "./bin/remark-admonitions.js";
import { defineConfig } from "astro/config";
import remarkDescription from "astro-remark-description";
import remarkDirective from "remark-directive";
import sitemap from "@astrojs/sitemap";
import { transformerNotationDiff, transformerNotationFocus, transformerMetaHighlight } from "@shikijs/transformers";
import react from "@astrojs/react";
import tailwind from "@astrojs/tailwind";

// https://astro.build/config
export default defineConfig({
  site: "https://nextflow.io/",
  outDir: "./output",
  markdown: {
    shikiConfig: {
      // Choose from Shiki's built-in themes (or add your own)
      // https://shiki.style/themes
      theme: "light-plus",
      transformers: [transformerNotationDiff(), transformerNotationFocus(), transformerMetaHighlight()],
    },
    remarkPlugins: [
      remarkDirective,
      admonitionsPlugin,
      [
        remarkDescription,
        {
          name: "excerpt",
          node: (node, i, parent) => {
            // Directly find and return true for the first paragraph that doesn't start with an image
            const firstParagraphIndex = parent?.children.findIndex(
              (child) => child.type === "paragraph" && !child.children[0]?.value?.startsWith("<img "),
            );

            if (i === firstParagraphIndex) {
              return true;
            }

            // If the first paragraph is not what we're looking for, then check for 'end of excerpt' comment
            if (
              parent?.children?.some(
                (child) => child.type === "html" && child.value === "<!-- end-archive-description -->",
              )
            ) {
              const sibling = parent?.children[i + 1];
              return sibling?.type === "html" && sibling?.value === "<!-- end-archive-description -->";
            }

            return false;
          },
        },
      ],
    ],
  },
  integrations: [
    react(),
    tailwind(),
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
