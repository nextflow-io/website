import { defineConfig } from "astro/config";
import remarkDescription from "astro-remark-description";
import sitemap from "@astrojs/sitemap";

// https://astro.build/config
export default defineConfig({
  site: "https://nextflow.io/",
  outDir: "./output",
  markdown: {
    shikiConfig: {
      // Choose from Shiki's built-in themes (or add your own)
      // https://shiki.style/themes
      theme: "light-plus",
    },
    remarkPlugins: [
      [
        remarkDescription,
        {
          name: "excerpt",
          // transform: (desc) => desc.split(' ').length < 150 ? desc : `${desc.split(' ').slice(150).join(' ')}...`,
          node: (node, i, parent) => {
            // check if parent has a child that is an html comment with the text 'end of excerpt'
            if (
              parent?.children?.some(
                (child) => child.type === "html" && child.value === "<!-- end-archive-description -->",
              )
            ) {
              const sibling = parent?.children[i + 1];
              return sibling?.type === "html" && sibling?.value === "<!-- end-archive-description -->";
            } else {
              // return the first paragraph otherwise

              // get the index of the first paragraph and check that it doesn't start with an image
              const firstParagraphIndex = parent?.children.findIndex((child) => child.type === "paragraph" && !child.children[0].value?.startsWith('<img '));
              
              // if the node is the first paragraph, return true
              return i === firstParagraphIndex;
            }
          },
        },
      ],
    ],
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
