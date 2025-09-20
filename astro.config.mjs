import admonitionsPlugin from "./bin/remark-admonitions.js";
import { defineConfig } from "astro/config";
import remarkDescription from "astro-remark-description";
import remarkDirective from "remark-directive";
import sitemap from "@astrojs/sitemap";
import expressiveCode from "astro-expressive-code";
import react from "@astrojs/react";
import tailwind from "@astrojs/tailwind";

// https://astro.build/config
export default defineConfig({
  site: "https://nextflow.io/",
  outDir: "./output",
  markdown: {
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
    expressiveCode({
      // Use light theme for consistency with site design
      themes: ["github-light"],
      // Enable advanced features
      plugins: ["text-markers", "frames"],
      // Comprehensive styling to match site's clean, minimal design
      styleOverrides: {
        // Core background and text styling - pure white like site
        codeBackground: "#ffffff",
        codeForeground: "#24292f", // GitHub light theme text color

        // Typography - match site's exact monospace stack
        codeFontFamily: "Menlo, Monaco, Consolas, 'Courier New', monospace",
        codeFontSize: "0.875rem", // 14px - readable code size
        codeLineHeight: "1.5",

        // Borders - subtle gray matching site's container styling
        borderColor: "#e5e7eb", // rgb(229, 231, 235) - light gray
        borderWidth: "1px",
        borderRadius: "0.375rem", // Tailwind rounded-md

        // Text markers - use site's green color scheme for highlighting
        textMarkers: {
          // Use the site's actual light green color for highlighting
          markBackground: "var(--nextflow-light-green)", // Direct use of site color
          markBorderColor: "transparent", // Clean, no borders

          // Make highlighting more visible but still clean
          backgroundOpacity: "0.4", // More visible highlighting
          borderOpacity: "0", // No border opacity
        },

        // Frames - clean, minimal styling
        frames: {
          // Remove all shadows and heavy styling
          shadowColor: "transparent",
          frameBoxShadowCssValue: "none",

          // Clean frame styling matching site containers
          editorBackground: "#ffffff",
          editorActiveTabBackground: "#f9fafb", // Very light gray for active tab
          editorActiveTabBorderColor: "#e5e7eb", // Match border color
          editorTabBarBackground: "#ffffff",
          editorTabBarBorderColor: "#e5e7eb",

          // Terminal-style frame styling
          terminalBackground: "#ffffff",
          terminalTitlebarBackground: "#f9fafb",
          terminalTitlebarForeground: "#6b7280", // Subtle gray text
          terminalTitlebarBorderColor: "#e5e7eb",
        },

        // Focus and selection states
        focusBorder: "var(--nextflow-green)", // Use site's green for focus
        codeSelectionBackground: "var(--nextflow-light-green)", // Use site's light green

        // Scrollbar styling
        scrollbarThumbColor: "#d1d5db", // Light gray
        scrollbarThumbHoverColor: "#9ca3af", // Slightly darker on hover

        // Ensure clean, minimal appearance throughout
        uiSelectionBackground: "var(--nextflow-light-green)",
        uiSelectionForeground: "#1f2937", // Dark text for contrast
      },
    }),
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
