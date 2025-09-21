import { defineEcConfig } from "astro-expressive-code";
import { pluginLineNumbers } from "@expressive-code/plugin-line-numbers";

export default defineEcConfig({
  // Use light theme for consistency with site design
  themes: ["github-light"],

  // Enable advanced features
  plugins: [pluginLineNumbers()],

  // Disable Expressive Code's built-in copy button and enable line numbers
  defaultProps: {
    showCopyToClipboardButton: true,
    showLineNumbers: true,
  },

  // Configuration with proper spacing for text marker labels
  styleOverrides: {
    // Increase inline padding to prevent label overlap with code
    codePaddingInline: "3rem", // Increased from default 1.35rem to accommodate labels

    // Use defaults for text markers
    textMarkers: {},
  },
});
