import { defineEcConfig } from "astro-expressive-code";
import { pluginLineNumbers } from "@expressive-code/plugin-line-numbers";

export default defineEcConfig({
  // Use light theme for consistency with site design
  themes: ["github-light"],

  // Enable advanced features
  plugins: [pluginLineNumbers()],

  // Disable Expressive Code's built-in copy button and enable line numbers
  defaultProps: {
    showCopyToClipboardButton: false,
    showLineNumbers: true,
  },

  // Comprehensive styling to match site's clean, minimal design
  styleOverrides: {
    // Core background and text styling - pure white like site
    codeBackground: "#ffffff",
    codeForeground: "#24292f", // GitHub light theme text color

    // Typography - match site's exact monospace stack
    codeFontFamily: "Menlo, Monaco, Consolas, 'Courier New', monospace",
    codeFontSize: "1.25rem", // 20px - readable code size, smaller than before
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

    // Larger UI text for better readability of frame titles
    uiFontSize: "1.2rem", // Increased from default 0.9rem for more prominent titles
  },
});
