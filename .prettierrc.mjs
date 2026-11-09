// .prettierrc.mjs
/** @type {import("prettier").Config} */
export default {
  printWidth: 120,
  // Keep in sync with `compressHTML` in astro.config.mjs
  astroCompressHTML: true,
  plugins: ["prettier-plugin-astro"],
  overrides: [
    {
      files: "*.astro",
      options: {
        parser: "astro",
      },
    },
  ],
};
