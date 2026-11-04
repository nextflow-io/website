const space = 4;

const getSpace = (increment) => `${increment * space}px`;

const getSpacing = () => {
  const increments = 112;
  const values = [];
  for (let i = 1; i < increments; i++) {
    values.push({ [i]: getSpace(i) });
  }
  return values.reduce((acc, curr) => ({ ...acc, ...curr }), {});
};

export const screens = {
  "2xs": "200px",
  xs: "340px",
  sm: "640px",
  md: "768px",
  lg: "1024px",
  xl: "1240px",
  xxl: "1500px",
};

export default {
  important: true,
  content: ["./src/**/*.{astro,html,js,jsx,ts,tsx,vue,svelte}"],
  theme: {
    screens,
    fontFamily: {
      display: ["Degular", "Helvetica", "sans-serif"],
      body: ["Inter", "-apple-system", "BlinkMacSystemFont", "Segoe UI", "Roboto", "Helvetica", "Arial", "sans-serif"],
      sans: ["Inter", "-apple-system", "BlinkMacSystemFont", "Segoe UI", "Roboto", "Helvetica", "Arial", "sans-serif"],
      mono: ["Menlo", "Monaco", "monospace", "sans-serif"],
    },
    fontSize: {
      xs: ["12px", "18px"],
      sm: ["14px", "21px"],
      base: ["16px", "24px"],
      lg: ["18px", "26px"],
      xl: ["24px", "32px"],
      "2xl": ["28px", "1"],
      "3xl": ["32px", "1"],
      "4xl": ["40px", "1"],
      "5xl": ["48px", "1"],
      "6xl": ["56px", "1"],
      "7xl": ["60px", "1"],
      "8xl": ["64px", "1"],
      "9xl": ["68px", "1"],
    },
    extend: {
      colors: {
        nextflow: {
          100: "#E2F7F3",
          200: "#B6ECE2",
          300: "#86E0CE",
          400: "#56D3BA",
          500: "#31C9AC",
          600: "#0DC09D",
          700: "#0CAE8E",
          800: "#0A967B",
          900: "#087F68",
          1000: "#065647",
          DEFAULT: "#0DC09D",
          "100-opacity": "rgba(226, 247, 243, 0.20)",
          "200-opacity": "rgba(182, 236, 226, 0.20)",
          "300-opacity": "rgba(134, 224, 206, 0.20)",
          "400-opacity": "rgba(86, 211, 186, 0.20)",
          "500-opacity": "rgba(49, 201, 172, 0.20)",
          "600-opacity": "rgba(13, 192, 157, 0.20)",
          "700-opacity": "rgba(12, 174, 142, 0.20)",
          "800-opacity": "rgba(10, 150, 123, 0.20)",
          "900-opacity": "rgba(8, 127, 104, 0.20)",
          "1000-opacity": "rgba(6, 86, 71, 0.20)",
        },
        brand: {
          100: "#F3F3F4",
          200: "#E8E7E9",
          300: "#D0CFD4",
          400: "#B9B7BE",
          500: "#A29FA8",
          600: "#8A8792",
          700: "#736F7D",
          800: "#5C5767",
          900: "#453F51",
          1000: "#2D273C",
          1100: "#160F26",
          DEFAULT: "#160F26",
          "100-opacity": "rgba(243, 243, 244, 0.20)",
          "200-opacity": "rgba(232, 231, 233, 0.20)",
          "300-opacity": "rgba(208, 207, 212, 0.20)",
          "400-opacity": "rgba(185, 183, 190, 0.20)",
          "500-opacity": "rgba(162, 159, 168, 0.20)",
          "600-opacity": "rgba(138, 135, 146, 0.20)",
          "700-opacity": "rgba(115, 111, 125, 0.20)",
          "800-opacity": "rgba(92, 87, 103, 0.20)",
          "900-opacity": "rgba(69, 63, 81, 0.20)",
          "1000-opacity": "rgba(45, 39, 60, 0.20)",
          "1100-opacity": "rgba(22, 15, 38, 0.20)",
        },
        gray: {
          100: "#F7F7F7",
          200: "#EAEBEB",
          300: "#DDDEDE",
          400: "#CFD0D1",
          500: "#C4C6C7",
          600: "#BABCBD",
          700: "#A8AAAB",
          800: "#919393",
          900: "#7B7B7B",
          1000: "#545555",
        },
        blue: {
          100: "#E8EBFC",
          200: "#C6CCF8",
          300: "#A1ABF3",
          400: "#7B89EE",
          500: "#5E6FEB",
          600: "#4256E7",
          700: "#3C4ED1",
          800: "#3443B4",
          900: "#2C3999",
          1000: "#1E2768",
          DEFAULT: "#4256E7",
        },
      },
      borderColor: {
        "brand-opacity": "var(--color-brand-900-opacity)",
      },
      borderRadius: {
        sm: "4px",
        md: "8px",
        lg: "16px",
        full: "9999px",
      },
      margin: {
        "1/12": "8.333333%",
        "2/12": "16.666667%",
        "3/12": "25%",
        "4/12": "33.333333%",
        "5/12": "41.666667%",
        "6/12": "50%",
        "7/12": "58.333333%",
        "8/12": "66.666667%",
        "9/12": "75%",
        "10/12": "83.333333%",
        "11/12": "91.666667%",
        "-1/12": "-8.333333%",
        "-2/12": "-16.666667%",
        "-3/12": "-25%",
        "-4/12": "-33.333333%",
        "-5/12": "-41.666667%",
        "-6/12": "-50%",
        "-7/12": "-58.333333%",
        "-8/12": "-66.666667%",
        "-9/12": "-75%",
        "-10/12": "-83.333333%",
        "-11/12": "-91.666667%",
      },
      flex: {
        "1/1": "0 0 100%",
        "1/12": "0 0 8.333333%",
        "2/12": "0 0 16.666667%",
        "3/12": "0 0 25%",
        "4/12": "0 0 33.333333%",
        "5/12": "0 0 41.666667%",
        "6/12": "0 0 50%",
        "7/12": "0 0 58.333333%",
        "8/12": "0 0 66.666667%",
        "9/12": "0 0 75%",
        "10/12": "0 0 83.333333%",
        "11/12": "0 0 91.666667%",
        "1/5": "0 0 20%",
      },
      height: {
        18: "4.5rem",
      },
      minWidth: {
        "1/12": "8.333333%",
        "2/12": "16.666667%",
        "3/12": "25%",
        "4/12": "33.333333%",
        "5/12": "41.666667%",
        "6/12": "50%",
        "7/12": "58.333333%",
        "8/12": "66.666667%",
        "9/12": "75%",
        "10/12": "83.333333%",
        "11/12": "91.666667%",
        "1/5": "20%",
      },
      width: {
        "1/7": "14.2857143%",
        "1/8": "12.5%",
      },
      spacing: getSpacing(),
      zIndex: {
        header: 100,
        sidebar: 150,
      },
      container: {
        center: true,
        screens: {
          sm: "1240px",
          md: "1240px",
          lg: "1240px",
          xl: "1240px",
        },
      },
    },
  },
  plugins: [
    function ({ addBase, theme }) {
      function extractColorVars(colorObj, colorGroup = "") {
        return Object.keys(colorObj).reduce((vars, colorKey) => {
          const value = colorObj[colorKey];

          let extension = `-${colorKey}`;
          if (colorKey === "DEFAULT") extension = "";

          const newVars =
            typeof value === "string"
              ? { [`--color${colorGroup}${extension}`]: value }
              : extractColorVars(value, `-${colorKey}`);

          return { ...vars, ...newVars };
        }, {});
      }

      const customVars = extractColorVars(theme("colors"));

      Object.entries(screens).forEach(([key, value]) => {
        customVars[`--${key}`] = value;
      });

      Object.keys(theme("fontSize")).forEach((key) => {
        customVars[`--size-${key}`] = theme("fontSize")[key][0];
        customVars[`--leading-${key}`] = theme("fontSize")[key][1];
      });
      addBase({
        ":root": customVars,
        p: {
          fontSize: theme("fontSize.base[0]"),
          lineHeight: theme("fontSize.base[1]"),
        },
      });
    },
  ],
};
