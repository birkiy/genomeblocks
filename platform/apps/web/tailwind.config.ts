import type { Config } from "tailwindcss";

// Tailwind is used for layout utilities only; all colors/typography come from
// the Genomeblocks design tokens in styles/tokens.css. We expose the brand
// token names as Tailwind colors so `bg-brand` etc. work, but the source of
// truth is the CSS variables.
const config: Config = {
  content: [
    "./app/**/*.{ts,tsx,mdx}",
    "./components/**/*.{ts,tsx}",
    "./content/**/*.{md,mdx}",
  ],
  theme: {
    extend: {
      colors: {
        brand: "var(--brand)",
        "brand-hover": "var(--brand-hover)",
        "brand-subtle": "var(--brand-subtle)",
        bg: "var(--color-bg)",
        "bg-raised": "var(--color-bg-raised)",
        "bg-sunken": "var(--color-bg-sunken)",
        border: "var(--color-border)",
        "border-strong": "var(--color-border-strong)",
        ink: "var(--color-text)",
        "ink-muted": "var(--color-text-muted)",
        "ink-subtle": "var(--color-text-subtle)",
      },
      fontFamily: {
        sans: ["Nunito", "Helvetica Neue", "sans-serif"],
        mono: ["Varela Round", "Helvetica Neue", "sans-serif"],
      },
      borderRadius: {
        sm: "6px",
        md: "10px",
        lg: "16px",
      },
      boxShadow: {
        sm: "var(--shadow-sm)",
        md: "var(--shadow-md)",
        lg: "var(--shadow-lg)",
      },
    },
  },
  plugins: [],
};

export default config;
