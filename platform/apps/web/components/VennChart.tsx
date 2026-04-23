"use client";

/**
 * Static 3-set Venn SVG. We don't do area-proportional geometry — just the
 * three overlapping circles with counts in each region. Keeps bundle size
 * tiny; zero chart lib dependency.
 *
 * The backend gives us disjoint counts (only_A, only_B, only_C, ab_only,
 * ac_only, bc_only, abc). Each corresponds to one labeled region below.
 */

export type VennCounts = {
  only_a: number;
  only_b: number;
  only_c: number;
  ab_only: number;
  ac_only: number;
  bc_only: number;
  abc: number;
  a_total: number;
  b_total: number;
  c_total: number;
};

type Props = {
  counts: VennCounts;
  labels: { a: string; b: string; c: string };
};

export function VennChart({ counts, labels }: Props) {
  return (
    <svg
      viewBox="0 0 440 360"
      style={{ width: "100%", height: "auto", maxWidth: 560 }}
      role="img"
      aria-label="Three-set Venn diagram of region overlaps"
    >
      {/* Circles — positioned by eye for a clean 3-set layout. */}
      <circle cx="160" cy="150" r="110" fill="var(--venn-a)" fillOpacity="0.35" />
      <circle cx="280" cy="150" r="110" fill="var(--venn-b)" fillOpacity="0.35" />
      <circle cx="220" cy="240" r="110" fill="var(--venn-c)" fillOpacity="0.35" />

      {/* Set labels */}
      <text x="70"  y="70"  style={labelStyle} fill="var(--venn-a)">{labels.a}</text>
      <text x="330" y="70"  style={labelStyle} fill="var(--venn-b)" textAnchor="end">{labels.b}</text>
      <text x="220" y="350" style={labelStyle} fill="var(--venn-c)" textAnchor="middle">{labels.c}</text>

      {/* Region counts — 7 disjoint regions */}
      <text x="110" y="130" style={countStyle}>{counts.only_a}</text>
      <text x="330" y="130" style={countStyle}>{counts.only_b}</text>
      <text x="220" y="290" style={countStyle}>{counts.only_c}</text>
      <text x="220" y="105" style={countStyle}>{counts.ab_only}</text>
      <text x="155" y="230" style={countStyle}>{counts.ac_only}</text>
      <text x="285" y="230" style={countStyle}>{counts.bc_only}</text>
      <text x="220" y="185" style={countStyle}>{counts.abc}</text>
    </svg>
  );
}

const labelStyle: React.CSSProperties = {
  fontFamily: "var(--font-sans)",
  fontWeight: 800,
  fontSize: 16,
};

const countStyle: React.CSSProperties = {
  fontFamily: "var(--font-sans)",
  fontWeight: 700,
  fontSize: 18,
  fill: "var(--color-text)",
  textAnchor: "middle",
};
