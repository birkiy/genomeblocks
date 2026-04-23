import Link from "next/link";

export default function HomePage() {
  return (
    <div
      style={{
        maxWidth: 1100,
        margin: "0 auto",
        padding: "var(--space-16) var(--space-6)",
      }}
    >
      <h1 style={{ fontSize: "var(--text-5xl)" }}>
        Fluent building blocks for regulatory genomics.
      </h1>
      <p
        style={{
          fontSize: "var(--text-lg)",
          color: "var(--color-text-muted)",
          maxWidth: 620,
        }}
      >
        From peaks to chromatin networks in a handful of expressive chained
        calls. Read the docs, then open the playground and run the same
        pipelines in your browser.
      </p>

      <div
        className="flex gap-3"
        style={{ marginTop: "var(--space-8)" }}
      >
        <Link href="/docs" style={primaryBtn}>Read the docs</Link>
        <Link href="/app" style={secondaryBtn}>Open the dashboard</Link>
      </div>

      <section
        className="grid md:grid-cols-3 gap-6"
        style={{ marginTop: "var(--space-16)" }}
      >
        <Card
          title="Static docs"
          body="Markdown + MDX. Pre-rendered at build time, served from the edge."
        />
        <Card
          title="Live playground"
          body="React dashboard backed by real Python compute on a sandboxed per-session drive."
        />
        <Card
          title="One domain"
          body="Docs and dashboard share the same origin, theme, and navigation. Seamless switch."
        />
      </section>
    </div>
  );
}

function Card({ title, body }: { title: string; body: string }) {
  return (
    <div
      style={{
        background: "var(--color-bg-raised)",
        border: "1px solid var(--color-border)",
        borderRadius: "var(--radius-md)",
        padding: "var(--space-6)",
        boxShadow: "var(--shadow-sm)",
      }}
    >
      <h3>{title}</h3>
      <p style={{ color: "var(--color-text-muted)", margin: 0 }}>{body}</p>
    </div>
  );
}

const primaryBtn: React.CSSProperties = {
  display: "inline-block",
  padding: "var(--space-3) var(--space-5)",
  background: "var(--brand)",
  color: "white",
  borderRadius: "var(--radius-md)",
  fontWeight: 600,
  boxShadow: "var(--shadow-sm)",
};

const secondaryBtn: React.CSSProperties = {
  display: "inline-block",
  padding: "var(--space-3) var(--space-5)",
  background: "var(--color-bg-raised)",
  color: "var(--brand)",
  border: "1px solid var(--color-border-strong)",
  borderRadius: "var(--radius-md)",
  fontWeight: 600,
};
