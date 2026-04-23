import Link from "next/link";

/**
 * Dashboard shell. Distinct left-rail from the docs nav — this is "the tool"
 * rather than "the reference". Shares the global Header, so brand and
 * cross-surface nav stay consistent.
 */
export default function DashboardLayout({ children }: { children: React.ReactNode }) {
  return (
    <div
      className="grid"
      style={{
        gridTemplateColumns: "220px 1fr",
        maxWidth: 1200,
        margin: "0 auto",
      }}
    >
      <aside
        style={{
          padding: "var(--space-8) var(--space-4)",
          borderRight: "1px solid var(--color-border)",
          minHeight: "100vh",
        }}
      >
        <div
          style={{
            fontSize: "var(--text-xs)",
            textTransform: "uppercase",
            letterSpacing: "0.08em",
            color: "var(--color-text-subtle)",
            marginBottom: "var(--space-3)",
          }}
        >
          Dashboard
        </div>
        <ul style={{ listStyle: "none", margin: 0, padding: 0 }}>
          <NavLink href="/app">Overview</NavLink>
          <NavLink href="/app/playground/venn">Venn — region overlap</NavLink>
        </ul>
      </aside>
      <section style={{ padding: "var(--space-10) var(--space-8)" }}>
        {children}
      </section>
    </div>
  );
}

function NavLink({ href, children }: { href: string; children: React.ReactNode }) {
  return (
    <li style={{ margin: "2px 0" }}>
      <Link
        href={href}
        style={{
          display: "block",
          padding: "var(--space-2) var(--space-3)",
          borderRadius: "var(--radius-sm)",
          fontSize: "var(--text-sm)",
          fontWeight: 500,
          color: "var(--color-text)",
        }}
      >
        {children}
      </Link>
    </li>
  );
}
