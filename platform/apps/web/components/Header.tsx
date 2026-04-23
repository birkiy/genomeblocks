import Link from "next/link";

/**
 * Shared header across docs and dashboard. Keeps the two surfaces feeling
 * like one product even though they render different content below.
 */
export function Header() {
  return (
    <header
      style={{
        borderBottom: "1px solid var(--color-border)",
        background: "var(--color-bg-raised)",
      }}
    >
      <div
        className="flex items-center justify-between"
        style={{
          maxWidth: 1200,
          margin: "0 auto",
          padding: "var(--space-4) var(--space-6)",
        }}
      >
        <Link
          href="/"
          style={{
            fontWeight: 800,
            fontSize: "var(--text-xl)",
            color: "var(--brand)",
            letterSpacing: "-0.02em",
          }}
        >
          genomeblocks
        </Link>
        <nav className="flex items-center gap-6">
          <Link href="/docs" style={navLinkStyle}>Docs</Link>
          <Link href="/app" style={navLinkStyle}>Dashboard</Link>
          <a
            href="https://github.com/birkiy/genomeblocks"
            style={navLinkStyle}
            target="_blank"
            rel="noreferrer"
          >
            GitHub
          </a>
        </nav>
      </div>
    </header>
  );
}

const navLinkStyle: React.CSSProperties = {
  fontWeight: 600,
  fontSize: "var(--text-sm)",
  color: "var(--color-text-muted)",
};
