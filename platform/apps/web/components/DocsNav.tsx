import Link from "next/link";
import type { DocEntry } from "@/lib/docs";

export function DocsNav({ entries, current }: { entries: DocEntry[]; current: string }) {
  return (
    <nav
      style={{
        position: "sticky",
        top: 0,
        alignSelf: "start",
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
        Documentation
      </div>
      <ul style={{ listStyle: "none", margin: 0, padding: 0 }}>
        {entries.map((entry) => {
          const active = entry.href === current;
          return (
            <li key={entry.href} style={{ margin: "2px 0" }}>
              <Link
                href={entry.href}
                style={{
                  display: "block",
                  padding: "var(--space-2) var(--space-3)",
                  borderRadius: "var(--radius-sm)",
                  fontSize: "var(--text-sm)",
                  fontWeight: active ? 700 : 500,
                  color: active ? "var(--brand)" : "var(--color-text)",
                  background: active ? "var(--brand-subtle)" : "transparent",
                }}
              >
                {entry.title}
              </Link>
            </li>
          );
        })}
      </ul>
    </nav>
  );
}
