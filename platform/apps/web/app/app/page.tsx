import Link from "next/link";

export default function DashboardHome() {
  return (
    <div>
      <h1>Dashboard</h1>
      <p style={{ color: "var(--color-text-muted)", maxWidth: 620 }}>
        Interactive playgrounds that run real compute against a sandboxed
        session drive. Nothing you upload is visible to anyone else; files are
        swept after your session expires.
      </p>

      <div
        className="grid md:grid-cols-2 gap-6"
        style={{ marginTop: "var(--space-6)" }}
      >
        <PlaygroundCard
          href="/app/playground/venn"
          title="Venn — region overlap"
          body="Upload three BED/narrowPeak files; see a 3-way overlap count rendered as a Venn diagram. Computed server-side with Python."
        />
      </div>
    </div>
  );
}

function PlaygroundCard({
  href,
  title,
  body,
}: {
  href: string;
  title: string;
  body: string;
}) {
  return (
    <Link
      href={href}
      style={{
        display: "block",
        background: "var(--color-bg-raised)",
        border: "1px solid var(--color-border)",
        borderRadius: "var(--radius-md)",
        padding: "var(--space-6)",
        boxShadow: "var(--shadow-sm)",
        color: "inherit",
      }}
    >
      <h3 style={{ marginTop: 0 }}>{title}</h3>
      <p style={{ color: "var(--color-text-muted)", margin: 0 }}>{body}</p>
    </Link>
  );
}
