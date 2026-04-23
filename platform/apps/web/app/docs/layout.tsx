import { listDocs } from "@/lib/docs";
import { DocsNav } from "@/components/DocsNav";
import { headers } from "next/headers";

export default async function DocsLayout({ children }: { children: React.ReactNode }) {
  const entries = await listDocs();
  // next/headers gives us the current URL for highlighting the active link.
  const pathname = headers().get("x-invoke-path") ?? "/docs";

  return (
    <div
      className="grid"
      style={{
        gridTemplateColumns: "240px 1fr",
        maxWidth: 1200,
        margin: "0 auto",
      }}
    >
      <DocsNav entries={entries} current={pathname} />
      <article
        className="prose"
        style={{ padding: "var(--space-10) var(--space-8)" }}
      >
        {children}
      </article>
    </div>
  );
}
