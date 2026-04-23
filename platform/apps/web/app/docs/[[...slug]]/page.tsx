import { notFound } from "next/navigation";
import { listDocs } from "@/lib/docs";

// Statically generate a page per MDX file at build time.
export async function generateStaticParams() {
  const entries = await listDocs();
  return entries.map((e) => ({ slug: e.slug.length ? e.slug : undefined }));
}

type Params = { slug?: string[] };

export default async function DocPage({ params }: { params: Params }) {
  const slug = params.slug ?? [];
  const segments = slug.length === 0 ? ["index"] : slug;

  // Dynamic MDX import so Next's MDX loader wraps the module in our
  // configured remark/rehype pipeline.
  let Mod: { default: React.ComponentType };
  try {
    Mod = await import(`@/content/docs/${segments.join("/")}.mdx`);
  } catch {
    notFound();
  }

  const MDX = Mod.default;
  return <MDX />;
}
