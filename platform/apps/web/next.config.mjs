import createMDX from "@next/mdx";
import remarkGfm from "remark-gfm";
import rehypeSlug from "rehype-slug";

/**
 * Same-domain API: rewrite /api/* to the orchestrator URL from env.
 * In dev, point NEXT_PUBLIC_API_ORIGIN at http://localhost:8080.
 * In prod, it's the Fly.io URL (https://<app>.fly.dev).
 *
 * Because this is an edge rewrite, the browser only ever sees one origin —
 * the Next.js app. Cookies/JWT flow cleanly with no CORS configuration.
 */
const API_ORIGIN =
  process.env.NEXT_PUBLIC_API_ORIGIN?.replace(/\/$/, "") ?? "http://localhost:8080";

/** @type {import('next').NextConfig} */
const nextConfig = {
  pageExtensions: ["ts", "tsx", "md", "mdx"],
  reactStrictMode: true,
  async rewrites() {
    return [
      {
        source: "/api/:path*",
        destination: `${API_ORIGIN}/:path*`,
      },
    ];
  },
};

const withMDX = createMDX({
  extension: /\.mdx?$/,
  options: {
    remarkPlugins: [remarkGfm],
    rehypePlugins: [rehypeSlug],
  },
});

export default withMDX(nextConfig);
