import type { MDXComponents } from "mdx/types";

/**
 * Map MDX elements to our styled equivalents. Docs MDX authors write plain
 * markdown; this file decides how each element looks.
 */
export function useMDXComponents(components: MDXComponents): MDXComponents {
  return { ...components };
}
