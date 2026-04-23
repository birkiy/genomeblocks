import type { Metadata } from "next";
import "./globals.css";
import { Header } from "@/components/Header";

export const metadata: Metadata = {
  title: {
    default: "genomeblocks",
    template: "%s — genomeblocks",
  },
  description:
    "Fluent building blocks for regulatory genomics — from peaks to chromatin networks in a handful of expressive chained calls.",
};

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="en">
      <body>
        <Header />
        <main>{children}</main>
      </body>
    </html>
  );
}
