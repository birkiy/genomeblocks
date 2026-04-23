"use client";

import { useState } from "react";
import { FileDrop } from "@/components/FileDrop";
import { VennChart, type VennCounts } from "@/components/VennChart";
import { ensureSession, uploadFiles, runJob } from "@/lib/api";

type Status =
  | { kind: "idle" }
  | { kind: "working"; step: string }
  | { kind: "done"; counts: VennCounts; labels: { a: string; b: string; c: string }; durationMs: number }
  | { kind: "error"; message: string };

export default function VennPlayground() {
  const [a, setA] = useState<File | null>(null);
  const [b, setB] = useState<File | null>(null);
  const [c, setC] = useState<File | null>(null);
  const [status, setStatus] = useState<Status>({ kind: "idle" });

  const ready = a && b && c && status.kind !== "working";

  async function onRun() {
    if (!a || !b || !c) return;
    try {
      setStatus({ kind: "working", step: "Starting session…" });
      await ensureSession();

      setStatus({ kind: "working", step: "Uploading files to your session drive…" });
      const uploaded = await uploadFiles([a, b, c]);
      const [fa, fb, fc] = uploaded.files;

      setStatus({ kind: "working", step: "Running Python overlap job…" });
      const res = await runJob<VennCounts>("venn", {
        a: fa.name,
        b: fb.name,
        c: fc.name,
      });

      setStatus({
        kind: "done",
        counts: res.result,
        labels: { a: a.name, b: b.name, c: c.name },
        durationMs: res.durationMs,
      });
    } catch (err) {
      setStatus({ kind: "error", message: err instanceof Error ? err.message : String(err) });
    }
  }

  return (
    <div>
      <h1>Venn — region overlap</h1>
      <p style={{ color: "var(--color-text-muted)", maxWidth: 620 }}>
        Upload three BED-format files (or narrowPeak — we only read the first
        three columns). The Python compute sidecar computes the 3-way overlap
        and returns disjoint region counts.
      </p>

      <div
        className="grid md:grid-cols-3 gap-4"
        style={{ marginTop: "var(--space-6)" }}
      >
        <FileDrop label="Set A" accept=".bed,.narrowPeak,.txt" onChange={setA} />
        <FileDrop label="Set B" accept=".bed,.narrowPeak,.txt" onChange={setB} />
        <FileDrop label="Set C" accept=".bed,.narrowPeak,.txt" onChange={setC} />
      </div>

      <div style={{ marginTop: "var(--space-6)" }}>
        <button
          onClick={onRun}
          disabled={!ready}
          style={{
            padding: "var(--space-3) var(--space-6)",
            background: ready ? "var(--brand)" : "var(--color-bg-sunken)",
            color: ready ? "white" : "var(--color-text-muted)",
            border: 0,
            borderRadius: "var(--radius-md)",
            fontWeight: 700,
            fontSize: "var(--text-base)",
            cursor: ready ? "pointer" : "not-allowed",
            boxShadow: ready ? "var(--shadow-sm)" : "none",
            transition: "all var(--duration-fast) var(--ease-out)",
          }}
        >
          {status.kind === "working" ? status.step : "Compute overlap"}
        </button>
      </div>

      {status.kind === "error" && (
        <div
          style={{
            marginTop: "var(--space-6)",
            padding: "var(--space-4)",
            borderRadius: "var(--radius-md)",
            background: "#fef2f2",
            color: "#991b1b",
            fontFamily: "var(--font-code)",
            fontSize: "var(--text-sm)",
          }}
        >
          {status.message}
        </div>
      )}

      {status.kind === "done" && (
        <div style={{ marginTop: "var(--space-8)" }}>
          <h2>Result</h2>
          <p
            style={{
              color: "var(--color-text-muted)",
              fontSize: "var(--text-sm)",
              marginBottom: "var(--space-4)",
            }}
          >
            Computed in {status.durationMs} ms. A={status.counts.a_total} regions,
            B={status.counts.b_total}, C={status.counts.c_total}.
          </p>
          <VennChart counts={status.counts} labels={status.labels} />
        </div>
      )}
    </div>
  );
}
