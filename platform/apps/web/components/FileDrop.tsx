"use client";

import { useCallback, useRef, useState } from "react";

type Props = {
  label: string;
  accept?: string;
  onChange: (file: File | null) => void;
};

/**
 * Single-file drop zone. The venn playground uses three of these, one per
 * set. Kept dumb — the parent owns the file state.
 */
export function FileDrop({ label, accept, onChange }: Props) {
  const inputRef = useRef<HTMLInputElement>(null);
  const [file, setFile] = useState<File | null>(null);
  const [dragging, setDragging] = useState(false);

  const onFiles = useCallback(
    (list: FileList | null) => {
      const f = list?.[0] ?? null;
      setFile(f);
      onChange(f);
    },
    [onChange]
  );

  return (
    <div
      onClick={() => inputRef.current?.click()}
      onDragOver={(e) => {
        e.preventDefault();
        setDragging(true);
      }}
      onDragLeave={() => setDragging(false)}
      onDrop={(e) => {
        e.preventDefault();
        setDragging(false);
        onFiles(e.dataTransfer.files);
      }}
      style={{
        border: `2px dashed ${dragging ? "var(--brand)" : "var(--color-border-strong)"}`,
        borderRadius: "var(--radius-md)",
        padding: "var(--space-6)",
        textAlign: "center",
        cursor: "pointer",
        background: dragging ? "var(--brand-subtle)" : "var(--color-bg-raised)",
        transition: "all var(--duration-fast) var(--ease-out)",
      }}
    >
      <div style={{ fontWeight: 700, marginBottom: 4 }}>{label}</div>
      <div
        style={{
          fontSize: "var(--text-sm)",
          color: file ? "var(--color-text)" : "var(--color-text-muted)",
        }}
      >
        {file ? file.name : "Drop a .bed / .narrowPeak file or click to choose"}
      </div>
      <input
        ref={inputRef}
        type="file"
        accept={accept}
        style={{ display: "none" }}
        onChange={(e) => onFiles(e.target.files)}
      />
    </div>
  );
}
