import type { FastifyInstance } from "fastify";
import fs from "node:fs";
import path from "node:path";
import { pipeline } from "node:stream/promises";
import { readSession, ensureDrive, mintSession } from "../session";

const MAX_BYTES = 50 * 1024 * 1024; // 50 MB per request

/**
 * POST /upload — multipart, fields named "files".
 *
 * Writes each uploaded file into the session's "Little Drive". Filenames are
 * sanitized so nothing can escape the session directory via traversal
 * (`../`) or absolute paths.
 *
 * If the caller doesn't have a session yet, we mint one and set the cookie
 * in the response — saves the client from having to make two round trips.
 */
export async function uploadRoutes(app: FastifyInstance) {
  app.post("/upload", async (req, reply) => {
    let session = await readSession(req.headers.cookie);
    if (!session) {
      const minted = await mintSession();
      reply.header("set-cookie", minted.cookie);
      session = { sid: minted.sessionId, iat: 0, exp: 0 };
    }

    const driveDir = await ensureDrive(session.sid);
    const saved: Array<{ name: string; size: number }> = [];

    const parts = req.parts();
    for await (const part of parts) {
      if (part.type !== "file") continue;

      const safeName = sanitizeFilename(part.filename);
      const dest = path.join(driveDir, safeName);

      // Stream to disk. Fastify's file stream exposes `truncated` after the
      // byte limit is hit; we check that at the end and refuse the request.
      const writeStream = fs.createWriteStream(dest);
      await pipeline(part.file, writeStream);

      if (part.file.truncated) {
        await fs.promises.unlink(dest).catch(() => {});
        return reply.status(413).send({ error: `file '${safeName}' exceeds ${MAX_BYTES} bytes` });
      }

      const stat = await fs.promises.stat(dest);
      saved.push({ name: safeName, size: stat.size });
    }

    if (saved.length === 0) {
      return reply.status(400).send({ error: "no files uploaded" });
    }

    return { files: saved };
  });
}

/** Strip directory parts and anything non-printable. Keeps the original extension. */
function sanitizeFilename(raw: string): string {
  const base = path.basename(raw);
  const cleaned = base.replace(/[^\w.\- ]+/g, "_").trim();
  return cleaned.length ? cleaned : `upload_${Date.now()}`;
}

export const UPLOAD_LIMIT_BYTES = MAX_BYTES;
