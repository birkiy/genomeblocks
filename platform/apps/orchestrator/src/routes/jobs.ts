import type { FastifyInstance } from "fastify";
import path from "node:path";
import fs from "node:fs/promises";
import { readSession, drivePath } from "../session";
import { runCompute, ComputeError } from "../runner";

/**
 * Whitelist of compute jobs. Each entry declares:
 *   - the required input fields (filenames relative to the session drive)
 *   - the Python-side job name (passed to compute.cli)
 *
 * Anything not listed here returns 404 — we never pass user input through
 * to the Python dispatcher as-is.
 */
const JOBS: Record<string, { fileFields: string[]; timeoutMs: number }> = {
  venn: { fileFields: ["a", "b", "c"], timeoutMs: 30_000 },
};

export async function jobRoutes(app: FastifyInstance) {
  app.post("/jobs/:name", async (req, reply) => {
    const { name } = req.params as { name: string };
    const spec = JOBS[name];
    if (!spec) return reply.status(404).send({ error: `unknown job: ${name}` });

    const session = await readSession(req.headers.cookie);
    if (!session) return reply.status(401).send({ error: "no session" });

    const body = (req.body ?? {}) as Record<string, unknown>;
    const drive = drivePath(session.sid);

    // Resolve every declared file field against the session drive, ensuring
    // the resolved path still lives under it.
    const resolved: Record<string, string> = {};
    for (const field of spec.fileFields) {
      const value = body[field];
      if (typeof value !== "string") {
        return reply.status(400).send({ error: `missing field: ${field}` });
      }
      const abs = path.resolve(drive, value);
      if (!abs.startsWith(drive + path.sep)) {
        return reply.status(400).send({ error: `bad path for ${field}` });
      }
      try {
        await fs.access(abs);
      } catch {
        return reply.status(404).send({ error: `file not found in session drive: ${value}` });
      }
      resolved[field] = abs;
    }

    try {
      const run = await runCompute(name, resolved, {
        cwd: drive,
        timeoutMs: spec.timeoutMs,
      });
      if (run.stderr) req.log.debug({ stderr: run.stderr }, "compute stderr");
      return { job: name, result: run.result, durationMs: run.durationMs };
    } catch (err) {
      if (err instanceof ComputeError) {
        req.log.error({ stderr: err.stderr, code: err.exitCode }, err.message);
        return reply.status(500).send({
          error: err.message,
          // Forward the stderr tail to help users debug malformed inputs.
          detail: err.stderr.slice(-400),
        });
      }
      throw err;
    }
  });
}
