import { spawn } from "node:child_process";
import path from "node:path";

/**
 * Node → Python bridge. We don't pipe structured input over stdin — instead
 * we pass the job name + the caller-provided JSON as CLI args, because it
 * makes every job reproducible from a shell:
 *
 *   python -m compute.cli venn '{"a":"x.bed","b":"y.bed","c":"z.bed"}'
 *
 * Python's stdout must be exactly one JSON document. Anything on stderr is
 * treated as log output (piped into our pino logger by the route handler).
 *
 * The `cwd` passed in scopes all relative paths to the session drive — so
 * Python can't accidentally read outside it unless the user explicitly
 * supplies an absolute path (which we sanitize out in the upload/path
 * validation layer before getting here).
 */

const PYTHON_BIN = process.env.PYTHON_BIN ?? "python3";

// Path to the `compute` package. In the container, WORKDIR is /app and the
// compute package sits at /app/compute, so `python -m compute.cli` works
// when PYTHONPATH includes /app/compute.
const COMPUTE_DIR =
  process.env.COMPUTE_DIR ?? path.resolve(process.cwd(), "../../compute");

export type RunResult<T = unknown> = {
  result: T;
  durationMs: number;
  stderr: string;
};

export class ComputeError extends Error {
  constructor(message: string, public readonly stderr: string, public readonly exitCode: number) {
    super(message);
    this.name = "ComputeError";
  }
}

export async function runCompute<T = unknown>(
  jobName: string,
  inputs: Record<string, unknown>,
  opts: { cwd: string; timeoutMs?: number }
): Promise<RunResult<T>> {
  const timeoutMs = opts.timeoutMs ?? 60_000;

  return new Promise((resolve, reject) => {
    const started = Date.now();
    const args = ["-m", "compute.cli", jobName, JSON.stringify(inputs)];

    const env = {
      ...process.env,
      // Ensure the compute package is importable even when pip install wasn't
      // done (e.g. local dev without editable install).
      PYTHONPATH: [COMPUTE_DIR, process.env.PYTHONPATH].filter(Boolean).join(":"),
      // Keep job cwd inside the session drive for any relative reads.
      PWD: opts.cwd,
    };

    const child = spawn(PYTHON_BIN, args, {
      cwd: opts.cwd,
      env,
      stdio: ["ignore", "pipe", "pipe"],
    });

    let stdout = "";
    let stderr = "";
    child.stdout.on("data", (chunk: Buffer) => { stdout += chunk.toString(); });
    child.stderr.on("data", (chunk: Buffer) => { stderr += chunk.toString(); });

    const timer = setTimeout(() => {
      child.kill("SIGKILL");
      reject(new ComputeError(`compute job '${jobName}' timed out after ${timeoutMs}ms`, stderr, -1));
    }, timeoutMs);

    child.on("error", (err) => {
      clearTimeout(timer);
      reject(new ComputeError(`failed to spawn python: ${err.message}`, stderr, -1));
    });

    child.on("close", (code) => {
      clearTimeout(timer);
      const durationMs = Date.now() - started;
      if (code !== 0) {
        return reject(
          new ComputeError(
            `compute job '${jobName}' exited with code ${code}`,
            stderr,
            code ?? -1
          )
        );
      }
      try {
        const result = JSON.parse(stdout) as T;
        resolve({ result, durationMs, stderr });
      } catch (err) {
        reject(
          new ComputeError(
            `compute job '${jobName}' returned invalid JSON: ${(err as Error).message}`,
            stderr || stdout,
            0
          )
        );
      }
    });
  });
}
