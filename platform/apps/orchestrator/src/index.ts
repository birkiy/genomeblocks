import Fastify from "fastify";
import cors from "@fastify/cors";
import cookie from "@fastify/cookie";
import multipart from "@fastify/multipart";
import { sessionRoutes } from "./routes/session";
import { uploadRoutes } from "./routes/upload";
import { jobRoutes } from "./routes/jobs";
import { sweepExpired, DRIVES_ROOT } from "./session";
import { UPLOAD_LIMIT_BYTES } from "./routes/upload";

const PORT = Number(process.env.PORT ?? 8080);
const HOST = process.env.HOST ?? "0.0.0.0";

async function main() {
  const app = Fastify({
    logger: {
      level: process.env.LOG_LEVEL ?? "info",
      transport:
        process.env.NODE_ENV !== "production"
          ? { target: "pino-pretty", options: { singleLine: true } }
          : undefined,
    },
  });

  // In prod, requests arrive via the Vercel edge rewrite (same-origin from
  // the browser's perspective). In dev, the Next.js rewrite does the same.
  // CORS only matters if you hit the orchestrator directly from another
  // origin — kept permissive-with-credentials for that case.
  await app.register(cors, {
    origin: process.env.CORS_ORIGIN ?? true,
    credentials: true,
  });
  await app.register(cookie);
  await app.register(multipart, {
    limits: { fileSize: UPLOAD_LIMIT_BYTES, files: 10 },
  });

  app.get("/healthz", async () => ({ ok: true, drivesRoot: DRIVES_ROOT }));

  await app.register(sessionRoutes);
  await app.register(uploadRoutes);
  await app.register(jobRoutes);

  // Periodic sweeper. Runs every 10 minutes; cheap (one readdir + stat per
  // session). Fly's scale-to-zero means a stopped machine can't sweep — but
  // that's fine: on cold start the next sweep picks up any stragglers.
  const sweeperMs = 10 * 60 * 1000;
  const sweeper = setInterval(async () => {
    try {
      const { removed } = await sweepExpired();
      if (removed > 0) app.log.info({ removed }, "swept expired drives");
    } catch (err) {
      app.log.error({ err }, "sweep failed");
    }
  }, sweeperMs);
  sweeper.unref();

  await app.listen({ port: PORT, host: HOST });
  app.log.info(`orchestrator listening on http://${HOST}:${PORT}`);
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
