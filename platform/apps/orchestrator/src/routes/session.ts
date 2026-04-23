import type { FastifyInstance } from "fastify";
import { mintSession, readSession, ensureDrive } from "../session";

/**
 * POST /session — idempotent "give me a session" endpoint.
 *
 * If the caller already has a valid cookie, we just touch the drive and
 * echo the existing session. If not, mint a new one. Either way, the
 * response shape is the same so the client doesn't have to branch.
 */
export async function sessionRoutes(app: FastifyInstance) {
  app.post("/session", async (req, reply) => {
    const existing = await readSession(req.headers.cookie);
    if (existing) {
      await ensureDrive(existing.sid);
      return {
        sessionId: existing.sid,
        expiresAt: new Date(existing.exp * 1000).toISOString(),
      };
    }

    const session = await mintSession();
    reply.header("set-cookie", session.cookie);
    return {
      sessionId: session.sessionId,
      expiresAt: session.expiresAt.toISOString(),
    };
  });
}
