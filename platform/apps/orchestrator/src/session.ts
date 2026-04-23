import { randomUUID } from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";
import { SignJWT, jwtVerify } from "jose";

/**
 * Session model: anonymous, JWT-in-cookie.
 *
 * On first contact we mint a sessionId (UUID), sign it into a JWT with a TTL
 * matching SESSION_TTL_MINUTES, set it as an HttpOnly cookie, and create
 * <DRIVES_ROOT>/<sessionId>/ — the "Little Drive". Every subsequent request
 * from that browser tag shares that drive.
 *
 * The cookie is the only state. Server restarts don't wipe sessions (the
 * JWT is self-contained and drives live on a mounted volume). A periodic
 * sweeper in index.ts removes drives older than the TTL.
 */

const COOKIE_NAME = "gb_session";
const TTL_MIN = Number(process.env.SESSION_TTL_MINUTES ?? 120);
const DRIVES_ROOT =
  process.env.DRIVES_ROOT ?? path.join(process.cwd(), "tmp", "drives");

// In dev we don't want to force HTTPS-only cookies.
const SECURE_COOKIE = process.env.NODE_ENV === "production";

function getSecret(): Uint8Array {
  const raw = process.env.SESSION_JWT_SECRET;
  if (!raw) {
    if (process.env.NODE_ENV === "production") {
      throw new Error("SESSION_JWT_SECRET must be set in production");
    }
    // Dev fallback — unstable across restarts but keeps local work moving.
    return new TextEncoder().encode("dev-only-insecure-secret-change-me");
  }
  return new TextEncoder().encode(raw);
}

export type SessionPayload = {
  sid: string;
  iat: number;
  exp: number;
};

export async function mintSession(): Promise<{
  sessionId: string;
  token: string;
  cookie: string;
  expiresAt: Date;
}> {
  const sid = randomUUID();
  const now = Math.floor(Date.now() / 1000);
  const exp = now + TTL_MIN * 60;

  const token = await new SignJWT({ sid })
    .setProtectedHeader({ alg: "HS256" })
    .setIssuedAt(now)
    .setExpirationTime(exp)
    .sign(getSecret());

  const cookieParts = [
    `${COOKIE_NAME}=${token}`,
    "Path=/",
    "HttpOnly",
    "SameSite=Lax",
    SECURE_COOKIE ? "Secure" : "",
    `Max-Age=${TTL_MIN * 60}`,
  ].filter(Boolean);

  await ensureDrive(sid);

  return {
    sessionId: sid,
    token,
    cookie: cookieParts.join("; "),
    expiresAt: new Date(exp * 1000),
  };
}

export async function readSession(cookieHeader: string | undefined): Promise<SessionPayload | null> {
  if (!cookieHeader) return null;
  const match = cookieHeader
    .split(";")
    .map((s) => s.trim())
    .find((s) => s.startsWith(`${COOKIE_NAME}=`));
  if (!match) return null;
  const token = match.slice(COOKIE_NAME.length + 1);
  try {
    const { payload } = await jwtVerify(token, getSecret());
    if (typeof payload.sid !== "string") return null;
    return payload as SessionPayload;
  } catch {
    return null;
  }
}

export function drivePath(sid: string): string {
  // Defence in depth: sid comes from a signed JWT, but we still refuse
  // anything that could escape the root.
  if (!/^[a-f0-9-]{36}$/.test(sid)) {
    throw new Error("invalid session id");
  }
  return path.join(DRIVES_ROOT, sid);
}

export async function ensureDrive(sid: string): Promise<string> {
  const dir = drivePath(sid);
  await fs.mkdir(dir, { recursive: true });
  return dir;
}

/**
 * Delete drives whose directory mtime is older than the TTL. Cheap and
 * correct: we update mtime on every write, so a quiet session expires when
 * its last activity falls off the horizon.
 */
export async function sweepExpired(): Promise<{ removed: number }> {
  let removed = 0;
  try {
    const entries = await fs.readdir(DRIVES_ROOT, { withFileTypes: true });
    const horizon = Date.now() - TTL_MIN * 60 * 1000;
    for (const entry of entries) {
      if (!entry.isDirectory()) continue;
      const full = path.join(DRIVES_ROOT, entry.name);
      const stat = await fs.stat(full);
      if (stat.mtimeMs < horizon) {
        await fs.rm(full, { recursive: true, force: true });
        removed++;
      }
    }
  } catch (err: unknown) {
    if ((err as NodeJS.ErrnoException).code !== "ENOENT") throw err;
  }
  return { removed };
}

export { DRIVES_ROOT, COOKIE_NAME };
