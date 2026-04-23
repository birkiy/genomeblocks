# Genomeblocks Platform

Unified web surface for [genomeblocks](https://github.com/birkiy/genomeblocks): static docs and an interactive dashboard/playground on the **same domain**, with a sandboxed per-session "Little Drive" and a Python compute sidecar.

## Architecture

Three-tier, decoupled, with Python used as a disposable compute worker.

```
┌──────────────────────────────────────────────────────────────┐
│                        genomeblocks.dev                      │
│                                                              │
│   /docs/*            /app/*                    /api/*        │
│   (MDX static)       (React dashboard)         (rewrite) ─┐  │
│                                                           │  │
│                       Next.js 14 on Vercel                │  │
└───────────────────────────────────────────────────────────┼──┘
                                                            │
                                                            ▼
                                          ┌─────────────────────────────┐
                                          │  Fastify orchestrator       │
                                          │  • JWT anon sessions        │
                                          │  • per-session "Little      │
                                          │    Drive" on volume          │
                                          │  • spawns Python compute    │
                                          └────────────┬────────────────┘
                                                       │ child_process.spawn
                                                       ▼
                                          ┌─────────────────────────────┐
                                          │  python -m compute.cli      │
                                          │  (pandas/numpy, stdout=JSON)│
                                          └─────────────────────────────┘
                                          (orchestrator + compute share
                                           one container image on Fly.io)
```

**Why this shape**
- **Static-first frontend.** Docs pre-render at build time; dashboard hydrates on demand. One Next.js app serves both, so they share the nav, theme, and URL space.
- **Node owns I/O and state.** Fastify handles the concurrent browser traffic, multipart uploads, JWT sessions, and the lifecycle of the scratch directory for each session.
- **Python owns compute.** Each job is a short-lived subprocess: CLI args in, JSON on stdout. No web framework overhead, no long-lived Python state.
- **Single domain.** `next.config.mjs` rewrites `/api/*` to the Fly.io orchestrator URL at the edge, so the browser never sees a second origin and there are no CORS dances.

## Repo layout

```
platform/
├── apps/
│   ├── web/              Next.js 14 app (docs + dashboard)
│   └── orchestrator/     Fastify server that spawns Python
├── compute/              Python package (jobs + CLI)
├── Dockerfile            Node + Python image for orchestrator+compute
├── fly.toml              Fly.io machine config (with volume for drives)
└── .github/workflows/    Deploy pipelines
```

## Local development

Prereqs: Node 20+, pnpm 9+, Python 3.11+.

```bash
cd platform
pnpm install
pip install -e ./compute

# Terminal 1
pnpm dev:orchestrator      # http://localhost:8080

# Terminal 2
pnpm dev:web               # http://localhost:3000 (proxies /api → :8080)
```

Open `http://localhost:3000`. The docs live at `/docs`; the playground lives at `/app/playground/venn`.

## Deployment

Every push to `main` triggers one or both workflows based on the paths changed:

| Change                                                    | Workflow                       | Target  |
| --------------------------------------------------------- | ------------------------------ | ------- |
| `apps/web/**`                                             | `deploy-web.yml`               | Vercel  |
| `apps/orchestrator/**`, `compute/**`, `Dockerfile`, `fly.toml` | `deploy-orchestrator.yml` | Fly.io  |

### One-time deploy setup

**Vercel**
1. Create a Vercel project pointing at this repo; set **Root Directory** to `platform/apps/web`.
2. Add env var `NEXT_PUBLIC_API_ORIGIN` = `https://<your-fly-app>.fly.dev`.
3. Capture `VERCEL_TOKEN`, `VERCEL_ORG_ID`, `VERCEL_PROJECT_ID` into GitHub repo secrets.

**Fly.io**
1. `flyctl launch --no-deploy` from `platform/` (pick a name, keep the provided `fly.toml`).
2. `flyctl volumes create drives --size 3` — persistent "Little Drive" storage.
3. `flyctl secrets set SESSION_JWT_SECRET=$(openssl rand -hex 32)`.
4. Put the Fly deploy token in GitHub secret `FLY_API_TOKEN` (`flyctl tokens create deploy`).

Push to `main` → both services redeploy automatically.

## Adding a new compute job

1. Drop a module in `compute/compute/jobs/<name>.py` exporting `def run(args: argparse.Namespace) -> dict`.
2. Register it in `compute/compute/cli.py`.
3. Add a route in `apps/orchestrator/src/routes/jobs.ts` (or extend the generic handler).
4. Build the dashboard view under `apps/web/app/app/playground/<name>/`.

The contract between Node and Python is deliberately narrow: **file paths in, JSON on stdout.** That makes every job independently testable from the command line.
