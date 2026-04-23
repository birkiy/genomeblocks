# syntax=docker/dockerfile:1.6
# -----------------------------------------------------------------------------
# Orchestrator + Compute image
# -----------------------------------------------------------------------------
# Ships Node (Fastify) and Python (compute jobs) together. Node spawns Python
# per-request via child_process; the two halves share the container filesystem
# so the per-session "Little Drive" is visible to both.
# -----------------------------------------------------------------------------

FROM python:3.12-slim AS base
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    NODE_ENV=production

# Install Node 20 (from NodeSource) + build essentials we'll need for any
# native deps. Keep this layer minimal — only what both stages rely on.
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates build-essential \
    && curl -fsSL https://deb.nodesource.com/setup_20.x | bash - \
    && apt-get install -y --no-install-recommends nodejs \
    && rm -rf /var/lib/apt/lists/* \
    && npm install -g pnpm@9.12.0

WORKDIR /app

# -----------------------------------------------------------------------------
# Python deps (compute) — installed first so TS changes don't bust this layer
# -----------------------------------------------------------------------------
COPY compute/pyproject.toml compute/
RUN pip install --upgrade pip && pip install -e ./compute
COPY compute/ compute/

# -----------------------------------------------------------------------------
# Node deps + build (orchestrator)
# -----------------------------------------------------------------------------
COPY pnpm-workspace.yaml package.json ./
COPY apps/orchestrator/package.json apps/orchestrator/
RUN pnpm install --filter orchestrator --frozen-lockfile=false

COPY apps/orchestrator/ apps/orchestrator/
RUN pnpm --filter orchestrator run build

# -----------------------------------------------------------------------------
# Runtime
# -----------------------------------------------------------------------------
ENV PORT=8080 \
    DRIVES_ROOT=/data/drives
EXPOSE 8080

# Fly.io mounts the volume here — orchestrator writes session dirs underneath.
VOLUME ["/data"]

CMD ["node", "apps/orchestrator/dist/index.js"]
