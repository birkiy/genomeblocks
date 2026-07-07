---
title: About this site
layout: default
nav_exclude: true
search_exclude: true
---

# About the documentation site

This `docs/` folder is a Jekyll site built with the [`just-the-docs`](https://just-the-docs.com) remote theme and deployed to GitHub Pages via the `.github/workflows/pages.yml` Actions workflow.

## Building locally

```bash
cd docs
bundle install
bundle exec jekyll serve --livereload
# open http://127.0.0.1:4000
```

## Enabling Pages on the repository

1. In the repo settings → **Pages**, set the source to **GitHub Actions**.
2. Push to `main`; the `Deploy GitHub Pages` workflow runs automatically on every change under `docs/`.
3. The site goes live at `https://<user>.github.io/<repo>/`.

## File layout

```
docs/
├── _config.yml          # Jekyll + just-the-docs config
├── Gemfile              # local-build dependencies
├── index.md             # landing page
├── installation.md
├── quickstart.md
├── concepts.md
├── release-notes.md
├── credits.md
├── guide/
│   ├── index.md         # user guide hub
│   ├── loci.md
│   ├── genes.md
│   ├── architecture.md
│   ├── signal.md
│   ├── browser.md
│   ├── bedpe.md
│   ├── motifs.md
│   └── atlas.md
├── api/
│   ├── index.md
│   ├── locus.md
│   ├── loci.md
│   ├── genes.md
│   ├── architecture.md
│   ├── signal.md
│   ├── browser.md      # documents the browserview module
│   ├── bedpe.md
│   ├── motifs.md
│   └── atlas.md
├── walkthrough/         # AR & FOXA1 real-data example
├── tools/               # in-browser Venn tool
└── assets/
```

Each page uses `just-the-docs` front matter:

- `nav_order`: position in the left sidebar.
- `parent`: nest inside a parent section.
- `has_children: true`: make the page a section root.
