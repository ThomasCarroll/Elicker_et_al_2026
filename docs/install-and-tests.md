---
title: Install And Tests
nav_order: 3
---

# Install And Tests

## Installation

The workflow is designed to run through Docker.

### Build the image

```bash
docker build -t fly-pipeline .
```

### Local repository contents

The repository tracks workflow code, configuration, and tests. Large reference files and generated outputs are intentionally excluded from Git and should remain local.

## Running the workflow directly

```bash
docker run --rm -it -v "$PWD:/work" -w /work fly-pipeline snakemake -j 4
```

## Smoke test

The smoke test checks that:

- the DAG builds
- early workflow stages run
- expected output files are created

Run it with:

```bash
bash tests/smoke_test.sh fly-pipeline
```

## Full test

The full test runs the whole configured workflow and verifies the expected final outputs.

Run it with:

```bash
bash tests/full_test.sh fly-pipeline
```

To measure a clean end-to-end run:

```bash
rm -rf results
time bash tests/full_test.sh fly-pipeline
```

## Updating the Pages site

This documentation site is intentionally simple.

- Each page is a Markdown file under `docs/`
- Navigation is controlled by `docs/_config.yml`
- To add a new page, create a new `.md` file in `docs/` and list it in `header_pages`

## Publishing with GitHub Pages

After pushing the repository:

1. Open the repository on GitHub.
2. Go to `Settings`.
3. Open `Pages`.
4. Set the source to `Deploy from a branch`.
5. Choose the `main` branch and the `/docs` folder.
6. Save.

GitHub will then publish the files from `docs/` as a Pages site.
