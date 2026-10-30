# PepSeMan

PepSeMan is a local-first prototype for exploring residue-level matched molecular pairs (MMPs) in linear and selected **cyclic peptides**. It imports molecular structures and assay rows, identifies unambiguous peptide backbones and bridges, and compares molecules that differ at exactly one residue or terminal cap within the same topology.

This is research software, not a validated prediction or clinical decision tool. It runs on your computer; no ChEMBL API or other network service is required at runtime.

## What it does

- Imports CSV, SMI, and SDF while retaining original rows, SMILES, measurement provenance, and QC exclusion reasons.
- Uses conservative RDKit backbone and bridge parsing for peptides of at least five residues. It supports head-to-tail closure, disulfides, thioethers (including cap–side-chain bridges), ring C=C bridges, 1,2,3-triazole bridges, and multiple mixed bridges. Ambiguous structures are marked unsupported rather than assigned a guessed sequence.
- Finds strict one-position residue/cap MMPs using indexed sequence contexts and a topology signature. Ring opening/closing, bridge relocation, and changes to bridge chemistry are not MMPs. Project data is isolated; residue decompositions and sequence abbreviations can be reused across projects.
- Displays natural L residues as uppercase letters, natural D residues as lowercase letters, acetyl caps as `[Ac]`, and other building blocks as editable bracketed abbreviations. The shared dictionary includes representative residue structure drawings.
- Lets you search QC records by ID, filter R decomposition by sequence and cyclization type, inspect matched pairs and activity/weight differences, and export the parent plus selected molecules as CSV with normalized sequence and bridge information.

## Requirements and quick start

- Python **3.12** and Node.js **22+** with npm.
- An RDKit-compatible platform. The prototype has been tested on macOS Apple Silicon; other platforms are not yet verified.

Run these commands from the repository root:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cd frontend
npm ci
npm run build
cd ..
```

If the private test files are present locally, run `python -m scripts.seed_prototype` to load `data/source/20261030_MDM2.csv` and `data/source/20261030_GLP-1.csv` into separate projects. These files are **not included in a fresh Git checkout**; pass your own permitted CSVs with `--file`. The seed skips an identical completed import on subsequent runs. Start the API and persistent job worker in **separate terminals**, with the virtual environment activated in each:

To create the combined cyclic test project from the local literature CSVs in `data/source/`, run `python -m scripts.seed_ring_peptide`. It stages one `ring-Peptide` snapshot and preserves the originating filename in each row's `Source File` field; the worker processes the queued job. The optional local `RAP1.csv` adds triazole examples when present.

```bash
python -m app
```

```bash
python -m app.worker
```

Then open <http://127.0.0.1:8000>. The built frontend is served by the API. The worker is needed for new imports and background analysis; closing the browser does not stop queued jobs.

On macOS, after setting up the environment and building the frontend, double-click `Launch PepSeMan.command` in the repository root. It runs the API and worker on `127.0.0.1:48001`, writes logs and PID files under `data/`, and opens the browser. The launcher uses the included `PepSeMan.icns` for its Finder icon (applied on first launch when Swift is available). It checks the repository's `.venv` first. If the port is already used by a server pointing at another data directory, it stops with a warning instead of showing the wrong project data.

For frontend development, run `cd frontend && npm run dev` separately; Vite serves <http://127.0.0.1:5173> and proxies `/api` to the API on port 8000.

## Interface

| Section | Purpose |
| --- | --- |
| **Project** | Create/select projects, import data, inspect QC records and job history, and edit the cross-project sequence dictionary. QC shows five rows per page with ID search and page jump. |
| **Analysis** | Choose a reference in R Decomposition; inspect its MMP table or scatter plot; select a pair to see the parent scaffold, changed building blocks, and property differences. Cyclic sequences show numbered bridge endpoints and head-to-tail display-origin controls. |
| **Licenses** | View separate code/data terms and credits for direct third-party dependencies. |

Project names normally follow `yyyymmdd_filename`; the combined cyclic test project is named `ring-Peptide`. In the MMP table, the changed building block is highlighted in yellow within the original parent sequence; the Original and Convert columns show the corresponding abbreviations. Linear scaffolds run N→C; head-to-tail structures are drawn as rings. Bridge bonds are color-coded and the changed residue remains yellow. Changing the head-to-tail display origin only rotates the visible sequence; it does not change stored identity or MMP results. The checkbox in the table header selects all pairs for the current parent, including pairs on other pages. CSV export includes the sequences, cyclization links, and active measurement data for the parent and checked partners.

## Prototype data and local storage

Source CSVs and the derived sequence-dictionary snapshot are local-only and Git-ignored. They are ChEMBL-derived data under separate [data terms](DATA_LICENSE.md), not the software's MIT license. `data/source/RAP1.csv` is an optional 25-row local triazole test file. Structure drawings are regenerated from the source molecules. GLP-1 has structures but no assay/pChEMBL columns, so activity values and activity deltas are unavailable for that project.

Runtime SQLite data, all source CSVs, uploads, generated assets, logs, and build/test artifacts are excluded by [`.gitignore`](.gitignore). The `ring-Peptide` project combines the local cyclic source files in one import snapshot; `oxitosin_SS` is byte-identical to `Oxytocin_S` and is not double-counted. Each combined row retains its original filename in `Source File`. Set `PMMP_DATA_DIR` or `PMMP_DB` to relocate runtime storage. After intentionally editing dictionary aliases, run `python -m scripts.dictionary_template` to refresh the local snapshot; use `python -m scripts.dictionary_template --check` to compare it with the current local database.

## Validation and limitations

```bash
python -m pytest
cd frontend
npm run build
npm run test:e2e
```

Browser tests require a Playwright browser installation. `npm audit --audit-level=high` can be run to check current JavaScript dependency advisories. See the [validation report](docs/VALIDATION_REPORT.md) for measured MDM2 prototype results and [architecture](docs/ARCHITECTURE.md) for the data flow.

Only exact `=` measurements with the same endpoint and assay produce a numeric activity delta. Inequality measurements remain visible but are not converted into exact values. Ambiguous/branched backbones, multiple fragments, residue insertions/deletions, and multi-site transformations remain outside the supported analysis scope. Cyclic assignment is conservative: unusual bridges and non-unique residue or bridge interpretations go to QC. The “olefin metathesis” category means a **ring C=C bridge inferred from structure**, not proof of the synthetic route; hydrogenated bridges are not automatically assigned this category. The synthetic cyclic fixtures are automated regression tests, not a measured classification-accuracy claim. Real cyclic data still require human-reviewed labels before accuracy can be reported. The parser may reject a structure that a human can interpret; a rejected or missing result must not be treated as evidence of absence.

Parser, MMP, and renderer versions are tracked separately. Existing linear revisions remain available; on database upgrade, projects with an accessible completed import are queued for reanalysis and the new result is activated in one database transaction after processing succeeds. Keep the worker running until these jobs complete. The shared dictionary still contains only editable bracketed non-natural residue aliases; it does not treat a disulfide-bound cysteine as a new unknown residue.

The API is intended for **loopback-only** use and does not implement authentication. Do not bind it to a public interface or place it behind an internet-facing proxy without adding appropriate access controls.

## License and data attribution

PepSeMan source code, accompanying documentation, and project-created icons and favicon are licensed under the [MIT License](LICENSE), with copyright attributed to keetane. This includes the visual assets in `assets/`, `frontend/public/`, and `PepSeMan.icns`. Local ChEMBL-derived CSVs and the derived sequence-dictionary snapshot are **not MIT-licensed** and are not tracked in this source-only release; their provenance and attribution limitations are documented in [DATA_LICENSE.md](DATA_LICENSE.md). Third-party software remains under its respective licenses; see the [publication license review](docs/PUBLICATION_LICENSE_REVIEW.md), [dependency credits](docs/DEPENDENCIES.md), and the in-app **Licenses** tab. The Licenses tab links to the ChEMBL records behind the local test datasets without bundling their CSV rows.

The original ChEMBL release and download details for the local CSVs were not recorded. ChEMBL also cautions that some calculated-property columns derive from commercial software. Review [the data notice](DATA_LICENSE.md) before republishing or reusing those fields. This public repository starts from a source-only snapshot; the prior repository was moved to a private archive. Copies held in forks, clones, or caches outside this repository may still remain. These setup instructions do not perform a GitHub push or release.

