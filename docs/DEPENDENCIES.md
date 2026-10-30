# Third-party dependencies and license scope

PepSeMan's own code, documentation, and project-created icons are under the repository [MIT license](../LICENSE). Local source data and the derived dictionary are outside that grant; see [DATA_LICENSE.md](../DATA_LICENSE.md). This page describes packages selected by the manifests, not a grant of rights to their code under PepSeMan's MIT license.

## Direct dependencies

| Dependency | Use | Upstream license |
| --- | --- | --- |
| FastAPI | Local HTTP API | [MIT](https://github.com/fastapi/fastapi/blob/master/LICENSE) |
| Pydantic | API validation | [MIT](https://github.com/pydantic/pydantic/blob/main/LICENSE) |
| RDKit | Molecular parsing and depiction | [BSD-3-Clause](https://github.com/rdkit/rdkit/blob/master/license.txt) |
| Uvicorn | Local server | [BSD-3-Clause](https://github.com/Kludex/uvicorn/blob/main/LICENSE.md) |
| React / React DOM | User interface | [MIT](https://github.com/facebook/react/blob/main/LICENSE) |
| TanStack Table | Data tables | [MIT](https://github.com/TanStack/table/blob/main/LICENSE) |
| Plotly.js / react-plotly.js | Scatter plots | [MIT](https://github.com/plotly/plotly.js/blob/main/LICENSE) / [MIT](https://github.com/plotly/react-plotly.js/blob/main/LICENSE) |
| HTTPX | Python test client | [BSD-3-Clause](https://github.com/encode/httpx/blob/master/LICENSE.md) |
| pytest / pytest-cov | Python tests and coverage | [MIT](https://github.com/pytest-dev/pytest/blob/main/LICENSE) / [MIT](https://github.com/pytest-dev/pytest-cov/blob/master/LICENSE) |
| Vite | Frontend build tooling | [MIT](https://github.com/vitejs/vite/blob/main/packages/vite/LICENSE.md) |
| TypeScript | Frontend type checking | [Apache-2.0](https://github.com/microsoft/TypeScript/blob/main/LICENSE.txt) |
| Playwright | Browser tests | [Apache-2.0](https://github.com/microsoft/playwright/blob/main/LICENSE) |

Python versions are pinned in `pyproject.toml` / `requirements.txt`; JavaScript versions, including transitive packages, are recorded in `frontend/package-lock.json`. A fresh installation downloads packages from their respective distributors; `node_modules`, Python environments, and built frontend bundles are **not** in the planned Git source tree.

## JavaScript transitive inventory

The lockfile was inspected on 2026-10-30. Its 355 package entries have these declared license identifiers: MIT 282, ISC 32, BSD-3-Clause 21, BSD-2-Clause 6, Apache-2.0 6, `(MIT OR Apache-2.0)` 1, CC-BY-4.0 1, BlueOak-1.0.0 1, Unlicense 1, Zlib 1, 0BSD 1, and two without a `license` field. No GPL/AGPL or non-commercial license identifier appeared in this lockfile. A metadata identifier is an audit clue, not a substitute for the package's license text.

- `caniuse-lite` (CC-BY-4.0) is a development/build-time dependency through Browserslist. Its upstream [README](https://github.com/browserslist/caniuse-lite/blob/main/README.md#license) requests attribution to caniuse.com if its data are redistributed. Its dataset is not copied into the planned Git source tree.
- `sax` (BlueOak-1.0.0) is an indirect dependency through Plotly.js. The [license](https://github.com/isaacs/sax-js/blob/main/LICENSE.md) requires the license text or a link to accompany copies of its code. This link is recorded here for distributions that include it.
- The two missing lockfile fields are `glsl-resolve/node_modules/xtend@2.2.0` and `stack-trace@0.0.9`. Their installed upstream package contents identify MIT in `xtend/README.md` and the full MIT text in `stack-trace/License`, respectively. Recheck them when updating the lockfile.

For the **planned source repository**, no third-party package files or prebuilt JavaScript are redistributed. If a future release publishes `frontend/dist`, a desktop binary, a container, or vendored dependencies, redo an artifact-level inventory and include the full applicable copyright/license/NOTICE materials in that distribution. A link table alone may not satisfy MIT, BSD, or Apache redistribution notice conditions. Likewise, review newly added fonts, photographs, figures, datasets, and contributor-provided code before including them.

