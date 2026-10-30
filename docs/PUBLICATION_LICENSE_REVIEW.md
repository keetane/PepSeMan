# Publication license review

Reviewed 2026-10-30 for a **source-code-only GitHub repository**. GitHub Pages publication was canceled. This is an engineering inventory, not a legal opinion or a guarantee about material outside the reviewed tree. ChEMBL-derived source data and its history are handled separately in [DATA_LICENSE.md](../DATA_LICENSE.md).

## What is actually published in this scope

| Material | Disposition |
| --- | --- |
| Application source, documentation, launcher, and project icon variants | Covered by the project's [MIT license](../LICENSE), with the user-selected `keetane` copyright attribution. The 8-bead motif was requested for this project; the PNG, favicon, and ICNS are variants of the project icon. No external font, photograph, paper figure, or downloaded icon is included. |
| Local test data and generated dictionary | Ignored and absent from this new public repository, which starts from a source-only snapshot. Source record links, not CSV rows, appear in the app's Licenses tab. Earlier copies outside this repository may remain in forks, clones, or caches. |
| Python and JavaScript dependencies | Referenced by pinned manifests/lockfile, not vendored or included in the Git source tree. Their upstream terms and the transitive JavaScript license categories are in [DEPENDENCIES.md](DEPENDENCIES.md). A user installing or building the app obtains those packages under their own licenses, not PepSeMan's MIT license. |
| External literature and web sources | Linked for background/provenance; no article figure, full text, or web image is included in the reviewed publication files. |

The direct dependencies have MIT, BSD-3-Clause, or Apache-2.0 terms. The 355 JavaScript lockfile entries include one CC-BY-4.0 *build-time dataset* (`caniuse-lite`), one BlueOak-1.0.0 *indirect package* (`sax`), and two missing metadata fields that were resolved from their installed package documents as MIT. No GPL/AGPL or non-commercial identifier was found in this lockfile. The package terms do not impose a copyleft license on this source repository merely because it references them; redistribution of package code or a built bundle would require a separate artifact-level notice review. See the upstream [MIT](https://opensource.org/license/mit), [BSD-3-Clause](https://opensource.org/license/bsd-3-clause), [Apache-2.0](https://www.apache.org/licenses/LICENSE-2.0), [`caniuse-lite`](https://github.com/browserslist/caniuse-lite/blob/main/README.md#license), and [`sax`](https://github.com/isaacs/sax-js/blob/main/LICENSE.md) terms.

## Release boundary

No non-ChEMBL license blocker was identified for the **reviewed source tree**. This finding does not cover a future desktop/container release, a published `frontend/dist` bundle, newly added assets or contributors, or unrecorded third-party provenance. Recheck those and include any required license/copyright/NOTICE files before extending the distribution. This new public repository contains no prior dataset commits; the earlier repository remains a private archive, and copies held elsewhere are not retracted. The license review itself is not legal advice.
