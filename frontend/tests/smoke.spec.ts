import { expect, test } from "@playwright/test";

test("creates a project and groups the workspace tabs", async ({ page }) => {
  await page.goto("/");
  await expect(page).toHaveTitle("PepSeMan");
  await expect(page.locator("header")).toHaveCount(0);
  await expect(page.locator(".tabs .app-brand")).toHaveText("PepSeMan");
  await expect(page.locator(".tabs .app-brand img")).toBeVisible();
  await expect(page.getByRole("navigation", { name: "Workspace sections" }).getByLabel("Active project")).toBeVisible();
  await page.getByLabel("Project filename").fill(`Playwright ${Date.now()}.csv`);
  const qcRequest = page.waitForRequest((request) => /\/api\/projects\/\d+\/qc$/.test(new URL(request.url()).pathname) && new URL(request.url()).searchParams.get("page_size") === "5");
  await page.getByRole("button", { name: "Create project" }).click();
  await qcRequest;
  await expect(page.getByLabel("Active project")).toHaveValue(/\d+/);
  const tabBar = await page.getByRole("navigation", { name: "Workspace sections" }).boundingBox();
  const projectSelect = await page.getByLabel("Active project").boundingBox();
  expect(tabBar && projectSelect).toBeTruthy();
  expect(tabBar!.x + tabBar!.width - projectSelect!.x - projectSelect!.width).toBeLessThan(45);
  for (const tab of ["Project", "Analysis", "Licenses"]) {
    await expect(page.getByRole("button", { name: tab, exact: true })).toBeVisible();
  }
  await expect(page.getByRole("button", { name: "Export", exact: true })).toHaveCount(0);
  await expect(page.getByRole("heading", { name: "Sequence dictionary" })).toBeVisible();
  await page.getByRole("button", { name: "Analysis" }).click();
  await expect(page.getByRole("heading", { name: "R Decomposition" })).toBeVisible();
  await expect(page.getByRole("tab", { name: "MMP Table" })).toBeVisible();
  await expect(page.locator(".analysis-upper .panel-heading")).toHaveCount(0);
  await expect(page.getByRole("tab", { name: "Scatter Plot" })).toBeVisible();
  const left = await page.locator(".analysis-left").boundingBox();
  const right = await page.locator(".analysis-right").boundingBox();
  const upper = await page.locator(".analysis-upper").boundingBox();
  const lower = await page.locator(".analysis-lower").boundingBox();
  expect(left && right && upper && lower).toBeTruthy();
  expect(left!.width / (left!.width + right!.width)).toBeGreaterThan(0.23);
  expect(left!.width / (left!.width + right!.width)).toBeLessThan(0.30);
  expect(lower!.y).toBeLessThan(upper!.y);
  expect(upper!.height / (upper!.height + lower!.height)).toBeGreaterThan(0.62);
  expect(upper!.height / (upper!.height + lower!.height)).toBeLessThan(0.70);
  await page.getByRole("button", { name: "Licenses" }).click();
  await expect(page.getByRole("heading", { name: "Licenses & credits" })).toBeVisible();
  await expect(page.locator(".license-panel > .license-table-scroll tbody tr")).toHaveCount(18);
  await expect(page.locator(".license-panel > .license-table-scroll tbody")).toContainText("CC BY-SA 3.0");
  await expect(page.locator(".license-panel > .license-table-scroll tbody")).toContainText("BlueOak-1.0.0");
  await expect(page.locator(".license-panel > .license-table-scroll tbody")).toContainText("caniuse-lite");
  await expect(page.locator(".license-panel a[href='https://github.com/rdkit/rdkit']")).toBeVisible();
  await expect(page.getByRole("heading", { name: "Test data sources" })).toBeVisible();
  await expect(page.locator(".license-sources tbody tr")).toHaveCount(8);
  await expect(page.locator(".license-sources a[href='https://www.ebi.ac.uk/chembl/explore/document/CHEMBL3751807']")).toBeVisible();
  await expect(page.locator(".license-sources")).toContainText("original selection filter not recorded");
});

test("QC ID search resets pagination and page jump opens the requested page", async ({ page }) => {
  const records = Array.from({ length: 12 }, (_, index) => ({ record_id: index + 1, source_molecule_id: `CHEMBL${index + 1}`, status: "accepted", exclusion_reason: null, original_smiles: "C", structure_hash: `hash-${index + 1}` }));
  await page.route("**/api/**", async (route) => {
    const url = new URL(route.request().url());
    if (url.pathname === "/api/projects") { await route.fulfill({ json: [{ project_id: 1, name: "20261028_Fixture", created_at: "2026-10-28", updated_at: "2026-10-28" }] }); return; }
    if (url.pathname === "/api/sequence-dictionary") { await route.fulfill({ json: [] }); return; }
    if (url.pathname.endsWith("/qc")) {
      const filtered = records.filter((row) => row.source_molecule_id.toLowerCase().includes((url.searchParams.get("id_filter") ?? "").toLowerCase()));
      const currentPage = Number(url.searchParams.get("page") ?? "1");
      const pageSize = Number(url.searchParams.get("page_size") ?? "5");
      await route.fulfill({ json: { rows: filtered.slice((currentPage - 1) * pageSize, currentPage * pageSize), total: filtered.length, page: currentPage, page_size: pageSize } });
      return;
    }
    await route.fulfill({ json: { rows: [], total: 0, page: 1, page_size: 50 } });
  });
  await page.goto("/?project=1&tab=Project");
  const qcTable = page.locator(".qc-toolbar + .table-shell");
  await expect(qcTable.locator("tbody tr")).toHaveCount(5);
  await page.getByLabel("Jump to page").fill("3");
  await page.getByRole("button", { name: "Go", exact: true }).click();
  await expect(qcTable).toContainText("CHEMBL11");
  await page.getByLabel("Search QC ID").fill("CHEMBL12");
  await expect(qcTable.locator("tbody tr")).toHaveCount(1);
  await expect(qcTable).toContainText("CHEMBL12");
  await expect(qcTable.locator(".pagination")).toContainText("1 / 1");
});

test("selected sequence filters table and scatter, then opens its pair detail", async ({ page }) => {
  const reference = "reference-hash";
  const residue = (type: string, r1: string) => ({ residue_type: type, r1, r2: null, r3: null, bridge: null, stereo_status: "L" });
  const original = { structure_hash: reference, n_cap: "H", c_cap: "OH", molecular_weight: 500, canonical_isomeric_smiles: "original-smiles", residues: Array.from({ length: 5 }, () => residue("amino acid", "C")), svg: "<svg xmlns='http://www.w3.org/2000/svg'></svg>", bb_svg: "<svg xmlns='http://www.w3.org/2000/svg'></svg>" };
  const converted = { ...original, structure_hash: "paired-hash", molecular_weight: 514, canonical_isomeric_smiles: "converted-smiles", residues: [residue("amino acid", "N"), ...original.residues.slice(1)] };
  const apiRequests: string[] = [];
  let exportPayload: { structure_hashes: string[]; mode: string } | undefined;
  await page.route("**/api/**", async (route) => {
    const url = new URL(route.request().url());
    apiRequests.push(`${url.pathname}${url.search}`);
    if (url.pathname.endsWith("/export")) {
      exportPayload = route.request().postDataJSON();
      await route.fulfill({ status: 200, contentType: "text/csv", body: "structure_hash\nreference-hash\npaired-hash\n" });
      return;
    }
    if (url.pathname.endsWith("/pairs/selection")) {
      await route.fulfill({ json: [1, 2].map((number) => ({ pair_key: `pair-${number}`, first_structure_hash: reference, second_structure_hash: "paired-hash" })) });
      return;
    }
    const body = url.pathname === "/api/projects" ? [{ project_id: 1, name: "Fixture", created_at: "2026-01-01", updated_at: "2026-01-01" }] :
      url.pathname.endsWith("/decompositions") ? { rows: [{ source_molecule_id: "ref", n_pairs: 1, value: 7.2, sequence: "[Ac]AAPAA", structure_hash: reference, status: "accepted" }], total: 1, page: 1, page_size: 1 } :
      /\/pairs\/pair-[12]$/.test(url.pathname) ? { pair_key: url.pathname.split("/").pop(), site_kind: "residue", site_position: 1, before: original, after: converted, scaffold_svg: "<svg xmlns='http://www.w3.org/2000/svg'></svg>", parent_sequence: "[Ac]AAPAA", parent_sequence_highlight_start: 4, parent_sequence_highlight_end: 5, properties: { before_reference: "ref", after_reference: "pair", before_chembl_value: 7.2, after_chembl_value: 7.8, delta_chembl_value: 0.6, before_mw: 500, after_mw: 514, delta_mw: 14, assay_id: "assay", endpoint: "IC50" }, alignment_warning: null } :
      url.pathname.endsWith("/pairs") ? { rows: [1, 2].map((number) => ({ pair_key: `pair-${number}`, first_structure_hash: reference, second_structure_hash: "paired-hash", parent_sequence: "[Ac]AAPAA", parent_sequence_highlight_start: 4, parent_sequence_highlight_end: 5, original: "A", convert: "[4FP1234]", site_kind: "residue", site_position: 1, delta_activity: 0.6, delta_mw: 14, before_mw: 500, after_mw: 514 })), total: 2, page: 1, page_size: 50 } :
      { rows: [], total: 0, page: 1, page_size: 50 };
    await route.fulfill({ json: body });
  });
  await page.goto("/?project=1&tab=Analysis");
  await expect(page.locator(".analysis-left thead th.sticky-column")).toHaveCount(3);
  await expect(page.locator(".analysis-left thead")).toContainText("nPairs");
  await expect(page.locator(".analysis-left thead")).toContainText("Value");
  await expect(page.locator(".analysis-left thead")).toContainText("sequence");
  await expect(page.locator(".analysis-left .panel-heading p:not(.eyebrow)")).toHaveCount(0);
  await expect(page.locator(".analysis-left .pagination").getByRole("button", { name: "Previous", exact: true })).toHaveCount(0);
  await expect(page.getByLabel("Value metric")).toHaveValue("pchembl");
  const mmpHeaders = await page.locator(".analysis-upper thead th").allTextContents();
  expect(mmpHeaders.slice(1, 7).map((item) => item.trim())).toEqual(["Parent sequence", "Original", "Convert", "Site position", "ΔValue", "ΔMW"]);
  await page.locator(".analysis-upper .sort-button").getByText("Convert", { exact: true }).click();
  await expect.poll(() => apiRequests.some((request) => new URL(`http://localhost${request}`).searchParams.get("sort") === "convert:asc")).toBeTruthy();
  await expect(page.locator(".analysis-upper tbody tr").first()).not.toContainText("-X-");
  await expect(page.locator(".analysis-upper tbody tr").first().locator("td").nth(1)).toContainText("[Ac]AAPAA");
  await expect(page.locator(".analysis-upper tbody tr").first().locator("td").nth(1).locator("mark")).toHaveText("A");
  const parentWidth = await page.locator(".analysis-upper thead th").nth(1).boundingBox();
  expect(parentWidth!.width).toBeGreaterThanOrEqual(299);
  const resize = await page.getByRole("separator", { name: "Resize parent_sequence column" }).boundingBox();
  await page.mouse.move(resize!.x + resize!.width / 2, resize!.y + resize!.height / 2);
  await page.mouse.down();
  await page.mouse.move(resize!.x + resize!.width / 2 + 80, resize!.y + resize!.height / 2);
  await page.mouse.up();
  await expect.poll(async () => (await page.locator(".analysis-upper thead th").nth(1).boundingBox())!.width).toBeGreaterThan(parentWidth!.width + 60);
  const headerPositions = () => page.locator(".analysis-left thead th").evaluateAll((cells) => cells.map((cell) => cell.getBoundingClientRect().left));
  const beforeScroll = await headerPositions();
  await page.locator(".analysis-left .table-scroll").evaluate((element) => { element.scrollLeft = 180; });
  const afterScroll = await headerPositions();
  expect(afterScroll.slice(0, 3)).toEqual(beforeScroll.slice(0, 3));
  expect(afterScroll[3]).toBeLessThan(beforeScroll[3]);
  await expect(page.locator(".analysis-left tbody tr.selected-row")).toContainText(reference);
  await expect(page.locator(".analysis-upper tbody tr.selected-row")).toContainText("pair-1");
  await page.locator(".analysis-upper").getByRole("button", { name: "Next row" }).click();
  await expect(page.locator(".analysis-upper tbody tr.selected-row")).toContainText("pair-2");
  await page.locator(".analysis-upper").getByRole("button", { name: "Previous row" }).click();
  await expect(page.locator(".analysis-upper tbody tr.selected-row")).toContainText("pair-1");
  await expect(page.getByRole("row").filter({ hasText: "pair-1" })).toBeVisible();
  await expect.poll(() => apiRequests.some((url) => url.includes(`/pairs?reference_hash=${reference}`))).toBe(true);
  await page.getByRole("checkbox", { name: "Select all MMP pairs" }).click();
  await expect(page.getByRole("checkbox", { name: "Select all MMP pairs" })).toBeChecked();
  await expect(page.getByRole("checkbox", { name: "Select pair-1" })).toBeChecked();
  await expect(page.getByRole("checkbox", { name: "Select pair-2" })).toBeChecked();
  await page.getByRole("checkbox", { name: "Select all MMP pairs" }).click();
  await expect(page.getByRole("checkbox", { name: "Select all MMP pairs" })).not.toBeChecked();
  await page.getByRole("checkbox", { name: "Select pair-1" }).check();
  await expect(page.getByRole("button", { name: "Export parent + selected molecules" })).toBeEnabled();
  const download = page.waitForEvent("download");
  await page.getByRole("button", { name: "Export parent + selected molecules" }).click();
  await download;
  expect(exportPayload).toEqual({ structure_hashes: [reference, "paired-hash"], mode: "selected_compounds" });
  await page.getByRole("row").filter({ hasText: "pair-1" }).click();
  await expect(page.getByRole("heading", { name: "Original BB", exact: true })).toBeVisible();
  await expect(page.getByRole("heading", { name: "Convert BB", exact: true })).toBeVisible();
  await expect(page.getByRole("heading", { name: "Scaffold", exact: true })).toBeVisible();
  await expect(page.locator(".scaffold-card .scaffold-sequence")).toContainText("[Ac]AAPAA");
  await expect(page.locator(".scaffold-card .scaffold-sequence mark")).toHaveText("A");
  await expect(page.locator(".scaffold-direction")).toHaveCount(0);
  await expect(page.getByRole("heading", { name: /Pair difference/ })).toBeVisible();
  await expect(page.locator(".analysis-lower svg")).toHaveCount(3);
  await expect(page.locator(".analysis-lower")).toContainText("ΔChEMBL Value");
  await page.getByRole("tab", { name: "Scatter Plot" }).click();
  await expect(page.getByText("2 pairs with numeric values for both axes.", { exact: false })).toBeVisible();
  await page.getByLabel("Value metric").selectOption("length");
  await expect.poll(() => apiRequests.some((url) => url.includes("value_metric=length"))).toBe(true);
  await page.getByLabel("Filter sequence").fill("PAA");
  await expect.poll(() => apiRequests.some((url) => url.includes("sequence_filter=PAA"))).toBe(true);
});

test("global sequence dictionary edits a shared abbreviation", async ({ page }) => {
  const hash = "a".repeat(64);
  let saved: string | undefined;
  await page.route("**/api/**", async (route) => {
    const url = new URL(route.request().url());
    if (url.pathname.endsWith("/svg")) { await route.fulfill({ contentType: "image/svg+xml", body: "<svg xmlns='http://www.w3.org/2000/svg'></svg>" }); return; }
    const entry = { variant_hash: hash, abbreviation: saved ? `[${saved}]` : "[XAAAA]", default_abbreviation: "[XAAAA]", custom_abbreviation: saved ?? null, residue_type: "Pro-like bridge", stereo_status: "L", r1: null, r2: null, r3: null, bridge: "CCC", structure_url: `/api/sequence-dictionary/${hash}/svg` };
    if (route.request().method() === "PUT") { saved = route.request().postDataJSON().abbreviation; await route.fulfill({ json: { ...entry, abbreviation: `[${saved}]`, custom_abbreviation: saved } }); return; }
    await route.fulfill({ json: url.pathname === "/api/sequence-dictionary" ? [entry] : url.pathname === "/api/projects" ? [] : { rows: [], total: 0, page: 1, page_size: 50 } });
  });
  await page.goto("/?tab=Project");
  await expect(page.locator(".dictionary-scroll tbody tr")).toHaveCount(1);
  await expect(page.locator(".dictionary-scroll img")).toBeVisible();
  await page.getByLabel(`Abbreviation for ${hash}`).fill("4FP");
  await page.locator(".dictionary-scroll").getByRole("button", { name: "Save" }).click();
  await expect(page.locator(".dictionary-scroll tbody tr")).toContainText("[4FP]");
  expect(saved).toBe("4FP");
});

test("row arrows cross MMP table page boundaries", async ({ page }) => {
  await page.route("**/api/**", async (route) => {
    const url = new URL(route.request().url());
    if (/\/pairs\/pair-/.test(url.pathname)) {
      await route.fulfill({ status: 404, json: { detail: "Not part of the navigation fixture" } });
      return;
    }
    const body = url.pathname === "/api/projects" ? [{ project_id: 1, name: "Fixture" }] :
      url.pathname.endsWith("/decompositions") ? { rows: [{ source_molecule_id: "ref", n_pairs: 3, value: 7, structure_hash: "reference-hash", status: "accepted", sequence: "AAAAA" }], total: 1, page: 1, page_size: 1 } :
      url.pathname.endsWith("/pairs") ? { rows: (url.searchParams.get("page") === "2" ? [3] : [1, 2]).map((number) => ({ pair_key: `pair-${number}`, first_structure_hash: "reference-hash", second_structure_hash: `partner-${number}`, parent_sequence: "AAAAA", original: "A", convert: "a" })), total: 3, page: Number(url.searchParams.get("page")), page_size: 2 } :
      { error: "detail unavailable in navigation fixture" };
    await route.fulfill({ json: body });
  });
  await page.goto("/?project=1&tab=Analysis");
  const table = page.locator(".analysis-upper");
  await expect(table.locator("tbody tr.selected-row")).toContainText("pair-1");
  await table.getByRole("button", { name: "Next row" }).click();
  await expect(table.locator("tbody tr.selected-row")).toContainText("pair-2");
  await table.getByRole("button", { name: "Next row" }).click();
  await expect(table.locator("tbody tr.selected-row")).toContainText("pair-3");
  await table.getByRole("button", { name: "Previous row" }).click();
  await expect(table.locator("tbody tr.selected-row")).toContainText("pair-2");
});

test("up and down keys navigate R decomposition and MMP rows without hijacking controls", async ({ page }) => {
  await page.route("**/api/**", async (route) => {
    const url = new URL(route.request().url());
    if (/\/pairs\/pair-/.test(url.pathname)) {
      await route.fulfill({ status: 404, json: { detail: "Not part of the navigation fixture" } });
      return;
    }
    const body = url.pathname === "/api/projects" ? [{ project_id: 1, name: "Fixture" }] :
      url.pathname.endsWith("/decompositions") ? { rows: [1, 2].map((number) => ({ source_molecule_id: `ref-${number}`, n_pairs: 3, value: 7, structure_hash: `reference-${number}`, status: "accepted", sequence: `AAAA${number}` })), total: 2, page: 1, page_size: 2 } :
      url.pathname.endsWith("/pairs") ? { rows: (url.searchParams.get("page") === "2" ? [3] : [1, 2]).map((number) => ({ pair_key: `pair-${number}`, first_structure_hash: url.searchParams.get("reference_hash"), second_structure_hash: `partner-${number}`, parent_sequence: "AAAAA", original: "A", convert: "a" })), total: 3, page: Number(url.searchParams.get("page")), page_size: 2 } :
      { rows: [], total: 0, page: 1, page_size: 50 };
    await route.fulfill({ json: body });
  });
  await page.goto("/?project=1&tab=Analysis");
  const decomposition = page.getByRole("region", { name: "R Decomposition rows" });
  const mmp = page.getByRole("region", { name: "MMP Table rows" });
  await expect(decomposition.locator("tr.selected-row")).toContainText("ref-1");
  await decomposition.focus();
  await page.keyboard.press("ArrowDown");
  await expect(decomposition.locator("tr.selected-row")).toContainText("ref-2");
  await page.keyboard.press("ArrowUp");
  await expect(decomposition.locator("tr.selected-row")).toContainText("ref-1");

  await expect(mmp.locator("tr.selected-row")).toContainText("pair-1");
  await mmp.locator("tbody tr").first().click();
  await expect(mmp).toBeFocused();
  await page.keyboard.press("ArrowDown");
  await expect(mmp.locator("tr.selected-row")).toContainText("pair-2");
  await page.keyboard.press("ArrowDown");
  await expect(mmp.locator("tr.selected-row")).toContainText("pair-3");
  await expect(mmp).toBeFocused();
  await page.keyboard.press("ArrowUp");
  await expect(mmp.locator("tr.selected-row")).toContainText("pair-2");

  await page.getByRole("checkbox", { name: "Select pair-2" }).focus();
  await page.keyboard.press("ArrowDown");
  await expect(mmp.locator("tr.selected-row")).toContainText("pair-2");
  await page.getByLabel("Filter sequence").focus();
  await page.keyboard.press("ArrowDown");
  await expect(decomposition.locator("tr.selected-row")).toContainText("ref-1");
});

test("cyclic topology filter, link markers, and display origin stay in UI state", async ({ page }) => {
  const links = [
    { kind: "head_to_tail", first: { kind: "residue", position: 6 }, second: { kind: "residue", position: 1 } },
    { kind: "disulfide", first: { kind: "residue", position: 2 }, second: { kind: "residue", position: 5 } },
  ];
  const residues = Array.from({ length: 6 }, (_, index) => ({ position: index + 1, residue_type: "amino acid", stereo_status: "L", r1: "C" }));
  const structure = { structure_hash: "cyclic-parent", canonical_isomeric_smiles: "C", n_cap: null, c_cap: null, residues, bb_svg: "<svg xmlns='http://www.w3.org/2000/svg'></svg>" };
  const requests: string[] = [];
  await page.route("**/api/**", async (route) => {
    const url = new URL(route.request().url());
    requests.push(`${url.pathname}${url.search}`);
    const body = url.pathname === "/api/projects" ? [{ project_id: 1, name: "Cyclic" }] :
      url.pathname.endsWith("/decompositions") ? { rows: [{ source_molecule_id: "ring", n_pairs: 1, value: 7, sequence: "ACDEFG", structure_hash: "cyclic-parent", status: "accepted", length: 6, backbone_kind: "head_to_tail", cyclization_types: "disulfide, head_to_tail", link_count: 2, cyclization_links: links }], total: 1, page: 1, page_size: 1 } :
      url.pathname.endsWith("/pairs/pair-ring") ? { pair_key: "pair-ring", site_kind: "residue", site_position: 3, partner_site_position: 3, before: structure, after: { ...structure, structure_hash: "cyclic-partner" }, scaffold_svg: "<svg xmlns='http://www.w3.org/2000/svg'></svg>", parent_sequence: "ACDEFG", parent_sequence_highlight_start: 2, parent_sequence_highlight_end: 3, backbone_kind: "head_to_tail", cyclization_links: links, properties: {} } :
      url.pathname.endsWith("/pairs") ? { rows: [{ pair_key: "pair-ring", first_structure_hash: "cyclic-parent", second_structure_hash: "cyclic-partner", parent_sequence: "ACDEFG", original: "D", convert: "V", site_position: 3, cyclization_links: links, cyclization_types: "disulfide, head_to_tail", residue_count: 6 }], total: 1, page: 1, page_size: 50 } :
      { rows: [], total: 0, page: 1, page_size: 50 };
    await route.fulfill({ json: body });
  });
  await page.goto("/?project=1&tab=Analysis");
  await expect(page.locator(".analysis-left .cycle-symbol")).toBeVisible();
  await expect(page.locator(".analysis-left .link-marker")).toHaveCount(2);
  await page.getByLabel("Filter cyclization type").selectOption("disulfide");
  await expect.poll(() => requests.some((request) => request.includes("cyclization_type=disulfide"))).toBe(true);
  await expect(page.getByLabel("Display origin residue")).toBeVisible();
  await page.getByLabel("Display origin residue").selectOption("3");
  await expect(page).toHaveURL(/origin=3/);
  await expect(page.locator(".scaffold-card .cycle-symbol")).toBeVisible();
  await expect(page.locator(".scaffold-card .sequence-residue").first()).toHaveText("D");
  await expect(page.locator(".scaffold-card .link-legend-item")).toHaveCount(2);
});

