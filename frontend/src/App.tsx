import { FormEvent, useCallback, useEffect, useRef, useState } from "react";
import Plot from "react-plotly.js";
import { api, Page, Row } from "./api";
import { DataTable } from "./DataTable";
import { HighlightedSequence } from "./HighlightedSequence";
import { CyclicSequence } from "./CyclicSequence";

const tabs = ["Project", "Analysis", "Licenses"] as const;
const pairLeadingColumns = ["parent_sequence", "original", "convert", "site_position", "delta_activity", "delta_mw"];
type Tab = (typeof tabs)[number];
type AnalysisView = "MMP Table" | "Scatter Plot";
type Project = { project_id: number; name: string; created_at: string; updated_at: string };
type DictionaryEntry = { variant_hash: string; abbreviation: string; default_abbreviation: string; custom_abbreviation: string | null; residue_type: string; stereo_status: string; r1: string | null; r2: string | null; r3: string | null; bridge: string | null; structure_url: string };

function readUrl() {
  const value = new URLSearchParams(location.search);
  return {
    projectId: Number(value.get("project")) || 0,
    tab: (tabs.includes(value.get("tab") as Tab) ? value.get("tab") :
      ["R Decomposition", "MMP Table", "Scatter Plot", "Pair Detail", "Export"].includes(value.get("tab") || "") ? "Analysis" : "Project") as Tab,
    pair: value.get("pair") || "",
    reference: value.get("reference") || "",
    origin: Math.max(1, Number(value.get("origin")) || 1),
  };
}

export default function App() {
  const initial = readUrl();
  const [projects, setProjects] = useState<Project[]>([]);
  const [projectId, setProjectId] = useState(initial.projectId);
  const [tab, setTab] = useState<Tab>(initial.tab);
  const [pairKey, setPairKey] = useState(initial.pair);
  const [refresh, setRefresh] = useState(0);
  const [reference, setReference] = useState(initial.reference);
  const [analysisView, setAnalysisView] = useState<AnalysisView>("MMP Table");
  const [selectedPairs, setSelectedPairs] = useState<Map<string, string>>(new Map());
  const referenceRef = useRef(reference);
  referenceRef.current = reference;
  const [message, setMessage] = useState("");
  const [valueMetric, setValueMetric] = useState("pchembl");
  const [sequenceSearch, setSequenceSearch] = useState("");
  const [cyclizationFilter, setCyclizationFilter] = useState("all");
  const [displayOrigin, setDisplayOrigin] = useState(initial.origin);

  const loadProjects = useCallback(() => {
    api<Project[]>("/api/projects").then((items) => {
      setProjects(items);
      if (!projectId && items.length) setProjectId(items[0].project_id);
    });
  }, [projectId]);
  useEffect(loadProjects, [loadProjects]);

  useEffect(() => {
    const params = new URLSearchParams();
    if (projectId) params.set("project", String(projectId));
    params.set("tab", tab);
    if (reference) params.set("reference", reference);
    if (pairKey) params.set("pair", pairKey);
    if (displayOrigin > 1) params.set("origin", String(displayOrigin));
    history.replaceState(null, "", `?${params}`);
  }, [projectId, tab, pairKey, reference, displayOrigin]);

  const choosePair = (row: Row) => {
    setPairKey(String(row.pair_key));
  };
  const updatePairSelection = (keys: Set<string>, rows: Row[]) => {
    setSelectedPairs((current) => {
      const next = new Map([...current].filter(([key]) => keys.has(key)));
      for (const row of rows) if (keys.has(String(row.pair_key))) {
        next.set(String(row.pair_key), String(row.first_structure_hash === reference ? row.second_structure_hash : row.first_structure_hash));
      }
      return next;
    });
  };
  const selectAllPairs = async (checked: boolean) => {
    if (!checked) { setSelectedPairs(new Map()); return; }
    if (!reference) return;
    const selectedReference = reference;
    try {
      const rows = await api<Row[]>(`/api/projects/${projectId}/pairs/selection?reference_hash=${encodeURIComponent(selectedReference)}`);
      if (referenceRef.current !== selectedReference) return;
      setSelectedPairs(new Map(rows.map((row) => [String(row.pair_key), String(row.first_structure_hash === selectedReference ? row.second_structure_hash : row.first_structure_hash)])));
    } catch { setMessage("Could not select all MMP pairs."); }
  };
  const exportSelectedPairs = async () => {
    if (!reference || !selectedPairs.size) return;
    const structureHashes = [...new Set([reference, ...selectedPairs.values()])];
    const response = await fetch(`/api/projects/${projectId}/export`, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ structure_hashes: structureHashes, mode: "selected_compounds" }) });
    if (!response.ok) { setMessage("Selected compounds could not be exported."); return; }
    const href = URL.createObjectURL(await response.blob());
    const anchor = document.createElement("a"); anchor.href = href; anchor.download = `project-${projectId}-selected-compounds.csv`; anchor.click(); URL.revokeObjectURL(href);
  };
  const selectedProject = projects.find((item) => item.project_id === projectId);

  return (
    <div className="app">
      <nav className="tabs" aria-label="Workspace sections">
        <div className="app-brand" aria-label="PepSeMan"><img src="/app-icon.png" alt="" /><span>PepSeMan</span></div>
        {tabs.map((item) => (
          <button type="button" key={item} className={tab === item ? "active" : ""} onClick={() => setTab(item)}>
            {item}
          </button>
        ))}
        <div className="project-switcher">
          <label htmlFor="project">Active project</label>
          <select id="project" value={projectId} onChange={(event) => { setProjectId(Number(event.target.value)); setReference(""); setPairKey(""); setSelectedPairs(new Map()); }}>
            <option value={0}>Select a project</option>
            {projects.map((project) => <option value={project.project_id} key={project.project_id}>{project.name}</option>)}
          </select>
        </div>
      </nav>
      <main>
        <div className="context-bar">
          <span className="status-dot" />
          {selectedProject ? `${selectedProject.name} · project ${selectedProject.project_id}` : "No project selected"}
          <span className="context-note">5+ residues · strict one-site MMP · N→C</span>
        </div>
        {message && <div className="notice">{message}</div>}
        {tab === "Project" && <><Projects projects={projects} reload={loadProjects} select={(id) => { setProjectId(id); setReference(""); setPairKey(""); setSelectedPairs(new Map()); }} /><SequenceDictionary /><ImportQc projectId={projectId} refresh={refresh} setRefresh={setRefresh} setMessage={setMessage} /><JobHistory projectId={projectId} refresh={refresh} setRefresh={setRefresh} setMessage={setMessage} /></>}
        {tab === "Analysis" && <div className="analysis-workspace">
          <Workspace className="analysis-left" title="R Decomposition" actions={<label className="value-selector"><select aria-label="Value metric" value={valueMetric} onChange={(event) => setValueMetric(event.target.value)}><option value="pchembl">pChEMBL Value</option><option value="molecular_weight">Molecular Weight</option><option value="length">Peptide Length</option></select></label>}>
            <div className="sequence-search"><input type="search" aria-label="Filter sequence" placeholder="Filter sequence… (case-sensitive)" value={sequenceSearch} onChange={(event) => { setSequenceSearch(event.target.value); setReference(""); setPairKey(""); setSelectedPairs(new Map()); }} /><select aria-label="Filter cyclization type" value={cyclizationFilter} onChange={(event) => { setCyclizationFilter(event.target.value); setReference(""); setPairKey(""); setSelectedPairs(new Map()); }}><option value="all">All topologies</option><option value="linear">Linear</option><option value="head_to_tail">Head-to-tail</option><option value="thioether">Thioether</option><option value="disulfide">Disulfide</option><option value="olefin_metathesis">Ring C=C (structure-inferred)</option><option value="triazole">Triazole bridge</option></select></div>
            <DataTable endpoint={`/api/projects/${projectId}/decompositions?value_metric=${valueMetric}&sequence_filter=${encodeURIComponent(sequenceSearch)}&cyclization_type=${cyclizationFilter}`} refreshKey={refresh} onRowClick={(row) => { if (row.status === "accepted") { setReference(String(row.structure_hash)); setPairKey(""); setDisplayOrigin(1); setSelectedPairs(new Map()); } }} activeRowKey="structure_hash" activeRowValue={reference} navigationLabel="R Decomposition rows" stickyColumns={["source_molecule_id", "n_pairs", "value"]} autoSelectFirst unpaginated />
          </Workspace>
          <div className="analysis-right">
            {reference && pairKey ? <PairDetail projectId={projectId} pairKey={pairKey} reference={reference} displayOrigin={displayOrigin} setDisplayOrigin={setDisplayOrigin} /> :
              <section className="panel analysis-lower"><div className="detail-topline"><h2>Pair detail</h2><span>Scaffold · Original BB · Convert BB · Properties</span></div><div className="empty">Select a pair from the table or scatter plot below.</div></section>}
            <section className="panel analysis-upper">
              <div className="subtabs" role="tablist" aria-label="MMP views"><button role="tab" aria-selected={analysisView === "MMP Table"} onClick={() => setAnalysisView("MMP Table")}>MMP Table</button><button role="tab" aria-selected={analysisView === "Scatter Plot"} onClick={() => setAnalysisView("Scatter Plot")}>Scatter Plot</button><div className="mmp-export"><span>{selectedPairs.size} selected</span><button type="button" disabled={!reference || !selectedPairs.size} onClick={exportSelectedPairs}>Export parent + selected molecules</button></div></div>
              {!reference ? <div className="empty">Select an accepted reference sequence on the left.</div> : analysisView === "MMP Table" ?
                <DataTable key={reference} endpoint={`/api/projects/${projectId}/pairs?reference_hash=${encodeURIComponent(reference)}`} refreshKey={refresh} onRowClick={choosePair} activeRowKey="pair_key" activeRowValue={pairKey} navigationLabel="MMP Table rows" selectionKey="pair_key" selectedKeys={new Set(selectedPairs.keys())} onSelectionChange={updatePairSelection} onSelectAll={selectAllPairs} leadingColumns={pairLeadingColumns} autoSelectFirst resizableColumns /> :
                <Scatter projectId={projectId} reference={reference} choosePair={choosePair} refresh={refresh} />}
            </section>
          </div>
        </div>}
        {tab === "Licenses" && <LicensesPanel />}
      </main>
    </div>
  );
}

const licenses = [
  { name: "PepSeMan", role: "Application code, docs, and icons", license: "MIT", credit: "keetane", source: "https://github.com/keetane/PepSeMan", terms: "https://github.com/keetane/PepSeMan/blob/main/LICENSE" },
  { name: "ChEMBL", role: "Locally supplied prototype data, when present", license: "CC BY-SA 3.0", credit: "EMBL-EBI ChEMBL", source: "https://www.ebi.ac.uk/chembl/", terms: "https://creativecommons.org/licenses/by-sa/3.0/" },
  { name: "RDKit", role: "Molecular parsing and depiction", license: "BSD-3-Clause", credit: "RDKit contributors", source: "https://github.com/rdkit/rdkit", terms: "https://github.com/rdkit/rdkit/blob/master/license.txt" },
  { name: "FastAPI", role: "Python API framework", license: "MIT", credit: "FastAPI contributors", source: "https://github.com/fastapi/fastapi", terms: "https://github.com/fastapi/fastapi/blob/master/LICENSE" },
  { name: "Pydantic", role: "API data validation", license: "MIT", credit: "Pydantic contributors", source: "https://github.com/pydantic/pydantic", terms: "https://github.com/pydantic/pydantic/blob/main/LICENSE" },
  { name: "Uvicorn", role: "Local ASGI server", license: "BSD-3-Clause", credit: "Uvicorn contributors", source: "https://github.com/Kludex/uvicorn", terms: "https://github.com/Kludex/uvicorn/blob/master/LICENSE.md" },
  { name: "React / React DOM", role: "User interface", license: "MIT", credit: "Meta and React contributors", source: "https://github.com/facebook/react", terms: "https://github.com/facebook/react/blob/main/LICENSE" },
  { name: "TanStack Table", role: "Interactive data tables", license: "MIT", credit: "TanStack contributors", source: "https://github.com/TanStack/table", terms: "https://github.com/TanStack/table/blob/main/LICENSE" },
  { name: "Plotly.js", role: "Scatter plots", license: "MIT", credit: "Plotly, Inc. and contributors", source: "https://github.com/plotly/plotly.js", terms: "https://github.com/plotly/plotly.js/blob/master/LICENSE" },
  { name: "react-plotly.js", role: "React chart integration", license: "MIT", credit: "Plotly, Inc. and contributors", source: "https://github.com/plotly/react-plotly.js", terms: "https://github.com/plotly/react-plotly.js/blob/main/LICENSE" },
  { name: "Vite", role: "Frontend build tooling", license: "MIT", credit: "VoidZero and Vite contributors", source: "https://github.com/vitejs/vite", terms: "https://github.com/vitejs/vite/blob/main/packages/vite/LICENSE.md" },
  { name: "TypeScript", role: "Frontend type checking", license: "Apache-2.0", credit: "Microsoft and contributors", source: "https://github.com/microsoft/TypeScript", terms: "https://github.com/microsoft/TypeScript/blob/main/LICENSE.txt" },
  { name: "Playwright", role: "Browser tests", license: "Apache-2.0", credit: "Microsoft and contributors", source: "https://github.com/microsoft/playwright", terms: "https://github.com/microsoft/playwright/blob/main/LICENSE" },
  { name: "pytest", role: "Python tests", license: "MIT", credit: "pytest contributors", source: "https://github.com/pytest-dev/pytest", terms: "https://github.com/pytest-dev/pytest/blob/main/LICENSE" },
  { name: "pytest-cov", role: "Python test coverage", license: "MIT", credit: "pytest-cov contributors", source: "https://github.com/pytest-dev/pytest-cov", terms: "https://github.com/pytest-dev/pytest-cov/blob/master/LICENSE" },
  { name: "HTTPX", role: "API test client", license: "BSD-3-Clause", credit: "Encode contributors", source: "https://github.com/encode/httpx", terms: "https://github.com/encode/httpx/blob/master/LICENSE.md" },
  { name: "sax", role: "Indirect Plotly.js dependency", license: "BlueOak-1.0.0", credit: "sax contributors", source: "https://github.com/isaacs/sax-js", terms: "https://github.com/isaacs/sax-js/blob/main/LICENSE.md" },
  { name: "caniuse-lite", role: "Build-time browser compatibility data", license: "CC BY 4.0", credit: "caniuse.com", source: "https://github.com/browserslist/caniuse-lite", terms: "https://github.com/browserslist/caniuse-lite/blob/main/README.md#license" },
];

const testDataSources = [
  { name: "MDM2", record: "Target CHEMBL5023", url: "https://www.ebi.ac.uk/chembl/explore/target/CHEMBL5023" },
  { name: "GLP-1", record: "Compounds; original selection filter not recorded", url: "https://www.ebi.ac.uk/chembl/explore/compounds/" },
  { name: "Grb2-SH2_S", record: "Document CHEMBL1144969", url: "https://www.ebi.ac.uk/chembl/explore/document/CHEMBL1144969" },
  { name: "MEGA-A4_S", record: "Document CHEMBL5154904", url: "https://www.ebi.ac.uk/chembl/explore/document/CHEMBL5154904" },
  { name: "Oxytocin_S", record: "Document CHEMBL3244056", url: "https://www.ebi.ac.uk/chembl/explore/document/CHEMBL3244056" },
  { name: "RGD_H2T", record: "Document CHEMBL4765475", url: "https://www.ebi.ac.uk/chembl/explore/document/CHEMBL4765475" },
  { name: "UbE2", record: "Document CHEMBL1153736", url: "https://www.ebi.ac.uk/chembl/explore/document/CHEMBL1153736" },
  { name: "RAP1", record: "Document CHEMBL3751807", url: "https://www.ebi.ac.uk/chembl/explore/document/CHEMBL3751807" },
];

function LicensesPanel() {
  return <section className="panel license-panel"><div className="panel-heading"><div><h2>Licenses &amp; credits</h2></div><p>Application, local data terms, direct dependencies, and notable transitive licenses.</p></div>
    <div className="license-table-scroll"><table><thead><tr><th>Component</th><th>Use</th><th>License</th><th>Credit</th><th>Links</th></tr></thead><tbody>{licenses.map((item) => <tr key={item.name}><td>{item.name}</td><td>{item.role}</td><td>{item.license}</td><td>{item.credit}</td><td><a href={item.source} target="_blank" rel="noreferrer noopener">Project</a><a href={item.terms} target="_blank" rel="noreferrer noopener">License</a></td></tr>)}</tbody></table></div>
    <div className="license-sources"><h3>Test data sources</h3><p>These links identify the ChEMBL records used for local testing; they do not download or reproduce the local CSV selections. GLP-1's exact source filter was not recorded.</p>
      <div className="license-table-scroll"><table><thead><tr><th>Dataset</th><th>ChEMBL record</th><th>Source</th></tr></thead><tbody>{testDataSources.map((item) => <tr key={item.name}><td>{item.name}</td><td>{item.record}</td><td><a href={item.url} target="_blank" rel="noreferrer noopener">Open ChEMBL</a></td></tr>)}</tbody></table></div>
    </div>
    <p className="license-note">PepSeMan code, documentation, and icons are MIT-licensed. ChEMBL-derived test data have separate terms and are not bundled in this source-only release. Earlier copies may persist in forks, clones, or caches. See LICENSE, DATA_LICENSE.md, and docs/DEPENDENCIES.md for scope and provenance. Third-party code retains its own license.</p>
  </section>;
}

function Workspace({ title, hint, actions, children, className = "" }: { title: string; hint?: string; actions?: React.ReactNode; children: React.ReactNode; className?: string }) {
  return <section className={`panel ${className}`}><div className="panel-heading"><div><p className="eyebrow">ANALYSIS WORKSPACE</p><h2>{title}</h2></div>{actions ?? (hint ? <p>{hint}</p> : null)}</div>{children}</section>;
}

function Projects({ projects, reload, select }: { projects: Project[]; reload: () => void; select: (id: number) => void }) {
  const [name, setName] = useState("");
  const create = async (event: FormEvent) => {
    event.preventDefault();
    const project = await api<Project>("/api/projects", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ name }) });
    setName(""); reload(); select(project.project_id);
  };
  return <Workspace title="Projects" hint="Projects isolate source records, measurements, revisions, contexts, and MMP results.">
    <form className="create-form" onSubmit={create}><input aria-label="Project filename" placeholder="Filename (e.g. MDM2.csv)" value={name} onChange={(e) => setName(e.target.value)} required /><button>Create project</button></form>
    <div className="project-grid">{projects.map((project) => <button type="button" className="project-card" key={project.project_id} onClick={() => select(project.project_id)}><span className="project-id">PROJECT {project.project_id}</span><strong>{project.name}</strong><span>Updated {new Date(project.updated_at).toLocaleString()}</span></button>)}</div>
  </Workspace>;
}

function SequenceDictionary() {
  const [entries, setEntries] = useState<DictionaryEntry[]>([]);
  const [filter, setFilter] = useState("");
  const [error, setError] = useState("");
  useEffect(() => {
    api<DictionaryEntry[]>("/api/sequence-dictionary").then(setEntries).catch((reason) => setError(reason.message));
  }, []);
  const update = async (variantHash: string, abbreviation: string) => {
    const entry = await api<DictionaryEntry>(`/api/sequence-dictionary/${variantHash}`, { method: "PUT", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ abbreviation }) });
    setEntries((current) => current.map((item) => item.variant_hash === variantHash ? entry : item));
  };
  const visible = entries.filter((item) => [item.abbreviation, item.default_abbreviation, item.residue_type, item.stereo_status, item.r1, item.r2, item.r3, item.bridge].some((value) => String(value ?? "").toLowerCase().includes(filter.toLowerCase())));
  return <Workspace title="Sequence dictionary" hint="Shared across all projects. Only bracketed non-natural building blocks are editable; natural L/D amino acids and [Ac] use fixed codes.">
    <div className="dictionary-toolbar"><input type="search" aria-label="Search dictionary" placeholder="Search abbreviation or residue…" value={filter} onChange={(event) => setFilter(event.target.value)} /><span>{visible.length} / {entries.length} building blocks</span></div>
    {error && <div className="notice error">{error}</div>}
    <div className="dictionary-scroll"><table><thead><tr><th>Structure</th><th>Abbreviation</th><th>Type / stereochemistry</th><th>Side-chain fragments</th><th>Edit shared abbreviation</th></tr></thead><tbody>{visible.map((item) => <DictionaryRow key={item.variant_hash} item={item} update={update} />)}</tbody></table></div>
  </Workspace>;
}

function DictionaryRow({ item, update }: { item: DictionaryEntry; update: (variantHash: string, abbreviation: string) => Promise<void> }) {
  const [draft, setDraft] = useState(item.custom_abbreviation ?? "");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");
  useEffect(() => { setDraft(item.custom_abbreviation ?? ""); }, [item.custom_abbreviation]);
  const save = async (value: string) => {
    setBusy(true); setError("");
    try { await update(item.variant_hash, value); } catch (reason) { setError(reason instanceof Error ? reason.message : "Could not save abbreviation"); }
    finally { setBusy(false); }
  };
  return <tr><td><img loading="lazy" src={item.structure_url} alt={`Structure of ${item.abbreviation}`} /></td><td><strong>{item.abbreviation}</strong><small>Default: {item.default_abbreviation}</small></td><td>{item.residue_type}<small>{item.stereo_status}</small></td><td>{[item.r1, item.r2, item.r3, item.bridge].filter(Boolean).join(" · ") || "—"}</td><td><div className="dictionary-edit"><input aria-label={`Abbreviation for ${item.variant_hash}`} value={draft} maxLength={18} placeholder="e.g. 4FP" onChange={(event) => setDraft(event.target.value)} /><button type="button" disabled={busy || draft === (item.custom_abbreviation ?? "")} onClick={() => save(draft)}>Save</button><button type="button" className="secondary" disabled={busy || !item.custom_abbreviation} onClick={() => save("")}>Reset</button></div>{error && <small className="dictionary-error">{error}</small>}</td></tr>;
}

function ImportQc({ projectId, refresh, setRefresh, setMessage }: { projectId: number; refresh: number; setRefresh: React.Dispatch<React.SetStateAction<number>>; setMessage: (value: string) => void }) {
  const [file, setFile] = useState<File>();
  const [qcIdFilter, setQcIdFilter] = useState("");
  const upload = async (event: FormEvent) => {
    event.preventDefault(); if (!file || !projectId) return;
    const extension = file.name.split(".").pop()?.toLowerCase();
    const sourceFormat = extension === "sdf" ? "sdf" : extension === "smi" ? "smi" : "csv";
    const result = await api<{ duplicate: boolean; job_id: number | null }>(`/api/projects/${projectId}/imports?source_name=${encodeURIComponent(file.name)}&source_format=${sourceFormat}`, { method: "POST", body: await file.arrayBuffer() });
    setMessage(result.duplicate ? "Identical content detected; structural analysis cache reused." : `Import queued as persistent job ${result.job_id}.`);
    setRefresh((value) => value + 1);
  };
  return <Workspace title="Import & quality control" hint="CSV (semicolon-aware), SMI, and SDF are accepted. Originals and exclusion reasons are preserved.">
    <form className="upload" onSubmit={upload}><label><span>Source file</span><input type="file" accept=".csv,.smi,.sdf" onChange={(event) => setFile(event.target.files?.[0])} /></label><button disabled={!file || !projectId}>Queue import</button></form>
    <h3>Import history</h3><DataTable endpoint={`/api/projects/${projectId}/imports`} refreshKey={refresh} />
    <div className="qc-toolbar"><h3>QC records</h3><input type="search" aria-label="Search QC ID" placeholder="Search ID…" value={qcIdFilter} onChange={(event) => setQcIdFilter(event.target.value)} /></div>
    <DataTable endpoint={`/api/projects/${projectId}/qc?id_filter=${encodeURIComponent(qcIdFilter)}`} refreshKey={refresh} pageSize={5} showPageJump />
  </Workspace>;
}

function Scatter({ projectId, reference, choosePair, refresh }: { projectId: number; reference: string; choosePair: (row: Row) => void; refresh: number }) {
  const [rows, setRows] = useState<Row[]>([]);
  const [x, setX] = useState("delta_mw");
  const [y, setY] = useState("after_mw");
  useEffect(() => {
    let cancelled = false;
    api<Page>(`/api/projects/${projectId}/pairs?reference_hash=${encodeURIComponent(reference)}&page=1&page_size=500`)
      .then((page) => { if (!cancelled) setRows(page.rows); })
      .catch(() => { if (!cancelled) setRows([]); });
    return () => { cancelled = true; };
  }, [projectId, reference, refresh]);
  const axes = ["before_mw", "after_mw", "delta_mw", "before_activity", "after_activity", "delta_activity"];
  const plotted = rows.filter((row) => typeof row[x] === "number" && typeof row[y] === "number");
  return <div className="scatter-panel"><p>{plotted.length} pair{plotted.length === 1 ? "" : "s"} with numeric values for both axes. Choose a point to show its pair detail.</p>
    <div className="axis-controls"><label>X <select value={x} onChange={(e) => setX(e.target.value)}>{axes.map((item) => <option key={item}>{item}</option>)}</select></label><label>Y <select value={y} onChange={(e) => setY(e.target.value)}>{axes.map((item) => <option key={item}>{item}</option>)}</select></label></div>
    <Plot data={[{ type: "scattergl", mode: "markers", x: plotted.map((row) => row[x] as number), y: plotted.map((row) => row[y] as number), text: plotted.map((row) => String(row.pair_key)), marker: { color: plotted.map((row) => Number(row.site_position)), colorscale: "Viridis", size: 10 } }]} layout={{ autosize: true, height: 360, margin: { t: 20, r: 20, b: 60, l: 70 }, xaxis: { title: { text: x } }, yaxis: { title: { text: y } }, paper_bgcolor: "transparent", plot_bgcolor: "#f7faf7" }} useResizeHandler style={{ width: "100%" }} onClick={(event) => choosePair(plotted[event.points[0].pointIndex])} />
  </div>;
}

function PairDetail({ projectId, pairKey, reference, displayOrigin, setDisplayOrigin }: { projectId: number; pairKey: string; reference: string; displayOrigin: number; setDisplayOrigin: (value: number) => void }) {
  const [detail, setDetail] = useState<any>();
  const [error, setError] = useState("");
  useEffect(() => {
    let cancelled = false;
    setDetail(undefined); setError("");
    api(`/api/projects/${projectId}/pairs/${pairKey}?reference_hash=${encodeURIComponent(reference)}`)
      .then((value) => { if (!cancelled) setDetail(value); })
      .catch((reason) => { if (!cancelled) setError(reason.message); });
    return () => { cancelled = true; };
  }, [projectId, pairKey, reference]);
  if (error) return <section className="panel analysis-lower"><div className="detail-topline"><h2>Pair detail</h2></div><div className="notice error">{error}</div></section>;
  if (!detail) return <section className="panel analysis-lower"><div className="detail-topline"><h2>Pair detail</h2></div><div className="empty">Loading…</div></section>;
  const before = detail.before.residues as Row[];
  const after = detail.after.residues as Row[];
  const site = Number(detail.site_position);
  const beforeSite = before[site - 1];
  const afterSite = after[Number(detail.partner_site_position) - 1];
  const properties = detail.properties as Row;
  const links = (detail.cyclization_links ?? []) as { kind: string; first: { kind: string; position: number }; second: { kind: string; position: number } }[];
  const cyclicBackbone = detail.backbone_kind === "head_to_tail";
  return <section className="panel analysis-lower">
    <div className="detail-topline"><h2>Pair detail</h2><span>{detail.site_kind === "residue" ? `Residue ${site}` : detail.site_kind} · {String(properties.before_reference ?? "reference")} → {String(properties.after_reference ?? "pair")}</span></div>
    <div className="detail-content">
      <article className="detail-card scaffold-card"><h3>Scaffold</h3><p>{cyclicBackbone ? "Head-to-tail ring · residue numbers follow the canonical sequence" : "Parent/reference molecular structure"}</p><div className="scaffold-sequence"><span>Parent sequence</span>{links.length ? <CyclicSequence sequence={String(detail.parent_sequence ?? "")} links={links} origin={displayOrigin} residueCount={before.length} highlightPosition={site} /> : <HighlightedSequence sequence={String(detail.parent_sequence ?? "")} start={Number(detail.parent_sequence_highlight_start)} end={Number(detail.parent_sequence_highlight_end)} />}</div>{cyclicBackbone && <label className="origin-control">Display origin <select aria-label="Display origin residue" value={Math.min(displayOrigin, before.length)} onChange={(event) => setDisplayOrigin(Number(event.target.value))}>{before.map((_, index) => <option value={index + 1} key={index + 1}>Residue {index + 1}</option>)}</select></label>}{links.length > 0 && <div className="link-legend">{links.map((link, index) => <span className={`link-legend-item link-${link.kind}`} key={index}>{index + 1}: {link.kind === "olefin_metathesis" ? "Ring C=C (route unproven)" : link.kind.replaceAll("_", " ")} ({link.first.kind === "residue" ? link.first.position : link.first.kind}–{link.second.kind === "residue" ? link.second.position : link.second.kind})</span>)}</div>}<div className="detail-svg" dangerouslySetInnerHTML={{ __html: detail.scaffold_svg }} /></article>
      <BuildingBlockCard title="Original BB" value={detail.before} residue={beforeSite} siteKind={detail.site_kind} />
      <BuildingBlockCard title="Convert BB" value={detail.after} residue={afterSite} siteKind={detail.site_kind} />
      <article className="detail-card property-card"><h3>Pair difference &amp; properties</h3><dl>
        <div><dt>ChEMBL Value (pChEMBL)</dt><dd>{formatProperty(properties.before_chembl_value)} → {formatProperty(properties.after_chembl_value)}</dd></div>
        <div className="property-emphasis"><dt>ΔChEMBL Value</dt><dd>{formatProperty(properties.delta_chembl_value)}</dd></div>
        <div><dt>MW</dt><dd>{formatProperty(properties.before_mw)} → {formatProperty(properties.after_mw)}</dd></div>
        <div><dt>ΔMW</dt><dd>{formatProperty(properties.delta_mw)}</dd></div>
        <div><dt>Assay / endpoint</dt><dd>{String(properties.assay_id ?? "—")} / {String(properties.endpoint ?? "—")}</dd></div>
        {detail.site_kind === "residue" && <div><dt>R1 / R2 / R3</dt><dd>{["r1", "r2", "r3"].map((field) => `${String(beforeSite?.[field] ?? "—")} → ${String(afterSite?.[field] ?? "—")}`).join(" · ")}</dd></div>}
      </dl></article>
    </div>
  </section>;
}

function formatProperty(value: unknown) {
  return typeof value === "number" ? value.toFixed(2) : "—";
}

function BuildingBlockCard({ title, value, residue, siteKind }: { title: string; value: any; residue?: Row; siteKind: string }) {
  const cap = siteKind === "n_cap" ? value.n_cap : value.c_cap;
  return <article className="detail-card bb-card"><h3>{title}</h3>{siteKind === "residue" ? <><p>{String(residue?.residue_type ?? "Residue")} · {String(residue?.stereo_status ?? "—")}</p><div className="detail-svg" dangerouslySetInnerHTML={{ __html: value.bb_svg }} /></> : <><p>No amino-acid building block change</p><div className="cap-value">{siteKind === "n_cap" ? "N cap" : "C cap"}: {String(cap ?? "—")}</div></>}</article>;
}

function JobHistory({ projectId, refresh, setRefresh, setMessage }: { projectId: number; refresh: number; setRefresh: React.Dispatch<React.SetStateAction<number>>; setMessage: (value: string) => void }) {
  const [revisions, setRevisions] = useState<Row[]>([]);
  useEffect(() => { if (projectId) api<Page>(`/api/projects/${projectId}/revisions?page=1&page_size=500`).then((value) => setRevisions(value.rows)); }, [projectId, refresh]);
  const restore = async (revisionId: number) => {
    const value = await api<{ job_id: number }>(`/api/projects/${projectId}/revisions/${revisionId}/restore`, { method: "POST" });
    setMessage(`Revision restore queued as job ${value.job_id}.`);
    setRefresh((item) => item + 1);
  };
  return <Workspace title="Persistent job & revision history" hint="The worker runs independently of the browser. Restores create a new revision without rewriting history.">
    <DataTable endpoint={`/api/projects/${projectId}/jobs`} refreshKey={refresh} />
    <h3>Analysis revisions</h3>
    <DataTable endpoint={`/api/projects/${projectId}/revisions`} refreshKey={refresh} />
    <div className="revision-list">{revisions.map((item) => <div key={String(item.revision_id)}><span>Revision {String(item.revision_id)} · input {String(item.input_revision)} · {String(item.created_at)}</span><button type="button" className="secondary" onClick={() => restore(Number(item.revision_id))}>Restore</button></div>)}</div>
  </Workspace>;
}

