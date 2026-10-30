import { useEffect, useMemo, useRef, useState } from "react";
import {
  ColumnDef,
  SortingState,
  flexRender,
  getCoreRowModel,
  useReactTable,
} from "@tanstack/react-table";
import { api, display, Page, Row } from "./api";
import { HighlightedSequence } from "./HighlightedSequence";
import { CyclicSequence } from "./CyclicSequence";

type Props = {
  endpoint: string;
  pageSize?: number;
  showPageJump?: boolean;
  refreshKey?: number;
  onRowClick?: (row: Row) => void;
  selectedKeys?: Set<string>;
  selectionKey?: string;
  onSelectionChange?: (keys: Set<string>, rows: Row[]) => void;
  onSelectAll?: (checked: boolean) => void | Promise<void>;
  activeRowKey?: string;
  activeRowValue?: string;
  navigationLabel?: string;
  stickyColumns?: string[];
  leadingColumns?: string[];
  autoSelectFirst?: boolean;
  unpaginated?: boolean;
  resizableColumns?: boolean;
};

export function DataTable({
  endpoint,
  pageSize = 50,
  showPageJump = false,
  refreshKey = 0,
  onRowClick,
  selectedKeys,
  selectionKey,
  onSelectionChange,
  onSelectAll,
  activeRowKey,
  activeRowValue,
  navigationLabel,
  stickyColumns = [],
  leadingColumns = [],
  autoSelectFirst = false,
  unpaginated = false,
  resizableColumns = false,
}: Props) {
  const [data, setData] = useState<Page>({ rows: [], total: 0, page: 1, page_size: pageSize });
  const [sorting, setSorting] = useState<SortingState>([]);
  const [page, setPage] = useState(1);
  const [jumpPage, setJumpPage] = useState("1");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");
  const [columnWidths, setColumnWidths] = useState<Record<string, number>>({});
  const resizeCleanup = useRef<(() => void) | null>(null);
  useEffect(() => () => resizeCleanup.current?.(), []);
  const selectionRef = useRef({ onRowClick, activeRowValue });
  const pendingNavigation = useRef<1 | -1 | null>(null);
  const tableScrollRef = useRef<HTMLDivElement>(null);
  useEffect(() => { selectionRef.current = { onRowClick, activeRowValue }; }, [onRowClick, activeRowValue]);
  const sort = sorting.map((item) => `${item.id}:${item.desc ? "desc" : "asc"}`).join(",");
  useEffect(() => { setPage(1); setJumpPage("1"); }, [endpoint]);
  useEffect(() => { setJumpPage(String(page)); }, [page]);

  useEffect(() => {
    const controller = new AbortController();
    setLoading(true);
    setError("");
    const separator = endpoint.includes("?") ? "&" : "?";
    api<Page>(`${endpoint}${separator}page=${page}&page_size=${pageSize}&sort=${encodeURIComponent(sort)}${unpaginated ? "&all_rows=true" : ""}`, {
      signal: controller.signal,
    })
      .then((result) => {
        setData(result);
        const target = pendingNavigation.current;
        pendingNavigation.current = null;
        const { onRowClick: choose, activeRowValue: active } = selectionRef.current;
        if (!choose || !result.rows.length) return;
        if (target) {
          choose(result.rows[target === 1 ? 0 : result.rows.length - 1]);
        } else if (autoSelectFirst && activeRowKey && !result.rows.some((row) => String(row[activeRowKey]) === active)) {
          choose(result.rows[0]);
        }
      })
      .catch((reason) => {
        if (reason.name !== "AbortError") setError(reason.message);
      })
      .finally(() => setLoading(false));
    return () => controller.abort();
  }, [endpoint, page, pageSize, sort, refreshKey, autoSelectFirst, activeRowKey, unpaginated]);

  const columns = useMemo<ColumnDef<Row>[]>(() => {
    const hidden = new Set(["cyclization_links", "topology_signature", "canonical_start", "partner_site_position", "residue_count"]);
    const names = data.rows.length ? Object.keys(data.rows[0]).filter((name) => !name.startsWith("parent_sequence_highlight_") && !hidden.has(name)) : [];
    const orderedNames = [...leadingColumns.filter((name) => names.includes(name)), ...names.filter((name) => !leadingColumns.includes(name))];
    const labels: Record<string, string> = { source_molecule_id: "ID", n_pairs: "nPairs", value: "Value", parent_sequence: "Parent sequence", original: "Original", convert: "Convert", site_position: "Site position", delta_activity: "ΔValue", delta_mw: "ΔMW" };
    const dynamic = orderedNames.map((name): ColumnDef<Row> => ({
      id: name,
      accessorFn: (row) => row[name],
      header: labels[name] ?? name.replaceAll("_", " "),
      cell: (info) => {
        const text = display(info.getValue());
        if (name === "cyclization_types") {
          return <span title={text}>{text.replaceAll("head_to_tail", "Head-to-tail").replaceAll("olefin_metathesis", "Ring C=C (route unproven)").replaceAll("disulfide", "Disulfide").replaceAll("thioether", "Thioether").replaceAll("triazole", "Triazole bridge")}</span>;
        }
        if (name === "scaffold_preview" && typeof info.getValue() === "string") {
          return <img className="scaffold-preview" src={String(info.getValue())} alt="Highlighted pair scaffold preview" />;
        }
        if (name === "sequence" && Array.isArray(info.row.original.cyclization_links) && info.row.original.cyclization_links.length) {
          return <CyclicSequence sequence={text} links={info.row.original.cyclization_links as any} residueCount={Number(info.row.original.length)} />;
        }
        if (name === "parent_sequence" && Array.isArray(info.row.original.cyclization_links) && info.row.original.cyclization_links.length) {
          return <CyclicSequence sequence={text} links={info.row.original.cyclization_links as any} residueCount={Number(info.row.original.residue_count)} highlightPosition={Number(info.row.original.site_position)} />;
        }
        if (name === "parent_sequence") {
          return <HighlightedSequence sequence={text} start={Number(info.row.original.parent_sequence_highlight_start)} end={Number(info.row.original.parent_sequence_highlight_end)} />;
        }
        return <span title={text}>{text}</span>;
      },
    }));
    if (!selectionKey || !onSelectionChange || !selectedKeys) return dynamic;
    return [
      {
        id: "select",
        enableSorting: false,
        header: () => {
          const keys = data.rows.map((row) => String(row[selectionKey]));
          const selectedCount = onSelectAll ? selectedKeys.size : keys.filter((key) => selectedKeys.has(key)).length;
          const totalCount = onSelectAll ? data.total : keys.length;
          return <input
            type="checkbox"
            aria-label={onSelectAll ? "Select all MMP pairs" : "Select all rows on this page"}
            title={onSelectAll ? "Select all MMP pairs" : "Select all rows on this page"}
            disabled={totalCount === 0}
            checked={totalCount > 0 && selectedCount === totalCount}
            ref={(input) => { if (input) input.indeterminate = selectedCount > 0 && selectedCount < totalCount; }}
            onClick={(event) => event.stopPropagation()}
            onChange={(event) => {
              if (onSelectAll) { void onSelectAll(event.target.checked); return; }
              const next = new Set(selectedKeys);
              for (const key of keys) event.target.checked ? next.add(key) : next.delete(key);
              onSelectionChange(next, data.rows);
            }}
          />;
        },
        cell: ({ row }) => {
          const key = String(row.original[selectionKey]);
          return (
            <input
              aria-label={`Select ${key}`}
              type="checkbox"
              checked={selectedKeys.has(key)}
              onClick={(event) => event.stopPropagation()}
              onChange={() => {
                const next = new Set(selectedKeys);
                next.has(key) ? next.delete(key) : next.add(key);
                onSelectionChange(next, data.rows);
              }}
            />
          );
        },
      },
      ...dynamic,
    ];
  }, [data.rows, data.total, leadingColumns, onSelectionChange, onSelectAll, selectedKeys, selectionKey]);

  const table = useReactTable({
    data: data.rows,
    columns,
    state: { sorting },
    onSortingChange: (next) => {
      setSorting(next);
      setPage(1);
    },
    manualSorting: true,
    enableMultiSort: true,
    getCoreRowModel: getCoreRowModel(),
  });
  const pages = Math.max(1, Math.ceil(data.total / data.page_size));
  const requestedPage = Number(jumpPage);
  const validJumpPage = Number.isInteger(requestedPage) && requestedPage >= 1 && requestedPage <= pages;
  const activeIndex = activeRowKey ? data.rows.findIndex((row) => String(row[activeRowKey]) === activeRowValue) : -1;
  const moveRow = (direction: 1 | -1) => {
    const nextIndex = activeIndex + direction;
    if (nextIndex >= 0 && nextIndex < data.rows.length) {
      onRowClick?.(data.rows[nextIndex]);
      requestAnimationFrame(() => tableScrollRef.current?.querySelectorAll("tbody tr")[nextIndex]?.scrollIntoView({ block: "nearest", inline: "nearest" }));
    } else if (direction === 1 && page < pages) {
      pendingNavigation.current = 1;
      setPage((value) => value + 1);
    } else if (direction === -1 && page > 1) {
      pendingNavigation.current = -1;
      setPage((value) => value - 1);
    }
  };
  const handleRowKeyDown = (event: React.KeyboardEvent<HTMLDivElement>) => {
    if (event.target !== event.currentTarget || !onRowClick || !activeRowKey || loading) return;
    if (event.key !== "ArrowUp" && event.key !== "ArrowDown") return;
    event.preventDefault();
    moveRow(event.key === "ArrowDown" ? 1 : -1);
  };
  const pinnedWidths: Record<string, number> = { source_molecule_id: 112, n_pairs: 82, value: 94 };
  const pinnedStyle = (columnId: string) => {
    const index = stickyColumns.indexOf(columnId);
    if (index < 0) return undefined;
    return { left: stickyColumns.slice(0, index).reduce((total, name) => total + (pinnedWidths[name] ?? 100), 0), width: pinnedWidths[columnId] ?? 100, minWidth: pinnedWidths[columnId] ?? 100 };
  };
  const defaultWidth = (columnId: string) => {
    if (columnId === "select") return 56;
    if (columnId === "parent_sequence") return Math.max(300, ...data.rows.map((row) => String(row.parent_sequence ?? "").length * 7 + 32));
    if (columnId === "scaffold_preview") return 210;
    if (columnId.endsWith("smiles")) return 320;
    if (["original", "convert"].includes(columnId)) return 160;
    return 150;
  };
  const widthOf = (columnId: string) => columnWidths[columnId] ?? defaultWidth(columnId);
  const columnStyle = (columnId: string) => resizableColumns ? { width: widthOf(columnId), minWidth: widthOf(columnId), maxWidth: widthOf(columnId) } : pinnedStyle(columnId);
  const startResize = (event: React.PointerEvent, columnId: string) => {
    event.preventDefault(); event.stopPropagation();
    resizeCleanup.current?.();
    const startX = event.clientX;
    const startWidth = widthOf(columnId);
    const move = (pointer: PointerEvent) => setColumnWidths((current) => ({ ...current, [columnId]: Math.max(56, startWidth + pointer.clientX - startX) }));
    const finish = () => { window.removeEventListener("pointermove", move); window.removeEventListener("pointerup", finish); resizeCleanup.current = null; };
    window.addEventListener("pointermove", move);
    window.addEventListener("pointerup", finish);
    resizeCleanup.current = finish;
  };
  const resizeWithKeyboard = (event: React.KeyboardEvent, columnId: string) => {
    if (event.key !== "ArrowLeft" && event.key !== "ArrowRight") return;
    event.preventDefault();
    setColumnWidths((current) => ({ ...current, [columnId]: Math.max(56, widthOf(columnId) + (event.key === "ArrowRight" ? 20 : -20)) }));
  };
  const tableWidth = resizableColumns ? table.getAllLeafColumns().reduce((sum, column) => sum + widthOf(column.id), 0) : undefined;

  return (
    <div className="table-shell">
      {error && <div className="notice error">{error}</div>}
      <div className="table-scroll" aria-busy={loading} aria-label={navigationLabel} role={navigationLabel ? "region" : undefined} tabIndex={onRowClick && activeRowKey ? 0 : undefined} onKeyDown={handleRowKeyDown} ref={tableScrollRef}>
        <table className={[stickyColumns.length ? "has-sticky-columns" : "", resizableColumns ? "resizable-table" : ""].join(" ")} style={resizableColumns ? { width: tableWidth, tableLayout: "fixed" } : undefined}>
          <thead>
            {table.getHeaderGroups().map((group) => (
              <tr key={group.id}>
                {group.headers.map((header) => (
                  <th key={header.id} className={stickyColumns.includes(header.column.id) ? "sticky-column" : ""} style={columnStyle(header.column.id)}>
                    {header.column.id === "select" ? flexRender(header.column.columnDef.header, header.getContext()) : <button
                      type="button"
                      className="sort-button"
                      disabled={!header.column.getCanSort()}
                      onClick={header.column.getToggleSortingHandler()}
                    >
                      {flexRender(header.column.columnDef.header, header.getContext())}
                      {{ asc: " ↑", desc: " ↓" }[header.column.getIsSorted() as string] ?? ""}
                    </button>}
                    {resizableColumns && <span className="column-resizer" role="separator" tabIndex={0} aria-label={`Resize ${header.column.id} column`} onPointerDown={(event) => startResize(event, header.column.id)} onKeyDown={(event) => resizeWithKeyboard(event, header.column.id)} />}
                  </th>
                ))}
              </tr>
            ))}
          </thead>
          <tbody>
            {table.getRowModel().rows.map((row) => (
              <tr key={row.id} onClick={() => { onRowClick?.(row.original); if (onRowClick && activeRowKey) tableScrollRef.current?.focus(); }} className={[onRowClick ? "clickable" : "", activeRowKey && activeRowValue && row.original[activeRowKey] === activeRowValue ? "selected-row" : ""].join(" ")}>
                {row.getVisibleCells().map((cell) => (
                  <td key={cell.id} className={stickyColumns.includes(cell.column.id) ? "sticky-column" : ""} style={columnStyle(cell.column.id)}>{flexRender(cell.column.columnDef.cell, cell.getContext())}</td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
        {!loading && !data.rows.length && <div className="empty">No rows</div>}
      </div>
      <div className="pagination">
        <span>{data.total.toLocaleString()} rows</span>
        {activeRowKey && onRowClick && <div className="row-navigation"><span>{activeIndex >= 0 ? (page - 1) * data.page_size + activeIndex + 1 : 0} / {data.total}</span><button type="button" aria-label="Previous row" title="Previous row" disabled={loading || !data.rows.length || (page === 1 && activeIndex <= 0)} onClick={() => moveRow(-1)}>↑</button><button type="button" aria-label="Next row" title="Next row" disabled={loading || !data.rows.length || (page >= pages && activeIndex >= data.rows.length - 1)} onClick={() => moveRow(1)}>↓</button></div>}
        {!unpaginated && <><button type="button" disabled={page <= 1} onClick={() => setPage((value) => value - 1)}>Previous</button>
        <span>{page} / {pages}</span>
        <button type="button" disabled={page >= pages} onClick={() => setPage((value) => value + 1)}>Next</button></>}
        {showPageJump && !unpaginated && <form className="page-jump" onSubmit={(event) => { event.preventDefault(); if (validJumpPage) setPage(requestedPage); }}><input type="number" aria-label="Jump to page" min={1} max={pages} value={jumpPage} onChange={(event) => setJumpPage(event.target.value)} /><button type="submit" disabled={!validJumpPage}>Go</button></form>}
      </div>
    </div>
  );
}

