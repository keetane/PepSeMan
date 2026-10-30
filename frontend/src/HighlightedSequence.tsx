export function HighlightedSequence({ sequence, start, end }: { sequence: string; start: number; end: number }) {
  const hasHighlight = Number.isInteger(start) && Number.isInteger(end) && start >= 0 && end > start && end <= sequence.length;
  return <span className="parent-sequence-text" title={sequence}>
    {hasHighlight ? <>{sequence.slice(0, start)}<mark>{sequence.slice(start, end)}</mark>{sequence.slice(end)}</> : sequence}
  </span>;
}
