type Link = {
  kind: string;
  first: { kind: string; position: number };
  second: { kind: string; position: number };
};

const labels: Record<string, string> = {
  head_to_tail: "head-to-tail",
  disulfide: "disulfide",
  thioether: "thioether",
  olefin_metathesis: "ring C=C (structure-inferred)",
  triazole: "triazole bridge (structure-inferred)",
};

export function CyclicSequence({ sequence, links, origin = 1, highlightPosition, residueCount }: {
  sequence: string;
  links: Link[];
  origin?: number;
  highlightPosition?: number;
  residueCount?: number;
}) {
  const cyclic = links.some((link) => link.kind === "head_to_tail");
  const tokens = sequence.match(/\[[^\]]+\]|./g) ?? [];
  const start = cyclic && tokens.length ? Math.max(0, Math.min(tokens.length - 1, origin - 1)) : 0;
  const rotated = tokens.slice(start).concat(tokens.slice(0, start));
  const capOffset = cyclic ? 0 : /^\[(?:Ac|H|N[A-F0-9]+)\]/.test(sequence) ? 1 : 0;
  const shownPosition = (index: number) => cyclic ? ((index + start) % tokens.length) + 1 : index - capOffset + 1;
  return <span className="cyclic-sequence" title={sequence}>
    {cyclic && <span className="cycle-symbol" title={`Head-to-tail ring · display origin ${start + 1}`}>↻</span>}
    {rotated.map((token, index) => {
      const position = shownPosition(index);
      const isResidue = position >= 1 && position <= (residueCount ?? tokens.length);
      const markers = !isResidue ? [] : links.flatMap((link, linkIndex) => {
        if (link.kind === "head_to_tail") return [];
        return [link.first, link.second].some((end) => end.kind === "residue" && end.position === position)
          ? [{ link, linkIndex }] : [];
      });
      return <span className="sequence-residue" key={`${position}-${index}`}>
        {highlightPosition === position ? <mark>{token}</mark> : token}
        {markers.map(({ link, linkIndex }) => <sup key={linkIndex} className={`link-marker link-${link.kind}`} title={`${labels[link.kind] ?? link.kind} · bridge ${linkIndex + 1}`}>{linkIndex + 1}</sup>)}
      </span>;
    })}
  </span>;
}
