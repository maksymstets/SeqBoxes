#!/usr/bin/env python3
"""
seqboxes - highlight user-defined subsequences of a FASTA/FASTQ record in
coloured boxes .

Usage
-----
  GUI:       python seqboxes.py proteins.fasta   (or run without arguments to be prompted)
  Headless:  python seqboxes.py proteins.fasta [--record Q9SMJ4] \
                 --enzyme Trypsin --enzyme "Pepsin (pH 1.3)=green" \
                 --rule "LAG=#d62728" --rule "GF=gold:fill" \
                 --png-dir out_pngs --csv out.csv
  Restriction sites in DNA records:
             python seqboxes.py plasmid.fasta --restriction EcoRI --restriction BsaI=blue
  Circular DNA (plasmids): add --circular, or tick the box in the GUI
  Lists:     python seqboxes.py --list-enzymes        (proteases)
             python seqboxes.py --list-restriction    (REBASE enzymes in Bio.Restriction)
  Site counts (valid numbers for --site, per sequence):
             python seqboxes.py proteins.fasta --enzyme Trypsin --count-sites
  Limit sites per sequence (numbering restarts for every sequence and protease):
             ... --enzyme ProteinaseK --max-sites 5      first 5 sites
             ... --enzyme Trypsin --site 3               only site number 3
             ... --enzyme Trypsin --site 1,4,7-9         several sites

Multi-sequence files: the same rules are applied to every record. The GUI shows
them stacked in one preview; saving writes one PNG per sequence (<SeqID>.png).
Use --record (repeatable) to restrict to some records.

Protease mode boxes the two residues flanking every cleaved bond (P1 and P1'),
using the Expasy PeptideCutter rules shipped with Pyteomics. Restriction mode
boxes the recognition site of each Bio.Restriction enzyme on either strand (cut
positions go to the CSV). Proteases act on protein records, restriction enzymes
on DNA records (detected automatically). The legend under the last row of each
PNG lists every enzyme (and motif) with its box colour.

Rule syntax for --rule:  SUBSEQUENCE=COLOUR[:box|fill]
Dependencies: biopython, pyteomics, matplotlib  
"""
import argparse
import csv
import os
import re
import sys
from dataclasses import dataclass
from typing import Optional

import matplotlib
from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_letters
from Bio.Restriction import AllEnzymes, CommOnly
from Bio.Seq import Seq
from matplotlib.colors import is_color_like, to_hex
from matplotlib.patches import Rectangle
from pyteomics import parser as pyparser

# ----------------------------------------------------------------- parsing


def read_records(path):
    """Read every record of a FASTA or FASTQ file with Biopython.

    The format is chosen from the first character of the file: ``>`` means FASTA,
    ``@`` means FASTQ (4-line records).

    Args:
        path: Path of the file to read.

    Returns:
        A list of ``(seq_id, sequence)`` tuples in file order. ``seq_id`` is
        Biopython's ``record.id``, i.e. the first word of the header line.

    Raises:
        ValueError: If the file is neither FASTA nor FASTQ, is malformed, or
            contains no records.
    """
    with open(path, encoding="utf-8-sig") as fh:
        first = fh.read(1)
    fmt = {">": "fasta", "@": "fastq"}.get(first)
    if fmt is None:
        raise ValueError("Not a FASTA ('>') or FASTQ ('@') file.")
    try:
        records = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(path, fmt)]
    except ValueError as e:  # Biopython raises ValueError on malformed input
        raise ValueError(f"Could not parse {fmt.upper()} file: {e}") from e
    if not records:
        raise ValueError("No records found in file.")
    return records


# ------------------------------------------------------------------ enzymes
# Proteases:           Pyteomics - pyteomics.parser.expasy_rules holds the Expasy
#                      PeptideCutter cleavage rules (regular expressions).
# Restriction enzymes: Biopython - Bio.Restriction (REBASE database).

PROTEASES = sorted(k for k in pyparser.expasy_rules if not k.endswith("_exception"))
RESTRICTION_ENZYMES = {str(e): e for e in AllEnzymes}
RESTRICTION_NAMES = sorted(RESTRICTION_ENZYMES, key=str.lower)
COMMERCIAL_NAMES = sorted((str(e) for e in CommOnly), key=str.lower)


END_PAD = "XXXX"

_DNA_LETTERS = frozenset(ambiguous_dna_letters)


def is_dna(seq):
    """Guess whether a sequence is DNA (True) or protein (False).

    A sequence counts as DNA when it consists only of IUPAC DNA letters and at least
    90 % of its letters are A, C, G, T or N. Case is ignored. A protein made up
    only of such letters would be misclassified; real proteins virtually never are.

    Args:
        seq: Sequence string.

    Returns:
        True for DNA, False for anything else, including the empty string.
    """
    s = seq.upper()
    if not s or not set(s) <= _DNA_LETTERS:
        return False
    return sum(s.count(c) for c in "ACGTN") / len(s) >= 0.9


def resolve_name(text, names, what="enzyme"):
    """Map user input to a canonical enzyme name.

    Tried in order: exact match, case-insensitive match, unique case-insensitive
    prefix.

    Args:
        text: Name as typed by the user.
        names: The valid canonical names.
        what: Word used in error messages, for example ``"protease"``.

    Returns:
        The canonical name.

    Raises:
        ValueError: If nothing matches, or the prefix is ambiguous (the message
            then lists the candidates).
    """
    if text in names:
        return text
    t = text.strip().lower()
    exact = [n for n in names if n.lower() == t]
    if len(exact) == 1:
        return exact[0]
    hits = [n for n in names if n.lower().startswith(t)] if t else []
    if len(hits) == 1:
        return hits[0]
    if not hits:
        raise ValueError(f"Unknown {what} '{text}'.")
    more = " ..." if len(hits) > 12 else ""
    raise ValueError(f"Ambiguous {what} '{text}': " + "; ".join(hits[:12]) + more)


def protease_sites(seq, protease):
    """Find every cleaved bond of a protease in a protein sequence.

    The rules are Pyteomics' Expasy PeptideCutter rules (``pyparser.icleave``,
    together with the ``<name>_exception`` rule when one exists, e.g. for trypsin).
    The sequence is padded with ``END_PAD`` before matching so that residues missing
    beyond a terminus count as allowed; only bonds whose P1 and P1' are both real
    residues are returned.

    Args:
        seq: Protein sequence (any case).
        protease: Canonical name from ``PROTEASES``.

    Returns:
        A list of ``(start, end, cut_after)`` tuples in N- to C-terminal order.
        ``start`` and ``end`` are 0-based half-open indices covering P1 and P1';
        ``cut_after`` is the 1-based number of P1, i.e. the cut falls after that
        residue.
    """
    s = seq.upper()
    exception = protease + "_exception"
    if exception not in pyparser.expasy_rules:
        exception = None
    pad = len(END_PAD)
    sites = []
    for start, _ in list(pyparser.icleave(END_PAD + s + END_PAD, protease,
                                          exception=exception))[1:]:
        p1_next = start - pad  # 0-based index of P1'
        if 1 <= p1_next <= len(s) - 1:
            sites.append((p1_next - 1, p1_next + 1, p1_next))
    return sites


def restriction_sites(seq, enzyme, circular=False):
    """Find the sites of a Bio.Restriction enzyme in a DNA sequence.

    Recognition sites on either strand come from the enzyme's own ``compsite``
    regular expression, so degenerate sites such as ``GATNNNNATC`` work. Each site
    is paired with its cut position and kept only if Bio.Restriction itself reports
    that cut (so, as in the library, a site in a linear sequence whose cut would
    fall outside the sequence is dropped). For enzymes without a cut position in
    REBASE every recognition site is returned with ``cut_after`` set to None.

    In a circular sequence the first ``size - 1`` bases are appended before
    searching, so a site may cross the origin; its ``end`` is then larger than
    ``len(seq)`` and the site continues at the start.

    Args:
        seq: DNA sequence (any case).
        enzyme: Canonical name from ``RESTRICTION_NAMES``.
        circular: Treat the sequence as circular.

    Returns:
        A sorted list of ``(start, end, cut_after)`` tuples. ``start`` and ``end``
        are 0-based half-open indices of the recognition site; ``cut_after`` is the
        1-based number of the last base before the cut on the top strand (it wraps
        to ``len(seq)`` when a circular cut lies before base 1), or None.
    """
    enz = RESTRICTION_ENZYMES[enzyme]
    s = seq.upper()
    n = len(s)
    if n == 0:
        return []
    # a circular sequence is searched with its first size-1 bases appended, so that
    # recognition sites crossing the origin are found; each site is kept once
    text = s + s[:enz.size - 1] if circular else s
    matches = [m for m in enz.compsite.finditer(text) if m.start() < n]
    if enz.fst5 is None or enz.fst3 is None:
        # REBASE lists the recognition sequence but no cut position for this
        # enzyme: box every recognition site and leave cut_after empty.
        return sorted((m.start(), m.start() + enz.size, None) for m in matches)
    # 1-based number of the first base after a cut (Bio.Restriction convention)
    cuts = set(enz.search(Seq(s), linear=not circular))
    sites = []
    for m in matches:
        if m.lastgroup.endswith("_as"):  # site found on the bottom strand
            cut = m.start() + 1 - enz.fst3
        else:
            cut = m.start() + 1 + enz.fst5
        if circular:
            cut = (cut - 1) % n + 1  # wrap around the origin
        if cut in cuts:
            sites.append((m.start(), m.start() + enz.size, (cut - 2) % n + 1))
    return sorted(sites)


# ----------------------------------------------------------------- layout


def parse_site_numbers(text):
    """Parse a site-number specification such as ``3`` or ``1,4,7-9``.

    Args:
        text: Comma-separated numbers and inclusive ranges; spaces are ignored.

    Returns:
        A frozenset of 1-based site numbers.

    Raises:
        ValueError: If the text is empty, contains something that is not a number
            or range, uses a number below 1, a downward range, or a number above
            100000.
    """
    nums = set()
    for part in text.replace(" ", "").split(","):
        if not part:
            continue
        lo, sep, hi = part.partition("-")
        try:
            a = int(lo)
            b = int(hi) if sep else a
        except ValueError:
            raise ValueError(f"'{part}' is not a number or a range like 4-7") from None
        if a < 1 or b < a:
            raise ValueError(f"'{part}' is not valid (site numbers start at 1; "
                             "ranges must go upwards)")
        if b > 100000:
            raise ValueError(f"'{part}' is too large")
        nums.update(range(a, b + 1))
    if not nums:
        raise ValueError("no site numbers given")
    return frozenset(nums)


def format_site_numbers(nums):
    """Format site numbers compactly; the inverse of ``parse_site_numbers``.

    Runs of three or more consecutive numbers become ranges, so
    ``frozenset({1, 4, 7, 8, 9})`` gives ``"1,4,7-9"`` and ``{2, 3}`` gives
    ``"2,3"``.

    Args:
        nums: An iterable of positive integers.

    Returns:
        The comma-separated string.
    """
    runs, run = [], []
    for n in sorted(nums):
        if run and n == run[-1] + 1:
            run.append(n)
        else:
            if run:
                runs.append(run)
            run = [n]
    if run:
        runs.append(run)
    return ",".join(str(r[0]) if len(r) == 1 else f"{r[0]},{r[1]}" if len(r) == 2
                    else f"{r[0]}-{r[-1]}" for r in runs)


@dataclass
class Rule:
    """One thing to box: a literal subsequence or the sites of an enzyme.

    ``kind`` selects the behaviour: ``"motif"`` boxes every occurrence of
    ``subseq``; ``"protease"`` boxes P1 and P1' of every bond cleaved by
    ``enzyme`` (protein records only); ``"restriction"`` boxes every recognition
    site of the Bio.Restriction enzyme ``enzyme`` (DNA records only).

    Enzyme sites are numbered 1, 2, 3 ... from the N-terminus / 5' end, restarting
    for every sequence. ``only_sites`` keeps just those numbers; otherwise
    ``max_sites`` keeps the first N; with neither, every site is boxed.

    Attributes:
        subseq: Literal subsequence (motif rules only).
        color: Any matplotlib colour; stored as given.
        style: ``"box"`` (outline) or ``"fill"`` (highlight).
        enzyme: Canonical enzyme name (enzyme rules only).
        kind: ``"motif"``, ``"protease"`` or ``"restriction"``.
        max_sites: Keep only the first N sites of each sequence, or None.
        only_sites: Keep only these site numbers (a frozenset), or None.
    """
    subseq: str
    color: str = "#d62728"
    style: str = "box"  # "box" or "fill"
    enzyme: str = ""
    kind: str = "motif"  # "motif" | "protease" | "restriction"
    max_sites: Optional[int] = None
    only_sites: Optional[frozenset] = None

    @property
    def label(self):
        """Text shown in the legend and the status line.

        Motifs show their upper-case subsequence. Enzymes show their name, followed by
        the active site filter, e.g. ``trypsin (first 3 sites)`` or
        ``EcoRI (sites 1,4-6)``.

        Returns:
            The label string (empty for a motif rule without a subsequence).
        """
        if not self.enzyme:
            return self.subseq.strip().upper()
        if self.only_sites:
            plural = "s" if len(self.only_sites) > 1 else ""
            return f"{self.enzyme} (site{plural} {format_site_numbers(self.only_sites)})"
        if self.max_sites is not None:
            plural = "" if self.max_sites == 1 else "s"
            return f"{self.enzyme} (first {self.max_sites} site{plural})"
        return self.enzyme

    def applies_to(self, seq):
        """Tell whether this rule can act on a given sequence.

        Proteases act only on protein records, restriction enzymes only on DNA records
        (see ``is_dna``); motifs act on both.

        Args:
            seq: Sequence string.

        Returns:
            True if the rule is applicable.
        """
        if self.kind == "protease":
            return not is_dna(seq)
        if self.kind == "restriction":
            return is_dna(seq)
        return True

    def all_sites(self, seq, circular=False):
        """List every site of the enzyme, before any ``max_sites``/``only_sites`` limit.

        Args:
            seq: Sequence string.
            circular: Treat DNA as circular (ignored for proteases).

        Returns:
            ``[(start, end, cut_after)]`` as described in ``protease_sites`` and
            ``restriction_sites``; an empty list for a motif rule.
        """
        if self.kind == "protease":
            return protease_sites(seq, self.enzyme)
        if self.kind == "restriction":
            return restriction_sites(seq, self.enzyme, circular and is_dna(seq))
        return []

    def find(self, seq, circular=False):
        """List the spans to box in a sequence, with the site filter applied.

        Args:
            seq: Sequence string.
            circular: Treat DNA as circular; matches may then cross the origin. Has no
                effect on protein records.

        Returns:
            ``[(start, end, site_number, cut_after)]`` with 0-based half-open
            ``start``/``end``. ``site_number`` is the site's ordinal among ALL sites of
            this enzyme in ``seq`` (so it is unchanged by the filter); both it and
            ``cut_after`` are None for motifs. A span crossing the origin of a circular
            sequence has ``end > len(seq)``. Returns an empty list if the rule does
            not apply to ``seq``.
        """
        if not self.applies_to(seq):
            return []
        circular = circular and is_dna(seq)
        if self.kind == "motif":
            sub = self.subseq.strip()
            return [(m, m + len(sub), None, None)
                    for m in find_matches(seq, sub, circular)]
        numbered = [(a, b, n, c) for n, (a, b, c) in
                    enumerate(self.all_sites(seq, circular), 1)]
        if self.only_sites:
            return [t for t in numbered if t[2] in self.only_sites]
        if self.max_sites is not None:
            return numbered[:self.max_sites]
        return numbered


def site_count_table(records, rules, circular=False):
    """Count the sites of every enzyme rule in every record.

    All sites are counted, before any ``max_sites``/``only_sites`` limit, so the
    result tells the user which site numbers exist. Rows are produced only for
    enzyme rules that apply to the record (proteases to proteins, restriction
    enzymes to DNA).

    Args:
        records: ``[(seq_id, sequence)]``.
        rules: Rules to evaluate; motif rules are ignored.
        circular: Treat DNA records as circular.

    Returns:
        ``[(seq_id, rule, positions)]`` where ``positions`` holds the 1-based
        first residue/base of each boxed span; the k-th entry is site number k.
    """
    return [(rid, r, [a + 1 for a, _, _ in r.all_sites(seq, circular)])
            for rid, seq in records for r in rules if r.enzyme and r.applies_to(seq)]


PALETTE = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf",
           "#8c564b", "#e377c2"]


LINE_LEN = 50      # residues per row
GROUP = 10         # residues per block
GAP = 1.0          # width of the space between blocks (in cell units)
ROW_H = 2.4        # row height (cell units)
LEFT_MARGIN = 1.5
TOP_MARGIN = 3.0


def cell_x(i):
    """Return the left edge, in cell units, of residue ``i`` in the drawing.

    Residues are laid out ``LINE_LEN`` per row in blocks of ``GROUP`` separated by a
    gap of ``GAP`` cells.

    Args:
        i: 0-based residue index.

    Returns:
        The x coordinate in cell units.
    """
    col = i % LINE_LEN
    return LEFT_MARGIN + col + (col // GROUP) * GAP


def cell_row(i):
    """Return the 0-based drawing row of residue ``i`` (``LINE_LEN`` per row)."""
    return i // LINE_LEN


def find_matches(seq, sub, circular=False):
    """Find every occurrence of a subsequence, overlapping ones included.

    Matching is case-insensitive and literal (no wildcards).

    Args:
        seq: Sequence to search.
        sub: Subsequence to look for.
        circular: Also find matches that cross the origin.

    Returns:
        The 0-based start of each match; empty if ``sub`` is empty. Each match is
        reported once, by its start inside ``seq``.
    """
    if not sub:
        return []
    s, sub = seq.upper(), sub.upper()
    text = s + s[:len(sub) - 1] if circular else s
    return [m.start() for m in re.finditer(f"(?=({re.escape(sub)}))", text)
            if m.start() < len(s)]


def build_segments(seq, rules, y_off=0.0, circular=False):
    """Turn rule matches into drawable box segments for one sequence.

    A match that runs past the end of a row is split into one segment per row, and
    a match crossing the origin of a circular sequence is split between the last
    and the first row. Overlapping rules are inset slightly so their boxes nest.

    Args:
        seq: Sequence string.
        rules: Rules to apply; those that do not apply to ``seq`` are skipped.
        y_off: Vertical offset of this sequence's block in the figure, in cell
            units.
        circular: Treat DNA as circular.

    Returns:
        A list of dicts, one per segment, with the keys ``rule``, ``rule_index``
        (position in ``rules``), ``start`` and ``end`` (1-based, inclusive; for a
        wrapped match ``end < start``), ``residues`` (the whole matched text),
        ``site_number``, ``cut_after``, ``seg_start``, ``seg_end``, ``row`` (all
        1-based) and the rectangle ``x0``, ``y0``, ``x1``, ``y1`` in cell units.
    """
    segs = []
    n = len(seq)
    active = [(i, r) for i, r in enumerate(rules) if r.applies_to(seq)]
    step = min(0.07, 0.4 / max(len(active) - 1, 1))  # inset per rule -> nested boxes
    for k, (r_idx, rule) in enumerate(active):
        for start, end, site_no, cut_after in rule.find(seq, circular):  # end exclusive
            wraps = end > n
            residues = ((seq + seq)[start:end] if wraps else seq[start:end]).upper()
            pieces = [(start, n), (0, end - n)] if wraps else [(start, end)]
            for p_start, p_end in pieces:
                pos = p_start
                while pos < p_end:
                    row = cell_row(pos)
                    seg_end = min(p_end, (row + 1) * LINE_LEN)
                    inset = step * k
                    x0 = cell_x(pos) + inset
                    x1 = cell_x(seg_end - 1) + 1 - inset
                    y0 = y_off + TOP_MARGIN + row * ROW_H + inset
                    y1 = y0 + 1.5 - 2 * inset
                    segs.append(dict(rule=rule, rule_index=r_idx, start=start + 1,
                                     end=(end - 1) % n + 1, residues=residues,
                                     site_number=site_no, cut_after=cut_after,
                                     seg_start=pos + 1, seg_end=seg_end, row=row + 1,
                                     x0=x0, y0=y0, x1=x1, y1=y1))
                    pos = seg_end
    return segs


# --------------------------------------------------------------- rendering

BLOCK_GAP = 1.5      # empty space between records (cell units)
MAX_PX = 60000       # keep PNG dimensions below matplotlib/Agg limits
LEGEND_ROW_H = 1.9
CAPTION_H = 1.3


def block_height(seq):
    """Return the vertical space one sequence needs, in cell units.

    Includes the title, every row of residues and the gap before the next block.

    Args:
        seq: Sequence string.
    """
    n_rows = max(1, -(-len(seq) // LINE_LEN))
    return TOP_MARGIN + n_rows * ROW_H + BLOCK_GAP


def figure_width():
    """Return the width of the drawing in cell units (the same for every figure)."""
    return LEFT_MARGIN + LINE_LEN + (LINE_LEN // GROUP - 1) * GAP + 1.5


def legend_layout(rules, width_cells):
    """Place the legend entries, wrapping them to the figure width.

    Rules with an empty label get no entry.

    Args:
        rules: Rules to list.
        width_cells: Figure width in cell units.

    Returns:
        ``([(rule, x, row)], n_rows)``: the entries with their left x and legend
        row, and the number of legend rows used.
    """
    items, x, row = [], LEFT_MARGIN, 0
    for r in rules:
        if not r.label:
            continue
        w = 2.4 + 0.55 * len(r.label)
        if x > LEFT_MARGIN and x + w > width_cells - 1.0:
            x, row = LEFT_MARGIN, row + 1
        items.append((r, x, row))
        x += w + 1.5
    return items, (row + 1 if items else 0)


def visible_rules(records, rules):
    """Select the rules that deserve a legend entry.

    Args:
        records: ``[(seq_id, sequence)]`` shown in the figure.
        rules: All rules.

    Returns:
        The rules with a label that can act on at least one of the records.
    """
    return [r for r in rules if r.label and any(r.applies_to(seq) for _, seq in records)]


def circular_dna_present(records, circular):
    """Tell whether the figure contains DNA that is drawn as circular.

    Args:
        records: ``[(seq_id, sequence)]``.
        circular: The user's circular setting.

    Returns:
        True if ``circular`` is set and at least one record is DNA.
    """
    return bool(circular) and any(is_dna(seq) for _, seq in records)


def captions(rules, circular_dna=False):
    """Build the grey explanatory lines printed under the legend.

    Args:
        rules: The visible rules.
        circular_dna: Whether circular DNA is shown.

    Returns:
        A list of caption strings: a P1/P1' note for proteases, a recognition-site
        note for restriction enzymes and a note about the origin for circular DNA.
    """
    lines = []
    if any(r.kind == "protease" for r in rules):
        lines.append("Boxed residues: P1 and P1' - the bond between them is cleaved")
    if any(r.kind == "restriction" for r in rules):
        lines.append("Boxed: restriction recognition sites (either strand); "
                     "cut positions are listed in the CSV")
    if circular_dna:
        lines.append("Circular DNA: a site crossing the origin continues at the start "
                     "of the sequence")
    return lines


def legend_height(rules, width_cells, circular_dna=False):
    """Return the vertical space taken by the legend and its captions.

    Args:
        rules: The visible rules.
        width_cells: Figure width in cell units.
        circular_dna: Whether circular DNA is shown.

    Returns:
        The height in cell units; 0.0 if there is nothing to list.
    """
    items, n_rows = legend_layout(rules, width_cells)
    if not items:
        return 0.0
    return n_rows * LEGEND_ROW_H + 0.5 + CAPTION_H * len(captions(rules, circular_dna))


def layout_size(records, rules, circular=False):
    """Compute the size of the whole figure in cell units.

    Args:
        records: ``[(seq_id, sequence)]`` to draw one under another.
        rules: All rules (only visible ones add legend height).
        circular: The user's circular setting.

    Returns:
        ``(width, height)``.
    """
    w = figure_width()
    vis = visible_rules(records, rules)
    circ = circular_dna_present(records, circular)
    return w, sum(block_height(seq) for _, seq in records) + 0.5 + legend_height(vis, w, circ)


def safe_dpi(records, rules, want=300, cell_px=14, circular=False):
    """Lower the export resolution if the PNG would become too large.

    Matplotlib's Agg backend cannot write images beyond roughly 65 000 pixels per
    side, so the dpi is reduced so that neither side exceeds ``MAX_PX``.

    Args:
        records: ``[(seq_id, sequence)]`` to be drawn.
        rules: All rules.
        want: Requested dpi.
        cell_px: Size of one cell in points.
        circular: The user's circular setting.

    Returns:
        The dpi to use: ``want``, or less for very large figures (never below 20).
    """
    w, h = layout_size(records, rules, circular)
    return max(20, min(want, int(MAX_PX / (max(w, h) * cell_px / 72))))


def render(records, rules, dpi=100, cell_px=14, circular=False):
    """Draw sequences with their boxes, legend and captions.

    Every record is drawn one under another with the same rules. Protein and DNA
    records can be mixed: each rule is applied only to records it fits. The legend
    sits under the last row and lists every rule that acts on at least one record.

    Args:
        records: ``[(seq_id, sequence)]``.
        rules: Rules to apply and list in the legend.
        dpi: Resolution of the figure; also scales the pixel coordinates.
        cell_px: Size of one cell in points.
        circular: Treat DNA records as circular (proteins are unaffected).

    Returns:
        ``(figure, segments)``. ``figure`` is a matplotlib ``Figure`` that is not
        attached to any backend. ``segments`` are the dicts of ``build_segments``
        plus ``record_id``, ``topology`` (``"linear"``, ``"circular"`` or ``""``
        for proteins) and ``px``, the box ``(x0, y0, x1, y1)`` in pixels of an image
        saved at ``dpi`` with the origin at the top-left.
    """
    from matplotlib.figure import Figure

    width_cells, height_cells = layout_size(records, rules, circular)
    fig = Figure(figsize=(width_cells * cell_px / 72, height_cells * cell_px / 72),
                 dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width_cells)
    ax.set_ylim(height_cells, 0)  # y grows downward, like image pixels
    ax.axis("off")

    fs = cell_px * 0.8  # font size in points ~ cell width
    segs, y_off = [], 0.0
    for rec_id, seq in records:
        ax.text(0.3, y_off + 0.3, rec_id, ha="left", va="top", fontsize=fs * 1.1,
                fontweight="bold")
        for i, aa in enumerate(seq):
            y = y_off + TOP_MARGIN + cell_row(i) * ROW_H
            ax.text(cell_x(i) + 0.5, y + 0.75, aa, ha="center", va="center",
                    family="monospace", fontsize=fs, zorder=3)
            if (i + 1) % GROUP == 0:
                ax.text(cell_x(i) + 1, y - 0.25, str(i + 1), ha="right", va="bottom",
                        family="monospace", fontsize=fs * 0.75, color="#444444")
        topology = ("circular" if circular else "linear") if is_dna(seq) else ""
        for s in build_segments(seq, rules, y_off, circular):
            s["record_id"] = rec_id
            s["topology"] = topology
            segs.append(s)
        y_off += block_height(seq)

    for s in segs:
        rule = s["rule"]
        fill = rule.style == "fill"
        ax.add_patch(Rectangle((s["x0"], s["y0"]), s["x1"] - s["x0"], s["y1"] - s["y0"],
                               facecolor=rule.color if fill else "none",
                               edgecolor=rule.color if not fill else "none",
                               linewidth=1.6, zorder=1 if fill else 2))

    # legend under the last row: coloured box + protease name / motif
    vis = visible_rules(records, rules)
    items, n_rows = legend_layout(vis, width_cells)
    for r, x, row in items:
        y = y_off + row * LEGEND_ROW_H
        fill = r.style == "fill"
        ax.add_patch(Rectangle((x, y + 0.15), 1.6, 1.2,
                               facecolor=r.color if fill else "none",
                               edgecolor="none" if fill else r.color, linewidth=1.6))
        ax.text(x + 2.2, y + 0.75, r.label, ha="left", va="center", fontsize=fs * 0.95)
    circ = circular_dna_present(records, circular)
    for i, line in enumerate(captions(vis, circ) if items else []):
        ax.text(LEFT_MARGIN, y_off + n_rows * LEGEND_ROW_H + 0.1 + i * CAPTION_H, line,
                ha="left", va="top", fontsize=fs * 0.7, color="#555555")

    # pixel coordinates (origin top-left) for the CSV
    sx = fig.get_figwidth() * dpi / width_cells
    sy = fig.get_figheight() * dpi / height_cells
    for s in segs:
        s["px"] = (round(s["x0"] * sx, 1), round(s["y0"] * sy, 1),
                   round(s["x1"] * sx, 1), round(s["y1"] * sy, 1))
    return fig, segs


def safe_filename(rec_id, used):
    """Turn a SeqID into a safe, unique file name stem.

    Runs of characters other than letters, digits, ``.``, ``_`` and ``-`` become one
    underscore, leading and trailing underscores/dots are removed and the result is
    cut to 100 characters. Names that collide, ignoring case, get ``_2``, ``_3`` ...

    Args:
        rec_id: The sequence identifier.
        used: Set of lower-case stems already taken; updated in place.

    Returns:
        The stem, e.g. ``"gi|565165|gb|AAB31199.1|"`` becomes
        ``"gi_565165_gb_AAB31199.1"``; ``"sequence"`` if nothing is left.
    """
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", rec_id).strip("_.") or "sequence"
    stem = stem[:100]
    name, n = stem, 2
    while name.lower() in used:  # duplicate SeqIDs or names differing only by case
        name, n = f"{stem}_{n}", n + 1
    used.add(name.lower())
    return name


def save_one_png(record, rules, path, dpi=300, circular=False):
    """Save one sequence as its own PNG image.

    Args:
        record: ``(seq_id, sequence)``.
        rules: Rules to apply.
        path: Output file path.
        dpi: Requested resolution; lowered by ``safe_dpi`` if necessary.
        circular: Treat DNA as circular.

    Returns:
        The dpi actually used.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    dpi = safe_dpi([record], rules, dpi, circular=circular)
    fig, _ = render([record], rules, dpi=dpi, circular=circular)
    FigureCanvasAgg(fig)
    fig.savefig(path, dpi=dpi, facecolor="white")
    return dpi


def save_pngs(records, rules, outdir, dpi=300, only_matches=False, circular=False):
    """Save one PNG per sequence, named after its SeqID.

    Args:
        records: ``[(seq_id, sequence)]``.
        rules: Rules to apply.
        outdir: Output folder; created if missing.
        dpi: Requested resolution.
        only_matches: Skip sequences in which no rule produces a box.
        circular: Treat DNA as circular.

    Returns:
        The list of written file paths.
    """
    os.makedirs(outdir, exist_ok=True)
    used, written = set(), []
    for rec in records:
        if only_matches and not build_segments(rec[1], rules, circular=circular):
            continue
        path = os.path.join(outdir, safe_filename(rec[0], used) + ".png")
        save_one_png(rec, rules, path, dpi, circular)
        written.append(path)
    return written


def save_csv(path, records, rules, dpi=300, circular=False):
    """Write the box coordinates of all sequences to one CSV file.

    There is one row per drawn segment (a match split over rows or across the
    origin gives several rows). Pixel coordinates are measured inside each
    sequence's own PNG, as written by ``save_pngs`` (origin at the top-left).

    Columns: record_id, enzyme, enzyme_type, site_number, cut_after, subsequence,
    colour, style, match_start, match_end, segment_start, segment_end, row, x0_px,
    y0_px, x1_px, y1_px, topology.

    Args:
        path: Output file path.
        records: ``[(seq_id, sequence)]``.
        rules: Rules to apply.
        dpi: Export resolution the pixel coordinates refer to.
        circular: Treat DNA as circular.

    Returns:
        The number of rows written (excluding the header).
    """
    rows = []
    for rec in records:
        _, segs = render([rec], rules, dpi=safe_dpi([rec], rules, dpi, circular=circular),
                         circular=circular)
        rows += segs
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["record_id", "enzyme", "enzyme_type", "site_number", "cut_after",
                    "subsequence", "colour", "style", "match_start", "match_end",
                    "segment_start", "segment_end", "row",
                    "x0_px", "y0_px", "x1_px", "y1_px", "topology"])
        for s in rows:
            r = s["rule"]
            site = "" if s["site_number"] is None else s["site_number"]
            cut = "" if s["cut_after"] is None else s["cut_after"]
            w.writerow([s["record_id"], r.enzyme, r.kind if r.enzyme else "", site, cut,
                        s["residues"], r.color, r.style, s["start"], s["end"],
                        s["seg_start"], s["seg_end"], s["row"], *s["px"], s["topology"]])
    return len(rows)


def positive_int(text):
    """argparse type: a whole number >= 1.

    Raises:
        argparse.ArgumentTypeError: For anything else.
    """
    try:
        n = int(text)
    except ValueError:
        n = 0
    if n < 1:
        raise argparse.ArgumentTypeError("must be a whole number >= 1")
    return n


def site_numbers_arg(text):
    """argparse type wrapping ``parse_site_numbers``.

    Raises:
        argparse.ArgumentTypeError: If the text is not a valid site specification.
    """
    try:
        return parse_site_numbers(text)
    except ValueError as e:
        raise argparse.ArgumentTypeError(str(e))


def enzyme_arg(names, what):
    """Build an argparse type for ``NAME[=COLOUR[:box|fill]]`` enzyme options.

    Args:
        names: Valid canonical enzyme names.
        what: Word used in error messages, e.g. ``"protease"``.

    Returns:
        A parser function returning ``(canonical_name, hex_colour_or_None, style)``.
    """
    def parse(text):
        """Parse one ``NAME[=COLOUR[:style]]`` command-line value.

        Raises:
            argparse.ArgumentTypeError: For an unknown or ambiguous name, an unknown
                colour or a style other than ``box``/``fill``.
        """
        name, _, rest = text.partition("=")
        color, _, style = rest.partition(":")
        try:
            name = resolve_name(name, names, what)
        except ValueError as e:
            raise argparse.ArgumentTypeError(str(e))
        if color and not is_color_like(color):
            raise argparse.ArgumentTypeError(f"Unknown colour '{color}'")
        if style and style not in ("box", "fill"):
            raise argparse.ArgumentTypeError("Style must be 'box' or 'fill'")
        return name, (to_hex(color) if color else None), (style or "box")
    return parse


def parse_rule(text):
    """argparse type for ``--rule SUBSEQUENCE=COLOUR[:box|fill]``.

    Returns:
        A motif ``Rule`` with the colour converted to hex.

    Raises:
        argparse.ArgumentTypeError: If there is no ``=``, the colour is unknown or
            the style is not ``box``/``fill``.
    """
    if "=" not in text:
        raise argparse.ArgumentTypeError("Rule must look like SEQ=COLOUR[:box|fill]")
    sub, rest = text.split("=", 1)
    color, _, style = rest.partition(":")
    style = style or "box"
    if not is_color_like(color):
        raise argparse.ArgumentTypeError(f"Unknown colour '{color}'")
    if style not in ("box", "fill"):
        raise argparse.ArgumentTypeError("Style must be 'box' or 'fill'")
    return Rule(sub.strip(), to_hex(color), style)


# --------------------------------------------------------------------- GUI


ALL = "(all records)"


def run_gui(records, initial_rules, circular=False):
    """Open the interactive window.

    The left panel holds the record selector, the circular-DNA checkbox, the
    protease, restriction-enzyme and subsequence rows, the site-selection controls
    (all / first N / only site numbers), the "Show site counts" window, and the
    Render, Save PNG(s) and Save CSV buttons. The right panel shows the figure
    (scrollable). The call blocks until the window is closed.

    Args:
        records: ``[(seq_id, sequence)]``. With more than one record the selector
            offers "(all records)" first.
        initial_rules: Rules used to pre-fill the rows and the site-selection
            controls (taken from the command line).
        circular: Initial state of the circular-DNA checkbox.
    """
    import tkinter as tk
    from tkinter import colorchooser, filedialog, messagebox, ttk

    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    root = tk.Tk()
    root.title("seqboxes - sequence box visualiser")
    root.geometry("1250x800")

    left = ttk.Frame(root, padding=8)
    left.pack(side="left", fill="y")
    right = ttk.Frame(root)
    right.pack(side="right", fill="both", expand=True)

    ttk.Label(left, text="Record(s) - rules apply to every shown sequence").pack(anchor="w")
    multi = len(records) > 1
    names = ([ALL] if multi else []) + [r[0] for r in records]
    rec_var = tk.StringVar(value=names[0])
    rec_box = ttk.Combobox(left, textvariable=rec_var, values=names, width=42,
                           state="readonly")
    rec_box.pack(anchor="w", pady=(0, 8))

    circ_var = tk.BooleanVar(value=circular)
    ttk.Checkbutton(left, text="Circular DNA (sites may cross the origin)",
                    variable=circ_var).pack(anchor="w", pady=(0, 8))

    ttk.Label(left, text="Proteases (protein sequences)").pack(anchor="w")
    enz_frame = ttk.Frame(left)
    enz_frame.pack(anchor="w", fill="x")
    ttk.Button(left, text="+ Add protease",
               command=lambda: add_row(kind="protease")).pack(anchor="w", pady=4)

    ttk.Label(left, text="Restriction enzymes (DNA sequences)").pack(anchor="w")
    res_frame = ttk.Frame(left)
    res_frame.pack(anchor="w", fill="x")
    ttk.Button(left, text="+ Add restriction enzyme",
               command=lambda: add_row(kind="restriction")).pack(anchor="w", pady=4)

    # which of the protease's cleavage sites to box (numbered per sequence)
    first_enz = next((r for r in initial_rules if r.enzyme), None)
    mode_var = tk.StringVar(value="all")
    n_var = tk.StringVar(value="3")
    site_var = tk.StringVar(value="1")
    if first_enz and first_enz.only_sites:
        mode_var.set("only")
        site_var.set(format_site_numbers(first_enz.only_sites))
    elif first_enz and first_enz.max_sites:
        mode_var.set("first")
        n_var.set(str(first_enz.max_sites))
    sites_box = ttk.LabelFrame(left, text="Which cleavage sites to box", padding=6)
    sites_box.pack(anchor="w", fill="x", pady=4)
    ttk.Radiobutton(sites_box, text="All sites", variable=mode_var,
                    value="all").grid(row=0, column=0, columnspan=3, sticky="w")
    ttk.Radiobutton(sites_box, text="First N sites:", variable=mode_var,
                    value="first").grid(row=1, column=0, sticky="w")
    ttk.Spinbox(sites_box, from_=1, to=9999, textvariable=n_var,
                width=6).grid(row=1, column=1, sticky="w")
    ttk.Radiobutton(sites_box, text="Only site no.:", variable=mode_var,
                    value="only").grid(row=2, column=0, sticky="w")
    ttk.Entry(sites_box, textvariable=site_var, width=14).grid(row=2, column=1,
                                                                sticky="w")
    ttk.Label(sites_box, text="e.g. 3  or  1,4,7-9").grid(row=2, column=2, sticky="w",
                                                         padx=6)
    ttk.Label(sites_box, foreground="#555", wraplength=330,
              text="Sites are counted from the N-terminus (5' end for DNA); numbering "
                   "restarts for every sequence and every enzyme."
              ).grid(row=3, column=0, columnspan=3, sticky="w", pady=(4, 0))
    ttk.Button(sites_box, text="Show site counts",
               command=lambda: show_counts()).grid(row=4, column=0, columnspan=3,
                                                   sticky="w", pady=(4, 0))

    ttk.Label(left, text="Subsequences to box").pack(anchor="w")
    rows_frame = ttk.Frame(left)
    rows_frame.pack(anchor="w", fill="x")
    rows = []  # each: dict(frame, sub, color, style, kind)

    def pick_color(var, swatch):
        """Open the colour chooser and apply the choice to a row.

        Args:
            var: The row's colour ``StringVar``.
            swatch: The row's colour swatch label.
        """
        c = colorchooser.askcolor(color=var.get(), parent=root)[1]
        if c:
            var.set(c)
            swatch.config(bg=c)

    def add_row(rule=None, kind="motif"):
        """Add one input row to the left panel.

        A row has a name field (a drop-down for proteases, an editable drop-down for
        restriction enzymes, a text box for subsequences), a colour swatch with a hex
        entry, a box/fill selector and a remove button.

        Args:
            rule: Rule used to pre-fill the row; if None a new row is created with the
                next palette colour and, for enzymes, the first unused enzyme.
            kind: ``"protease"``, ``"restriction"`` or ``"motif"``.
        """
        if rule is None:
            rule = Rule("", PALETTE[len(rows) % len(PALETTE)])
            if kind != "motif":
                chosen = {r["sub"].get() for r in rows if r["kind"] == kind}
                options = (["trypsin"] + PROTEASES if kind == "protease" else
                           ["EcoRI", "BamHI", "HindIII"] + COMMERCIAL_NAMES)
                rule.enzyme = next((n for n in options if n not in chosen), options[0])
        fr = ttk.Frame({"protease": enz_frame, "restriction": res_frame,
                        "motif": rows_frame}[kind])
        fr.pack(anchor="w", pady=2)
        sub = tk.StringVar(value=rule.enzyme if kind != "motif" else rule.subseq)
        col = tk.StringVar(value=rule.color)
        sty = tk.StringVar(value=rule.style)
        if kind == "protease":
            ttk.Combobox(fr, textvariable=sub, values=PROTEASES, width=30,
                         state="readonly").pack(side="left")
        elif kind == "restriction":  # editable: any of the ~1000 REBASE names can be typed
            ttk.Combobox(fr, textvariable=sub, values=COMMERCIAL_NAMES,
                         width=30).pack(side="left")
        else:
            ttk.Entry(fr, textvariable=sub, width=14).pack(side="left")
        swatch = tk.Label(fr, bg=rule.color, width=3, relief="solid", cursor="hand2")
        swatch.pack(side="left", padx=4)
        swatch.bind("<Button-1>", lambda e: pick_color(col, swatch))
        ttk.Entry(fr, textvariable=col, width=9).pack(side="left")

        def sync(*_):
            """Keep the swatch colour in step with the colour entry."""
            if is_color_like(col.get()):
                swatch.config(bg=col.get())
        col.trace_add("write", sync)
        ttk.Combobox(fr, textvariable=sty, values=["box", "fill"], width=5,
                     state="readonly").pack(side="left", padx=4)
        row = dict(frame=fr, sub=sub, color=col, style=sty, kind=kind)

        def remove():
            """Remove this row from the window and from the row list."""
            fr.destroy()
            rows.remove(row)
        ttk.Button(fr, text="x", width=2, command=remove).pack(side="left")
        rows.append(row)

    for r in initial_rules:
        add_row(r, kind=r.kind)
    if not any(r.kind == "motif" for r in initial_rules):
        add_row(Rule("", PALETTE[len(rows) % len(PALETTE)]))  # one empty motif row
    ttk.Button(left, text="+ Add subsequence", command=add_row).pack(anchor="w", pady=6)

    status = tk.StringVar()
    ttk.Label(left, textvariable=status, wraplength=330, foreground="#555").pack(
        anchor="w", pady=6)

    # scrollable canvas for the figure
    canvas = tk.Canvas(right, bg="white")
    vbar = ttk.Scrollbar(right, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vbar.set)
    vbar.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)
    state = dict(widget=None, fig=None)

    def current():
        """Read the window state into sequences and rules.

        Enzyme names are resolved (case-insensitively), colours are validated and the
        site-selection controls are parsed. Proteases come first, then restriction
        enzymes, then subsequences.

        Returns:
            ``(selected_records, rules, problems)``. ``problems`` is a list of
            messages about invalid input (unknown colour or enzyme, bad N or site
            numbers); the offending rows are left out of ``rules``.
        """
        idx = rec_box.current()
        if multi:
            sel = records if idx <= 0 else [records[idx - 1]]
        else:
            sel = records
        enzymes, motifs, bad = [], [], []
        mode, max_sites, only_sites = mode_var.get(), None, None
        if mode == "first":
            try:
                max_sites = int(n_var.get())
                if max_sites < 1:
                    raise ValueError
            except ValueError:
                bad.append("'First N sites': N must be a whole number >= 1")
        elif mode == "only":
            try:
                only_sites = parse_site_numbers(site_var.get())
            except ValueError as e:
                bad.append(f"'Only site no.': {e}")
        for r in rows:
            name = r["sub"].get().strip()
            if not name:
                continue
            c = r["color"].get().strip()
            if not is_color_like(c):
                bad.append(f"Unknown colour '{c}'")
                continue
            if r["kind"] != "motif":
                names = PROTEASES if r["kind"] == "protease" else RESTRICTION_NAMES
                try:
                    name = resolve_name(name, names, r["kind"].replace("restriction",
                                                                      "restriction enzyme"))
                except ValueError as e:
                    bad.append(str(e))
                    continue
                enzymes.append(Rule("", to_hex(c), r["style"].get(), enzyme=name,
                                    kind=r["kind"], max_sites=max_sites,
                                    only_sites=only_sites))
            else:
                motifs.append(Rule(name.upper(), to_hex(c), r["style"].get()))
        enzymes.sort(key=lambda r: r.kind != "protease")  # proteases, then restriction
        return sel, enzymes + motifs, bad

    def draw():
        """Render the current selection into the window and update the status line.

        The status line reports, for each rule, how many sites or matches were found
        in how many applicable sequences, or that the rule was skipped because no
        sequence of the right type is shown. Invalid input is reported in a dialog
        instead.
        """
        sel, rules, bad = current()
        if bad:
            messagebox.showwarning("Check input", "\n".join(bad))
            return
        fig, segs = render(sel, rules, circular=circ_var.get())
        if state["widget"] is not None:
            state["widget"].destroy()
        fc = FigureCanvasTkAgg(fig, master=canvas)
        w = fc.get_tk_widget()
        canvas.delete("all")
        canvas.create_window(0, 0, anchor="nw", window=w)
        fc.draw()
        canvas.configure(scrollregion=(0, 0, w.winfo_reqwidth(), w.winfo_reqheight()))
        state.update(widget=w, fig=fig)
        hits = {}  # rule index -> set of (record, match start)
        for sg in segs:
            hits.setdefault(sg["rule_index"], set()).add((sg["record_id"], sg["start"]))
        parts = []
        for i, r in enumerate(rules):
            usable = sum(1 for _, sq in sel if r.applies_to(sq))
            if usable == 0:
                kind = "protein" if r.kind == "protease" else "DNA"
                parts.append(f"{r.label}: skipped - no {kind} sequences shown")
                continue
            h = hits.get(i, set())
            n_rec = len({rid for rid, _ in h})
            what = "site(s)" if r.enzyme else "match(es)"
            parts.append(f"{r.label}: {len(h)} {what} in {n_rec}/{usable} sequence(s)")
        status.set(f"{len(sel)} sequence(s) shown. " + "; ".join(parts))

    def export_png():
        """Save the figure(s) as PNG.

        With one sequence selected a file dialog asks for a file name; with several, a
        folder is requested and one ``<SeqID>.png`` is written per sequence (honouring
        the "Skip sequences without matches" checkbox).
        """
        sel, rules, bad = current()
        if bad:
            messagebox.showwarning("Check input", "\n".join(bad))
            return
        if len(sel) == 1:
            p = filedialog.asksaveasfilename(
                defaultextension=".png", filetypes=[("PNG image", "*.png")],
                initialfile=safe_filename(sel[0][0], set()) + ".png")
            if p:
                save_one_png(sel[0], rules, p, circular=circ_var.get())
                status.set(f"Saved {p}")
            return
        d = filedialog.askdirectory(title="Folder for one PNG per sequence")
        if d:
            paths = save_pngs(sel, rules, d, only_matches=only_var.get(),
                              circular=circ_var.get())
            status.set(f"Saved {len(paths)} PNG file(s) to {d}")

    def export_csv():
        """Ask for a file name and write the box-coordinates CSV."""
        sel, rules, bad = current()
        if bad:
            messagebox.showwarning("Check input", "\n".join(bad))
            return
        p = filedialog.asksaveasfilename(defaultextension=".csv",
                                         filetypes=[("CSV file", "*.csv")])
        if p:
            n = save_csv(p, sel, rules, circular=circ_var.get())
            status.set(f"Saved {n} box segment(s) to {p}")

    def show_counts():
        """Open a window listing the number of sites per sequence and enzyme.

        For every sequence and every chosen enzyme it shows the total number of sites
        and ``#k@pos`` for each site: site number k and the position of its first
        residue or base. Numbering restarts in every sequence and for every enzyme, so
        the list tells which numbers can be used in "Only site no.". Enzymes that fit
        none of the shown sequences are listed at the end.
        """
        sel, rules, _ = current()
        enz = list({(r.kind, r.enzyme): r for r in rules if r.enzyme}.values())
        if not enz:
            messagebox.showinfo("Site counts", "Add at least one enzyme first.")
            return
        win = tk.Toplevel(root)
        win.title("Sites per sequence")
        txt = tk.Text(win, width=100, height=30, wrap="word", font=("Courier", 10))
        bar = ttk.Scrollbar(win, orient="vertical", command=txt.yview)
        txt.configure(yscrollcommand=bar.set)
        bar.pack(side="right", fill="y")
        txt.pack(side="left", fill="both", expand=True)
        txt.insert("end", "#k@pos = site number k at position pos (first residue / base "
                          "of the boxed span).\nNumbering restarts in every sequence "
                          "and for every enzyme.\n")
        last = None
        for rid, r, pos in site_count_table(sel, enz, circ_var.get()):
            if rid != last:
                txt.insert("end", f"\n{rid}\n")
                last = rid
            txt.insert("end", f"  {r.enzyme} ({r.kind}): {len(pos)} site(s)\n")
            if pos:
                txt.insert("end", "    " + ", ".join(
                    f"#{k}@{p}" for k, p in enumerate(pos, 1)) + "\n")
        unused = [r.enzyme for r in enz if not any(r.applies_to(sq) for _, sq in sel)]
        if unused:
            txt.insert("end", "\nNot applicable to any shown sequence (wrong sequence "
                              "type): " + ", ".join(unused) + "\n")
        txt.configure(state="disabled")

    only_var = tk.BooleanVar(value=False)
    if multi:
        ttk.Checkbutton(left, text="Skip sequences without matches when saving PNGs",
                        variable=only_var).pack(anchor="w", pady=(0, 4))
    btns = ttk.Frame(left)
    btns.pack(anchor="w")
    ttk.Button(btns, text="Render", command=draw).pack(side="left", padx=2)
    ttk.Button(btns, text="Save PNG(s)", command=export_png).pack(side="left", padx=2)
    ttk.Button(btns, text="Save CSV", command=export_csv).pack(side="left", padx=2)
    rec_box.bind("<<ComboboxSelected>>", lambda e: draw())

    draw()
    root.mainloop()


# --------------------------------------------------------------------- CLI


def main(argv=None):
    """Command-line entry point.

    Without a file argument the program asks for the path in the terminal. It then
    opens the GUI, unless ``--png-dir``/``--csv`` (headless export),
    ``--count-sites``, ``--list-enzymes`` or ``--list-restriction`` are given, in
    which case it prints or writes the results and returns without a window.

    Args:
        argv: Argument list; None means ``sys.argv[1:]``. Useful in tests.

    Raises:
        SystemExit: On invalid arguments (argparse), an unreadable file or an
            unknown ``--record``.
    """
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", nargs="?", help="FASTA or FASTQ file (prompted if omitted)")
    ap.add_argument("--record", action="append",
                    help="only use this record (SeqID or prefix); repeatable. "
                         "Default: all records in the file")
    ap.add_argument("--enzyme", action="append", default=[], metavar="NAME[=COLOUR[:box|fill]]",
                    type=enzyme_arg(PROTEASES, "protease"),
                    help="box the cleavage sites of a protease in protein sequences "
                         "(Expasy rules via Pyteomics); repeatable. Case-insensitive, "
                         "a unique prefix is enough")
    ap.add_argument("--restriction", action="append", default=[],
                    metavar="NAME[=COLOUR[:box|fill]]",
                    type=enzyme_arg(RESTRICTION_NAMES, "restriction enzyme"),
                    help="box the recognition sites of a restriction enzyme in DNA "
                         "sequences (Bio.Restriction); repeatable")
    ap.add_argument("--list-enzymes", action="store_true",
                    help="print the available proteases and exit")
    ap.add_argument("--list-restriction", action="store_true",
                    help="print the available restriction enzymes and exit")
    grp = ap.add_mutually_exclusive_group()
    grp.add_argument("--max-sites", type=positive_int, metavar="N",
                     help="box only the first N sites of each enzyme in each "
                          "sequence (counted from the N-terminus / 5' end)")
    grp.add_argument("--site", type=site_numbers_arg, metavar="NUMBERS",
                     help="box only these site numbers, e.g. 3 or 1,4,7-9. "
                          "Numbering restarts for every sequence and enzyme")
    ap.add_argument("--circular", action="store_true",
                    help="treat DNA records as circular: sites crossing the origin are "
                         "found and drawn split between the last and the first row")
    ap.add_argument("--count-sites", action="store_true",
                    help="print how many sites each enzyme has in each sequence "
                         "(and their positions), then exit")
    ap.add_argument("--rule", action="append", type=parse_rule, default=[],
                    help="SEQ=COLOUR[:box|fill]; repeatable")
    ap.add_argument("--png-dir", help="headless: write one <SeqID>.png per sequence "
                                      "into this folder")
    ap.add_argument("--only-matches", action="store_true",
                    help="with --png-dir: skip sequences that have no matches")
    ap.add_argument("--csv", help="headless: write box-coordinates CSV here")
    a = ap.parse_args(argv)

    if a.list_enzymes:
        print("\n".join(PROTEASES))
        return
    if a.list_restriction:
        print("\n".join(RESTRICTION_NAMES))
        return
    enzyme_rules = [Rule("", col or PALETTE[i % len(PALETTE)], style, enzyme=name,
                         kind=kind, max_sites=a.max_sites, only_sites=a.site)
                    for i, (kind, (name, col, style)) in enumerate(
                        [("protease", t) for t in a.enzyme]
                        + [("restriction", t) for t in a.restriction])]
    if (a.max_sites or a.site or a.count_sites) and not enzyme_rules:
        ap.error("--max-sites, --site and --count-sites need at least one "
                 "--enzyme or --restriction")
    a.rule = enzyme_rules + a.rule

    path = a.file
    if not path:
        if a.png_dir or a.csv:
            ap.error("a FASTA/FASTQ file is required in headless mode")
        while not path:
            path = input("Path to FASTA/FASTQ file: ").strip().strip("\"'")
            path = os.path.expanduser(path)
            if path and not os.path.isfile(path):
                print(f"File not found: {path}")
                path = None
    try:
        records = read_records(path)
    except (OSError, ValueError) as e:
        sys.exit(f"Error: {e}")
    if a.record:
        selected = []
        for q in a.record:
            hits = [r for r in records if r[0].startswith(q) or q in r[0].split("|")]
            if not hits:
                sys.exit(f"Record '{q}' not found. Available: "
                         + ", ".join(r[0][:30] for r in records))
            selected += [h for h in hits if h not in selected]
        records = selected

    if a.circular and not any(is_dna(sq) for _, sq in records):
        print("Note: --circular has no effect - no DNA sequences among the selected "
              "records.", file=sys.stderr)
    for r in enzyme_rules:  # tell the user when an enzyme cannot act on the records
        if not any(r.applies_to(sq) for _, sq in records):
            kind = "protein" if r.kind == "protease" else "DNA"
            print(f"Note: {r.enzyme} ({r.kind}) - no {kind} sequences among the "
                  "selected records, nothing will be boxed.", file=sys.stderr)

    if a.count_sites:
        print("record_id\tenzyme_type\tenzyme\tn_sites\t"
              "positions (k-th = site number k; first residue/base of the boxed span)")
        for rid, r, pos in site_count_table(records, enzyme_rules, a.circular):
            print(f"{rid}\t{r.kind}\t{r.enzyme}\t{len(pos)}\t{', '.join(map(str, pos))}")
        return

    if a.png_dir or a.csv:  # headless mode
        matplotlib.use("Agg")
        if a.png_dir:
            paths = save_pngs(records, a.rule, a.png_dir, only_matches=a.only_matches,
                              circular=a.circular)
            print(f"Wrote {len(paths)} PNG file(s) to {a.png_dir}")
        if a.csv:
            n = save_csv(a.csv, records, a.rule, circular=a.circular)
            print(f"Wrote {a.csv} ({n} segments)")
        return
    matplotlib.use("TkAgg")
    run_gui(records, a.rule, a.circular)


if __name__ == "__main__":
    main()
