#!/usr/bin/env python3
"""
Summarize mismatch statistics across SAM files.

High-level flow:
1. Load CDS bounds from a GTF.
2. For each aligned read, keep it only if shifted 5' position (with user offset)
   falls inside the corresponding CDS.
3. Parse MD+CIGAR to recover mismatch events.
4. Aggregate per-file totals and mutation spectra, then write one Excel sheet.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import zipfile
from collections import Counter
from typing import Dict, List, Tuple
from xml.sax.saxutils import escape


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
MD_RE = re.compile(r"(\d+|\^[A-Z]+|[A-Z])")
BASES = ("A", "C", "G", "T")
REF_SPAN_OPS = ("M", "D", "N", "=", "X")


def parse_tags(fields: List[str]) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    for tag_field in fields[11:]:
        parts = tag_field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_gtf_attributes(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_text.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        key_value = part.split(" ", 1)
        if len(key_value) != 2:
            continue
        key, value = key_value
        attrs[key] = value.strip().strip('"')
    return attrs


def normalize_ref_name(name: str) -> str:
    base = name.strip().split("|", 1)[0]
    if base.endswith("_transcript"):
        return base[: -len("_transcript")]
    return base


def load_cds_ranges(gtf_path: str) -> Dict[str, Tuple[int, int]]:
    """Build transcript/gene -> [min CDS start, max CDS end] from GTF CDS rows."""
    cds_ranges: Dict[str, List[int]] = {}
    with open(gtf_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue

            start = int(fields[3])
            end = int(fields[4])
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")

            names = set()
            if gene_id:
                names.add(normalize_ref_name(gene_id))
            if transcript_id:
                names.add(normalize_ref_name(transcript_id))

            for name in names:
                if not name:
                    continue
                if name not in cds_ranges:
                    cds_ranges[name] = [start, end]
                else:
                    cds_ranges[name][0] = min(cds_ranges[name][0], start)
                    cds_ranges[name][1] = max(cds_ranges[name][1], end)

    return {k: (v[0], v[1]) for k, v in cds_ranges.items()}


def reference_span(cigar: str) -> int:
    span = 0
    for length_str, op in CIGAR_RE.findall(cigar):
        if op in REF_SPAN_OPS:
            span += int(length_str)
    return span


def shifted_five_prime_position(pos: int, cigar: str, is_reverse: bool, offset: int) -> int:
    """
    Compute the offset-adjusted 5' genomic position for CDS gating.

    Forward read:  shifted = 5prime + offset
    Reverse read:  shifted = 5prime - offset
    """
    span = reference_span(cigar)
    if span <= 0:
        return -1
    five_prime = (pos + span - 1) if is_reverse else pos
    return (five_prime - offset) if is_reverse else (five_prime + offset)


def build_ref_to_read_map(cigar: str) -> Dict[int, int]:
    """
    Map reference-offset positions to read sequence indices for aligned positions.
    Offsets are 0-based from the alignment start.
    """
    ref_to_read: Dict[int, int] = {}
    read_pos = 0
    ref_pos = 0

    for length_str, op in CIGAR_RE.findall(cigar):
        length = int(length_str)
        if op in ("M", "=", "X"):
            for i in range(length):
                ref_to_read[ref_pos + i] = read_pos + i
            read_pos += length
            ref_pos += length
        elif op in ("I", "S"):
            read_pos += length
        elif op in ("D", "N"):
            ref_pos += length
        elif op in ("H", "P"):
            continue

    return ref_to_read


def extract_mismatches(seq: str, cigar: str, md: str) -> List[Tuple[str, str]]:
    """
    Return list of (ref_base, read_base) mismatches using MD + CIGAR + SEQ.
    """
    ref_to_read = build_ref_to_read_map(cigar)
    mismatches: List[Tuple[str, str]] = []
    ref_offset = 0
    seq = seq.upper()

    for token in MD_RE.findall(md.upper()):
        if token[0].isdigit():
            ref_offset += int(token)
        elif token.startswith("^"):
            ref_offset += len(token) - 1
        else:
            ref_base = token
            read_idx = ref_to_read.get(ref_offset)
            if read_idx is not None and 0 <= read_idx < len(seq):
                read_base = seq[read_idx]
                if ref_base in BASES and read_base in BASES and ref_base != read_base:
                    mismatches.append((ref_base, read_base))
            ref_offset += 1

    return mismatches


def analyze_sam(path: str, cds_ranges: Dict[str, Tuple[int, int]], offset: int) -> Dict[str, object]:
    """Return all per-file counters used in the output summary row."""
    total_aligned_reads = 0
    aligned_reads_in_cds = 0
    aligned_reads_with_mismatch = 0
    mismatch_events = Counter()
    aligned_reads_missing_md = 0
    aligned_reads_missing_cds = 0
    aligned_reads_outside_cds = 0

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            flag = int(fields[1])
            rname = normalize_ref_name(fields[2])
            pos = int(fields[3])
            cigar = fields[5]
            seq = fields[9]

            # Count only aligned reads as denominator.
            if (flag & 4) != 0 or cigar == "*":
                continue
            total_aligned_reads += 1
            is_reverse = (flag & 16) != 0

            # CDS inclusion filter: only analyze reads whose shifted 5' position is in CDS.
            cds = cds_ranges.get(rname)
            if cds is None:
                aligned_reads_missing_cds += 1
                continue

            shifted_pos = shifted_five_prime_position(pos, cigar, is_reverse, offset)
            if shifted_pos < cds[0] or shifted_pos > cds[1]:
                aligned_reads_outside_cds += 1
                continue
            aligned_reads_in_cds += 1

            tags = parse_tags(fields)
            md = tags.get("MD")

            if not md or cigar == "*" or seq == "*":
                aligned_reads_missing_md += 1
                continue

            # Count mismatch events after passing alignment + CDS + MD checks.
            mismatches = extract_mismatches(seq, cigar, md)
            if mismatches:
                aligned_reads_with_mismatch += 1
                for from_base, to_base in mismatches:
                    mismatch_events[(from_base, to_base)] += 1

    return {
        "sam_file": os.path.basename(path),
        "total_aligned_reads": total_aligned_reads,
        "aligned_reads_in_cds": aligned_reads_in_cds,
        "aligned_reads_with_mismatch": aligned_reads_with_mismatch,
        "aligned_reads_missing_md": aligned_reads_missing_md,
        "aligned_reads_missing_cds": aligned_reads_missing_cds,
        "aligned_reads_outside_cds": aligned_reads_outside_cds,
        "mismatch_events": dict(mismatch_events),
    }


def excel_col_name(index_1_based: int) -> str:
    result = ""
    index = index_1_based
    while index > 0:
        index, rem = divmod(index - 1, 26)
        result = chr(ord("A") + rem) + result
    return result


def xml_cell(value: object) -> str:
    if isinstance(value, (int, float)):
        return f'<c t="n"><v>{value}</v></c>'
    return f'<c t="inlineStr"><is><t>{escape(str(value))}</t></is></c>'


def write_xlsx(output_path: str, headers: List[str], rows: List[List[object]]) -> None:
    """Write a minimal single-sheet XLSX (no external Excel deps)."""
    last_col = excel_col_name(len(headers))
    last_row = len(rows) + 1
    sheet_rows = []
    sheet_rows.append(
        "<row r=\"1\">"
        + "".join(xml_cell(h) for h in headers)
        + "</row>"
    )
    for i, row in enumerate(rows, start=2):
        sheet_rows.append(
            f"<row r=\"{i}\">" + "".join(xml_cell(v) for v in row) + "</row>"
        )
    sheet_data = "".join(sheet_rows)

    worksheet_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<dimension ref="A1:{last_col}{last_row}"/>'
        f"<sheetData>{sheet_data}</sheetData>"
        "</worksheet>"
    )

    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets><sheet name="MismatchSummary" sheetId="1" r:id="rId1"/></sheets>'
        "</workbook>"
    )

    content_types_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        "</Types>"
    )

    root_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        "</Relationships>"
    )

    workbook_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet1.xml"/>'
        "</Relationships>"
    )

    with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types_xml)
        zf.writestr("_rels/.rels", root_rels_xml)
        zf.writestr("xl/workbook.xml", workbook_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", workbook_rels_xml)
        zf.writestr("xl/worksheets/sheet1.xml", worksheet_xml)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze mismatch statistics for all SAM files in a directory and write an Excel file."
    )
    parser.add_argument(
        "-d",
        "--directory",
        default=".",
        help="Directory containing SAM files (default: current directory)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="sam_mismatch_summary_all.xlsx",
        help="Output Excel filename (default: sam_mismatch_summary_all.xlsx)",
    )
    parser.add_argument(
        "--cds-gtf",
        default=None,
        help="GTF file containing CDS features (prompted if omitted)",
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=None,
        help="Offset from read 5' position used for CDS inclusion test (prompted if omitted)",
    )
    args = parser.parse_args()

    sam_paths = sorted(
        glob.glob(os.path.join(args.directory, "*.SAM"))
        + glob.glob(os.path.join(args.directory, "*.sam"))
    )
    if not sam_paths:
        raise SystemExit(f"No SAM files found in: {os.path.abspath(args.directory)}")
    cds_gtf = args.cds_gtf
    if cds_gtf is None:
        cds_gtf = input("Enter path to CDS GTF file: ").strip()
    if not cds_gtf:
        raise SystemExit("No CDS GTF path provided.")
    if not os.path.exists(cds_gtf):
        raise SystemExit(f"CDS GTF file not found: {os.path.abspath(cds_gtf)}")

    offset = args.offset
    if offset is None:
        offset = int(input("Enter offset for CDS inclusion filtering: ").strip())

    cds_ranges = load_cds_ranges(cds_gtf)
    if not cds_ranges:
        raise SystemExit(f"No CDS records found in: {os.path.abspath(cds_gtf)}")

    mutation_pairs = [
        (from_base, to_base)
        for from_base in BASES
        for to_base in BASES
        if from_base != to_base
    ]
    from_cols = [f"mismatches_from_{base}" for base in BASES]
    from_ratio_cols = [f"percent_from_{base}_vs_total_cds_aligned" for base in BASES]
    transition_ratio_cols = [f"percent_{f}_to_{t}_from_{f}" for f, t in mutation_pairs]

    headers = [
        "sam_file",
        "total_aligned_reads",
        "total_cds_aligned_reads",
        "aligned_reads_with_mismatch",
        "aligned_reads_missing_md",
        "aligned_reads_missing_cds",
        "aligned_reads_outside_cds",
        "total_mismatch_events",
    ] + from_cols + from_ratio_cols + [f"{f}_to_{t}" for f, t in mutation_pairs] + transition_ratio_cols

    # Each SAM contributes one row with counts and percentages.
    rows: List[List[object]] = []
    for sam_path in sam_paths:
        result = analyze_sam(sam_path, cds_ranges, offset)
        events = result["mismatch_events"]
        mismatches_from = {
            base: sum(events.get((base, to), 0) for to in BASES if to != base)
            for base in BASES
        }
        total_aligned = result["aligned_reads_in_cds"]
        total_mismatch_events = sum(mismatches_from.values())

        from_percents = [
            ((mismatches_from[base] / total_aligned) * 100.0) if total_aligned else 0.0
            for base in BASES
        ]
        transition_percents = []
        for from_base, to_base in mutation_pairs:
            from_total = mismatches_from[from_base]
            percent = (
                (events.get((from_base, to_base), 0) / from_total) * 100.0
                if from_total
                else 0.0
            )
            transition_percents.append(percent)

        row: List[object] = [
            result["sam_file"],
            result["total_aligned_reads"],
            result["aligned_reads_in_cds"],
            result["aligned_reads_with_mismatch"],
            result["aligned_reads_missing_md"],
            result["aligned_reads_missing_cds"],
            result["aligned_reads_outside_cds"],
            total_mismatch_events,
        ] + [mismatches_from[base] for base in BASES] + from_percents + [
            events.get((f, t), 0) for f, t in mutation_pairs
        ] + transition_percents
        rows.append(row)

    write_xlsx(args.output, headers, rows)
    print(f"Wrote {len(rows)} file summaries to: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
