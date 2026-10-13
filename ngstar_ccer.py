#!/usr/bin/env python3
"""NG-STAR minimum-spanning-tree-based clonal-complex assignment.

Implements the *custom CC designation rules* used in the Euro-GASP NG-STAR
workflow. It is not a general-purpose PHYLOViZ CC-grouping replacement.

Version 1.1.0: retains the original founder/CC and MST tie-breaking algorithms
for already-cleaned inputs, with validation, input conveniences, and audit tools.

No third-party Python dependencies are required.
"""

import argparse
import csv
import hashlib
import json
import platform
import re
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

VERSION = "1.1.0"
LOCI = ["penA", "mtrR", "porB", "ponA", "gyrA", "parC", "23S"]
UNGROUPED = "Ungrouped"
ST_ALIASES = {"st", "sequence type", "sequence_type", "ng-star st", "ng-star type",
              "ng_star_type", "ngstar_type"}
ISOLATE_ST_ALIASES = ("ST", "NG-STAR TYPE", "NG_STAR_TYPE", "NGSTAR_TYPE")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class UF:
    def __init__(self, n):
        self._id = list(range(n))
        self._rank = [0] * n

    def find(self, p):
        while p != self._id[p]:
            self._id[p] = self._id[self._id[p]]
            p = self._id[p]
        return p

    def connected(self, p, q):
        return self.find(p) == self.find(q)

    def union(self, p, q):
        i, j = self.find(p), self.find(q)
        if i == j:
            return
        if self._rank[i] < self._rank[j]:
            self._id[i] = j
        elif self._rank[i] > self._rank[j]:
            self._id[j] = i
        else:
            self._id[j] = i
            self._rank[i] += 1


def hamm(v1, v2):
    if len(v1) != len(v2):
        raise ValueError("Allele profiles have different numbers of loci")
    return sum(a != b for a, b in zip(v1, v2))


def calc_lvs(profiles):
    """Count SLVs, DLVs, etc. This is O(number_of_STs**2)."""
    n = len(profiles)
    maxlen = len(profiles[0])
    lvs = [[0] * maxlen for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = hamm(profiles[i], profiles[j])
            if 0 < d <= maxlen:
                lvs[i][d - 1] += 1
                lvs[j][d - 1] += 1
    return lvs, maxlen


def goeburst_mst(profiles, st_order):
    """Construct full MST with goeBURST-style ranking.

    Keeps the original ranking: LV counts, profile frequency and numerical ST;
    Kruskal's algorithm processes distances 1..7. Note that profile frequency
    in a one-row-per-ST database is *not* isolate observation frequency.

    Pairwise comparisons are repeated for each level to limit peak edge memory.
    This is CPU-intensive for thousands of STs. The graph depends on the full
    reference ST collection (not on a single study's subset).
    """
    n = len(profiles)
    if n <= 1:
        return []
    maxlen = len(profiles[0])
    freq_map = Counter(profiles)
    freq = [freq_map[profiles[i]] for i in range(n)]

    def st_num(st):
        try:
            return int(st)
        except ValueError:
            return float("inf")

    st_ids = [st_num(st) for st in st_order]
    print("  Computing LV counts ...", flush=True)
    lvs, _ = calc_lvs(profiles)
    uf = UF(n)
    tree = []

    for level in range(1, maxlen + 1):
        if len(tree) == n - 1:
            break
        print(f"  Processing level {level} edges (diff={level}) ...", flush=True)
        level_edges = [(i, j) for i in range(n) for j in range(i + 1, n)
                       if hamm(profiles[i], profiles[j]) == level]

        def edge_key(e):
            u, v = e
            key = []
            for k in range(maxlen):
                key.append(-max(lvs[u][k], lvs[v][k]))
                key.append(-min(lvs[u][k], lvs[v][k]))
            key.append(-max(freq[u], freq[v]))
            key.append(-min(freq[u], freq[v]))
            key.append(min(st_ids[u], st_ids[v]))
            key.append(max(st_ids[u], st_ids[v]))
            return key

        level_edges.sort(key=edge_key)
        for u, v in level_edges:
            if len(tree) == n - 1:
                break
            if not uf.connected(u, v):
                uf.union(u, v)
                tree.append((u, v, level))
    return tree


def detect_delimiter(filepath):
    """Prefer delimiter consistent with known column names in the header."""
    with open(filepath, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
    header = sample.splitlines()[0] if sample else ""
    candidates = (",", ";", "\t")
    scores = {}
    for delim in candidates:
        cells = [v.strip().lower() for v in next(csv.reader([header], delimiter=delim), [])]
        scores[delim] = (sum(v in {"id", "st", "cc"} or v in {x.lower() for x in LOCI}
                             or v in ST_ALIASES for v in cells), len(cells))
    chosen = max(candidates, key=lambda d: scores[d])
    if scores[chosen][1] > 1 and scores[chosen][0] > 0:
        return chosen
    try:
        return csv.Sniffer().sniff(sample, delimiters=",;\t").delimiter
    except csv.Error:
        # Header-only files with a nonstandard first-column name still work.
        counts = {d: header.count(d) for d in candidates}
        chosen = max(candidates, key=counts.get)
        if counts[chosen] == 0:
            raise ValueError(f"Could not detect comma, semicolon or tab delimiter: {filepath}")
        return chosen


def _find_column(fieldnames, names):
    return next((c for c in fieldnames if c.strip().lower() in names), None)


def load_alleles(filepath, first_column_is_st=False):
    """Read a 7-locus reference; require one *complete and unique* row per ST."""
    delimiter = detect_delimiter(filepath)
    st_alleles = {}
    st_order = []
    normalized_count = 0
    profiles_seen = {}
    with open(filepath, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        fieldnames = reader.fieldnames or []
        st_col = _find_column(fieldnames, ST_ALIASES)
        if st_col is None and first_column_is_st and fieldnames:
            st_col = fieldnames[0]
        if st_col is None:
            raise ValueError("ST column not recognized. Rename it 'ST' or pass "
                             "--first-column-is-st to use the first column.")
        columns = {loc: _find_column(fieldnames, {loc.lower()}) for loc in LOCI}
        missing = [loc for loc, col in columns.items() if col is None]
        if missing:
            raise ValueError(f"Reference file is missing NG-STAR loci: {', '.join(missing)}")

        for lineno, row in enumerate(reader, start=2):
            if None in row:  # extra fields, usually an incorrectly parsed delimiter
                raise ValueError(f"Too many fields in reference at record {lineno}; "
                                 "check CSV delimiter and quoting")
            st = (row.get(st_col) or "").strip()
            if not st:
                raise ValueError(f"Blank ST at reference record {lineno}")
            if st in st_alleles:
                raise ValueError(f"Duplicate NG-STAR ST {st!r} at record {lineno}")
            alleles = []
            for loc in LOCI:
                value = (row.get(columns[loc]) or "").strip()
                if not value or value.casefold() in ("na", "nan", "null", "none", "-"):
                    raise ValueError(f"Missing {loc} allele for ST {st} (record {lineno})")
                # European decimal comma in a semicolon-separated NG-STAR export.
                if re.fullmatch(r"\d+,\d+", value):
                    value = value.replace(",", ".")
                    normalized_count += 1
                alleles.append(value)
            profile = tuple(alleles)
            if profile in profiles_seen:
                other = profiles_seen[profile]
                raise ValueError(f"STs {other} and {st} share an identical 7-locus profile. "
                                 "A one-row-per-unique-ST reference is required; "
                                 "otherwise zero-distance edges are omitted by the legacy MST.")
            profiles_seen[profile] = st
            st_alleles[st] = profile
            st_order.append(st)
    if not st_order:
        raise ValueError("Reference contains no ST profiles")
    if normalized_count:
        print(f"  Normalised {normalized_count} decimal-comma allele values to dots.")
    return st_alleles, st_order


def load_json(filepath):
    with open(filepath, encoding="utf-8-sig") as fh:
        content = fh.read()
    try:
        return json.loads(content)
    except json.JSONDecodeError:
        # Compatibility with older PHYLOViZ export files containing trailing commas.
        cleaned = re.sub(r",\s*([\]}])", r"\1", content)
        if cleaned == content:
            raise
        print("  Warning: removed trailing commas from legacy PHYLOViZ JSON.")
        return json.loads(cleaned)


def write_output(filepath, st_order, st_alleles, cc_assignments):
    with open(filepath, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ST"] + LOCI + ["CC"])
        for st in st_order:
            writer.writerow([st] + list(st_alleles[st]) + [cc_assignments.get(st, UNGROUPED)])


def write_counts(filepath, st_order, st_to_deg):
    with open(filepath, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ST", "LinkCount"])
        for st in st_order:
            writer.writerow([st, st_to_deg.get(st, 0)])


def assign_ccs_to_isolates(filepath, cc_assignments, output_path):
    """Exact ST lookup; missing/unassigned statuses are not named CCs."""
    delimiter = detect_delimiter(filepath)
    stats = {"total": 0, "assigned": 0, "ungrouped": 0, "not_found": 0,
             "missing_st": 0}
    with open(filepath, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        fieldnames = reader.fieldnames or []
        id_col = _find_column(fieldnames, {"id"})
        if id_col is None:
            raise ValueError(f"Isolate file requires 'ID' column; found {fieldnames}")
        st_col = _find_column(fieldnames, {s.lower() for s in ISOLATE_ST_ALIASES})
        if st_col is None:
            raise ValueError(f"Isolate file requires ST or NG-STAR TYPE column; "
                             f"found {fieldnames}")
        with open(output_path, "w", newline="", encoding="utf-8") as out:
            writer = csv.writer(out)
            writer.writerow(["ID", "ST", "CC"])
            for lineno, row in enumerate(reader, start=2):
                if None in row:
                    raise ValueError(f"Too many fields in isolate file at record {lineno}")
                stats["total"] += 1
                isolate_id = (row.get(id_col) or "").strip()
                st = (row.get(st_col) or "").strip()
                if not st:
                    cc = "Missing_ST"
                    stats["missing_st"] += 1
                elif st not in cc_assignments:
                    cc = "ST_not_found"
                    stats["not_found"] += 1
                else:
                    cc = cc_assignments[st]
                    stats["ungrouped" if cc == UNGROUPED else "assigned"] += 1
                writer.writerow([isolate_id, st, cc])
    return stats


def allele_diff(a1, a2):
    if not a1 or not a2:
        return 99
    return hamm(a1, a2)


def assign_ccs_from_adj(st_alleles, st_order, adj, degree, id_to_st):
    """Original/custom founder designation logic (intentionally unchanged).

    Degree >29: own CC. Direct neighbour degree >29: search 1 hop.
    Otherwise search up to 2 edges of distance <=2.
    Candidates must differ by <=2 loci directly from query and have degree >2.
    Choose maximum-degree candidate. Ties follow BFS discovery order.

    This is a *custom CC designation*, not ordinary connected components,
    and is sensitive to the reference collection / MST topology.
    """
    st_to_id = {st: nid for nid, st in id_to_st.items()}
    cc_assignments = {}
    missing = []
    for st in st_order:
        nid = st_to_id.get(st)
        if nid is None:
            missing.append(st)
            cc_assignments[st] = UNGROUPED
            continue
        if degree.get(nid, 0) > 29:
            cc_assignments[st] = "CC" + st
            continue
        direct_major_founder = any(degree.get(nbr, 0) > 29
                                   for nbr, edge_d in adj[nid] if edge_d <= 2)
        max_hops = 1 if direct_major_founder else 2
        reachable = {nid: 0}
        queue = deque([(nid, 0)])
        while queue:
            cur, hops = queue.popleft()
            if hops >= max_hops:
                continue
            for nbr, edge_d in adj[cur]:
                if nbr not in reachable and edge_d <= 2:
                    reachable[nbr] = hops + 1
                    queue.append((nbr, hops + 1))
        query_alleles = st_alleles.get(st)
        valid = [n for n in reachable
                 if allele_diff(query_alleles, st_alleles.get(id_to_st[n])) <= 2
                 and degree.get(n, 0) > 2]
        if not valid:
            cc_assignments[st] = UNGROUPED
            continue
        best = max(valid, key=lambda n: degree.get(n, 0))
        cc_assignments[st] = "CC" + id_to_st[best]
    return cc_assignments, missing


def run_from_alleles(st_alleles, st_order):
    profiles = [st_alleles[st] for st in st_order]
    idx_to_st = {i: st for i, st in enumerate(st_order)}
    tree = goeburst_mst(profiles, st_order)
    print(f"  MST complete: {len(tree)} edges for {len(profiles)} STs.")
    adj = defaultdict(list)
    for u, v, d in tree:
        adj[u].append((v, d))
        adj[v].append((u, d))
    degree = {i: len(adj[i]) for i in range(len(profiles))}
    st_to_deg = {idx_to_st[i]: degree[i] for i in range(len(profiles))}
    return adj, degree, idx_to_st, st_to_deg


def run_from_json(st_alleles, json_path):
    data = load_json(json_path)
    nodes = data["nodes"]
    edges = data["edges"]
    print(f"  {len(nodes)} nodes, {len(edges)} edges.")
    id_to_st = {}
    for node in nodes:
        nid = node["id"]
        st = str(node["profile"])
        if nid in id_to_st:
            raise ValueError(f"Duplicate PHYLOViZ node id: {nid}")
        id_to_st[nid] = st
    if len(set(id_to_st.values())) != len(id_to_st):
        raise ValueError("Duplicate ST labels in PHYLOViZ JSON")
    adj = defaultdict(list)
    for edge in edges:
        u, v = edge["u"], edge["v"]
        if u not in id_to_st or v not in id_to_st:
            raise ValueError(f"PHYLOViZ edge references unknown node: {u}, {v}")
        d = allele_diff(st_alleles.get(id_to_st[u]), st_alleles.get(id_to_st[v]))
        if d == 99:
            raise ValueError("PHYLOViZ contains an ST missing from reference: "
                             f"{id_to_st[u]} or {id_to_st[v]}")
        adj[u].append((v, d))
        adj[v].append((u, d))
    degree = {nid: len(adj[nid]) for nid in id_to_st}
    st_to_deg = {st: degree[nid] for nid, st in id_to_st.items()}
    return adj, degree, id_to_st, st_to_deg


def summary(cc_assignments):
    counts = Counter(cc_assignments.values())
    return {"assigned": sum(v for k, v in counts.items() if k != UNGROUPED),
            "ungrouped": counts.get(UNGROUPED, 0),
            "unique_ccs": len([k for k in counts if k != UNGROUPED]),
            "top_10": sorted([(k, v) for k, v in counts.items() if k != UNGROUPED],
                             key=lambda item: -item[1])[:10]}


def compare_with_reference(path, cc_assignments):
    """Check this run against frozen, previously published ST->CC assignments."""
    delim = detect_delimiter(path)
    previous = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        cols = reader.fieldnames or []
        st_col = _find_column(cols, {"st"})
        cc_col = _find_column(cols, {"cc"})
        if st_col is None or cc_col is None:
            raise ValueError("--compare-to CSV must contain ST and CC columns")
        for row in reader:
            st = (row.get(st_col) or "").strip()
            if st in previous:
                raise ValueError(f"Duplicate ST in --compare-to: {st}")
            previous[st] = (row.get(cc_col) or "").strip()
    missing_old = sorted(set(cc_assignments) - set(previous))
    missing_new = sorted(set(previous) - set(cc_assignments))
    differences = [(st, previous[st], cc_assignments[st]) for st in previous.keys() & cc_assignments.keys()
                   if previous[st] != cc_assignments[st]]
    differences.sort(key=lambda item: (0, int(item[0])) if item[0].isdigit() else (1, item[0]))
    print(f"Reference comparison: {len(differences)} different CCs, "
          f"{len(missing_old)} STs only in new, {len(missing_new)} only in old")
    for st, old_cc, new_cc in differences[:12]:
        print(f"    ST {st}: previous={old_cc}, current={new_cc}")
    return not (missing_old or missing_new or differences)


def derived_path(path, suffix):
    p = Path(path)
    return str(p.with_name((p.stem if p.suffix.lower() == ".csv" else p.name) + suffix))


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Assign NG-STAR CCs using goeBURST-style MST and the custom "
                    "Golparian founder-selection rules."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--alleles", required=True, help="NG-STAR reference ST and 7-locus CSV/TSV")
    parser.add_argument("--json", help="Optional precomputed PHYLOViZ full-MST JSON")
    parser.add_argument("--output", default="ngstar_cc_output.csv", help="Reference ST-to-CC CSV")
    parser.add_argument("--isolates", help="Optional isolate ID,ST CSV/TSV")
    parser.add_argument("--isolate-output", help="Optional isolate CSV output path")
    parser.add_argument("--first-column-is-st", action="store_true",
                        help="Use first reference column as ST if header is not recognised")
    parser.add_argument("--compare-to", help="Check new ST-to-CC map against an old CSV; "
                        "exit with error if any ST/CC differs")
    parser.add_argument("--metadata", help="Optional JSON provenance report with hashes/results")
    parser.add_argument("--database-date", help="Download date (YYYY-MM-DD) for metadata")
    args = parser.parse_args(argv)
    if args.database_date and not re.fullmatch(r"\d{4}-\d{2}-\d{2}", args.database_date):
        parser.error("--database-date must be YYYY-MM-DD")

    print(f"Loading alleles from : {args.alleles}")
    st_alleles, st_order = load_alleles(args.alleles, args.first_column_is_st)
    print(f"  {len(st_order)} STs loaded.")
    if args.json:
        print(f"Loading MST from     : {args.json}")
        adj, degree, id_to_st, st_to_deg = run_from_json(st_alleles, args.json)
    else:
        print("Building goeBURST-style full MST from alleles ...")
        print("  (large databases require many pairwise comparisons)")
        adj, degree, id_to_st, st_to_deg = run_from_alleles(st_alleles, st_order)
    print("Assigning CCs ...")
    cc_assignments, missing = assign_ccs_from_adj(st_alleles, st_order, adj,
                                                   degree, id_to_st)
    if missing:
        print(f"  Warning: {len(missing)} reference STs absent from supplied MST: "
              f"{missing[:10]}")
    totals = summary(cc_assignments)
    print(f"  {totals['assigned']} STs assigned to {totals['unique_ccs']} unique CCs.")
    print(f"  {totals['ungrouped']} STs Ungrouped.")
    print("  Top 10 CCs:")
    for cc, count in totals["top_10"]:
        print(f"    {cc}: {count} STs")

    # Compare before writing, so a failed reproducibility check cannot replace
    # the previously published output.
    if args.compare_to and not compare_with_reference(args.compare_to, cc_assignments):
        raise ValueError("CC results disagree with --compare-to. Original files unchanged.")

    write_output(args.output, st_order, st_alleles, cc_assignments)
    print(f"Output written to    : {args.output}")
    counts_path = derived_path(args.output, "_counts.csv")
    write_counts(counts_path, st_order, st_to_deg)
    print(f"Link counts written  : {counts_path}")

    isolate_stats = None
    isolate_output = None
    if args.isolates:
        isolate_output = args.isolate_output or derived_path(args.output, "_isolates.csv")
        print(f"Assigning isolate CCs: {args.isolates}")
        isolate_stats = assign_ccs_to_isolates(args.isolates, cc_assignments,
                                                isolate_output)
        for label, key in (("Total isolates", "total"), ("Assigned CC", "assigned"),
                           ("Ungrouped ST", "ungrouped"), ("ST not found", "not_found"),
                           ("Missing ST", "missing_st")):
            print(f"  {label:14}: {isolate_stats[key]}")
        print(f"Isolate output written: {isolate_output}")

    if args.metadata:
        record = {
            "tool": "ngstar_ccer", "version": VERSION,
            "algorithm": "custom NG-STAR CC founder assignment on goeBURST-style full MST",
            "reference": {"path": str(args.alleles), "sha256": sha256_file(args.alleles),
                          "database_date": args.database_date, "num_st": len(st_order)},
            "script_sha256": sha256_file(__file__),
            "python_version": platform.python_version(),
            "mode": "PHYLOViZ_JSON" if args.json else "internal_MST",
            "phyloviz_json_sha256": sha256_file(args.json) if args.json else None,
            "reference_summary": totals,
            "isolate_summary": isolate_stats,
            "isolate_input_sha256": sha256_file(args.isolates) if args.isolates else None,
            "outputs": {"reference": str(args.output), "counts": counts_path,
                        "isolates": isolate_output},
        }
        with open(args.metadata, "w", encoding="utf-8") as fh:
            json.dump(record, fh, indent=2, ensure_ascii=False)
            fh.write("\n")
        print(f"Provenance written to  : {args.metadata}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
