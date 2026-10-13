# NG-STAR CCer

**Reproducible NG-STAR clonal-complex (CC) designation from reference allelic profiles, with optional isolate annotation.**

`ngstar_ccer.py` constructs a **goeBURST-style full minimum spanning tree (MST)** from the seven NG-STAR loci and applies the *custom founder-selection rules* used in the Euro-GASP workflow to produce an ST-to-CC lookup table. Alternatively, it can read an existing PHYLOViZ full-MST JSON export. A study-specific `ID,ST` table can then be annotated using that same reference mapping.

> **Scientific scope:** The script does **not** simply report the connected-component/CC groups output by PHYLOViZ and is not a generic or formally validated replacement for PHYLOViZ. It implements a specific CC **nomenclature heuristic** following full-MST construction. CC labels depend on the exact NG-STAR reference database snapshot and, in some ties, the MST topology/input order. Freeze the reference snapshot and check against previous assignments when comparing studies.

## Requirements

- Python 3.8 or later (standard library only; no `pip install` required)
- A complete NG-STAR reference-profile export with **one unique ST and one complete seven-locus allelic profile per row**
- Optional isolate input with an `ID` column and an `ST` column

The seven NG-STAR loci are `penA`, `mtrR`, `porB`, `ponA`, `gyrA`, `parC`, and `23S`. The script assigns CCs to *existing NG-STAR STs*; it **does not** perform NG-STAR sequence typing from genome assemblies. Use a separate typing program (e.g. `pyngSTar`) to obtain each isolate's ST first.

## Quick start

Download the full reference-profile table from the [NG-STAR database](https://ngstar.canada.ca/) and save the original export. The input may be comma-, semicolon-, or tab-separated, with a header.

Example reference (`profiles.csv`):

```csv
ST,penA,mtrR,porB,ponA,gyrA,parC,23S
1,34.001,5,1,1,1,1,1
2,34.001,5,1,1,1,1,2
3,34.001,5,1,1,1,2,1
```

Example isolate file (`isolates.csv`):

```csv
ID,ST
sample_A,1
sample_B,2
sample_C,999999
```

Run a complete analysis:

```bash
python ngstar_ccer.py \
  --alleles profiles.csv \
  --output profiles_to_cc.csv \
  --isolates isolates.csv \
  --isolate-output isolates_with_cc.csv \
  --database-date 2026-09-15 \
  --metadata run_metadata.json
```

This creates:

| File | Contents |
| --- | --- |
| `profiles_to_cc.csv` | Database ST, seven alleles, and assigned `CC` |
| `profiles_to_cc_counts.csv` | ST and number of MST edges (`LinkCount`) |
| `isolates_with_cc.csv` | Three columns: `ID,ST,CC` |
| `run_metadata.json` | Script version, SHA-256 checksums, reference download date, assignment summaries |

The isolate assignment statuses are `CC<founder ST>` for an assigned CC, `Ungrouped` for a reference ST without an assigned CC, `ST_not_found` for an ST absent from the reference snapshot, and `Missing_ST` for an empty ST.

### Direct import of NG-STAR exports

You **normally do not need** to manually convert a semicolon-separated export to commas or replace decimal commas with decimal points. The script detects `,`, `;`, or tab delimiters and normalises *allele values* such as `34,001` to `34.001`. In a comma-delimited CSV, embedded decimal commas must be properly quoted according to CSV standards.

If the first reference column is not called `ST`, `Sequence type`, or another recognised NG-STAR ST heading, explicitly identify the **first column** as ST:

```bash
python ngstar_ccer.py --alleles raw_ngstar_export.csv \
  --first-column-is-st --output profiles_to_cc.csv
```

**Do not blindly replace every comma or semicolon in the raw file.** Always verify that the number of imported STs agrees with the intended database snapshot.

### Reproduce the Euro-GASP analysis and compare against its original results

Use the *same* September 2026 profile export, *not* a newly downloaded database. Confirm that the program imports 8,282 STs (as in the recorded run) before drawing any conclusions about comparability.

```bash
python ngstar_ccer.py \
  --alleles ngstar_profiles_2026-09-15T07_51_00-05_00.csv \
  --output rerun_profiles_to_cc.csv \
  --isolates for_CC_designation.csv \
  --isolate-output rerun_CC_assigned.csv \
  --compare-to 260915_profiles_to_cc.csv \
  --database-date 2026-09-15 \
  --metadata eurogasp_cc_reproducibility.json
```

`--compare-to` compares **every ST and CC value** against a frozen previously generated `ST,CC` CSV. If any ST or CC differs, the program exits with an error **before writing new output files**. This is especially useful for verifying that code maintenance has not changed the classifications used in a manuscript.

### Use a precomputed PHYLOViZ full-MST JSON

```bash
python ngstar_ccer.py \
  --alleles profiles.csv \
  --json goeburst_full_mst.json \
  --output profiles_to_cc.csv
```

The JSON is expected to contain `nodes` (`id`, `profile`) and `edges` (`u`, `v`). For scientific comparability, the JSON and allele table should refer to the **same reference snapshot**. The script tolerates a limited legacy trailing-comma formatting defect in some PHYLOViZ JSON exports.

## Exact CC-designation logic

The graph is built using the original pipeline's full-MST construction: ascending 1–7-locus Hamming distance, followed by goeBURST-style edge priorities (numbers of SLVs/DLVs/etc., profile frequency, numeric ST identifier). After the tree has been built, its **MST degree** is counted for each ST.

For each reference ST:

1. If its MST degree is **>29**, designate it as its own major founder (`CC<ST>`).
2. Otherwise, if an immediate MST neighbour connected by an edge of **at most two allele differences** has degree **>29**, search at most **one MST edge**.
3. Otherwise, search up to **two MST edges**, traversing only edges with **at most two allele differences** each.
4. Among nodes discovered, retain only those with **more than two MST links** **and** at most **two allele differences directly relative to the queried ST** (not merely a short tree path).
5. Choose the candidate with the **highest MST degree** and assign `CC<candidate ST>`. If none qualifies, assign `Ungrouped`.
6. If several candidates tie on MST degree, retain the first discovered in breadth-first traversal, as in the original research script. This follows adjacency/MST edge ordering and is **not** a separate numerical-ST tiebreak.

The **>29 links** and **up-to-two-hop** rules are *custom implementation choices*. They should be explicitly stated when describing this pipeline; they are not a universal definition of goeBURST CCs. The same reference snapshot should be used for comparative historical datasets.

## Research/reproducibility notes

- Both database growth and differences between a PHYLOViZ-generated MST and the internal goeBURST-style MST can change CC assignments.
- A complete reference ST must have all seven allele values. Duplicated STs, duplicated seven-locus profiles, and incomplete allele profiles cause an error rather than silently influencing topology.
- This implementation is deliberately **CPU-oriented and O(N²)** in reference ST comparisons, with distance-specific candidate-edge sorting. For ~8,000 STs it can be resource-intensive. It is currently **single threaded**: providing many cores will not automatically parallelise the MST calculation.
- The MST `LinkCount` is the number of selected tree edges incident on an ST, **not** the number of observed isolates or the total count of possible single-locus neighbours.
- The manuscript must use the true counts reported by the frozen reference snapshot. The example Euro-GASP run supplied to this project recorded **8,282 STs, 743 named CCs, 7,982 STs assigned, and 300 ungrouped STs**. This should not be conflated with the number of isolates or with alternative database snapshots.
- The September log also reports 8,859 isolate rows processed, with 8,795 CC assignments, 36 ungrouped STs, 28 STs absent from the reference and zero missing STs. These represent the *combined input* and not the 2023 study alone.
- Avoid publishing any potentially sensitive isolate identifiers or unpublished study data in a public GitHub repository; the files in `examples/` are synthetic.

## Testing

```bash
python -m unittest -v test_ngstar_ccer.py
```

The tests include complete execution, the four isolate statuses, decimal comma input, header handling, detection of invalid references, PHYLOViZ JSON compatibility, the founder rules, and fail-fast comparison with previous CC assignments. Development regression testing against the previous script on three synthetic profiles sets (30, 75 and 120 STs) produced identical CCs. **The actual September 2026 database is not included in this repository and has not been re-tested here.**

## References

- Golparian D, Sánchez-Busó L, Cole M, Unemo M. *Neisseria gonorrhoeae* Sequence Typing for Antimicrobial Resistance (NG-STAR) clonal complexes are consistent with genomic phylogeny and provide simple nomenclature, rapid visualization and antimicrobial resistance (AMR) lineage predictions. *J Antimicrob Chemother*. 2021;76:940–944. https://doi.org/10.1093/jac/dkaa552
- Francisco AP, Bugalho M, Ramirez M, Carriço JA. Global optimal eBURST analysis of multilocus typing data using a graphic matroid approach. *BMC Bioinformatics*. 2009;10:152. https://doi.org/10.1186/1471-2105-10-152
- PHYLOViZ user guide: https://www.phyloviz.net/wiki/dataanalysis/

## License

Released under the [MIT License](LICENSE). Copyright (c) 2026 Daniel Golparian. See `LICENSE` for the full terms.
