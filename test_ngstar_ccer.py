"""Run: python -m unittest -v test_ngstar_ccer.py (standard-library only)."""
import csv
import json
import random
import tempfile
import unittest
from pathlib import Path

import ngstar_ccer as mod


class NGStarCCerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.alleles = self.root / "reference.csv"
        self.isolates = self.root / "isolates.csv"
        self.output = self.root / "results.csv"

    def write_profiles(self, rows, delimiter=",", header="ST", decimal_comma=False):
        with self.alleles.open("w", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh, delimiter=delimiter)
            wr.writerow([header] + mod.LOCI)
            for st, profile in rows:
                wr.writerow([st] + list(profile))

    def test_end_to_end_and_isolate_statuses(self):
        rows = [
            ("10", ("34.001", "1", "1", "1", "1", "1", "1")),
            ("11", ("34.001", "2", "1", "1", "1", "1", "1")),
            ("12", ("34.001", "3", "1", "1", "1", "1", "1")),
            ("13", ("34.001", "4", "1", "1", "1", "1", "1")),
            ("99", ("77.001", "9", "9", "9", "9", "9", "9")),
        ]
        self.write_profiles(rows)
        with self.isolates.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter=";")
            w.writerow(["ID", "ST"])
            w.writerows([["iso-1", "10"], ["iso-2", "99"], ["iso-3", "123456"],
                         ["iso-4", ""]])
        report = self.root / "audit.json"
        self.assertEqual(mod.main(["--alleles", str(self.alleles), "--output",
                                   str(self.output), "--isolates", str(self.isolates),
                                   "--metadata", str(report), "--database-date", "2026-09-15"]), 0)
        with self.output.open(newline="") as fh:
            mapping = {r["ST"]: r["CC"] for r in csv.DictReader(fh)}
        self.assertEqual(len(mapping), 5)
        self.assertEqual(mapping["99"], mod.UNGROUPED)
        with (self.root / "results_isolates.csv").open(newline="") as fh:
            statuses = {r["ID"]: r["CC"] for r in csv.DictReader(fh)}
        self.assertEqual(statuses["iso-1"], mapping["10"])
        self.assertEqual(statuses["iso-2"], mod.UNGROUPED)
        self.assertEqual(statuses["iso-3"], "ST_not_found")
        self.assertEqual(statuses["iso-4"], "Missing_ST")
        self.assertEqual(json.loads(report.read_text())["reference"]["num_st"], 5)
        # Frozen mapping matches on rerun.
        self.assertEqual(mod.main(["--alleles", str(self.alleles),
                                   "--compare-to", str(self.output),
                                   "--output", str(self.root / "rerun.csv")]), 0)

    def test_decimal_comma_and_header_alias(self):
        rows = [("1", ("34,001", "1", "1", "1", "1", "1", "1")),
                ("2", ("34,001", "2", "1", "1", "1", "1", "1"))]
        self.write_profiles(rows, delimiter=";", header="Sequence type")
        mapped, order = mod.load_alleles(self.alleles)
        self.assertEqual(order, ["1", "2"])
        self.assertEqual(mapped["1"][0], "34.001")

    def test_unknown_first_column(self):
        self.write_profiles([("1", ("10.001", "1", "1", "1", "1", "1", "1"))],
                            header="Profile")
        with self.assertRaisesRegex(ValueError, "ST column"):
            mod.load_alleles(self.alleles)
        mapped, _ = mod.load_alleles(self.alleles, first_column_is_st=True)
        self.assertIn("1", mapped)

    def test_duplicate_and_incomplete_reference_rejected(self):
        base = ("1", "1", "1", "1", "1", "1", "1")
        self.write_profiles([("1", base), ("1", ("2",) + base[1:])])
        with self.assertRaisesRegex(ValueError, "Duplicate NG-STAR ST"):
            mod.load_alleles(self.alleles)
        self.write_profiles([("1", base), ("2", base)])
        with self.assertRaisesRegex(ValueError, "identical"):
            mod.load_alleles(self.alleles)
        self.write_profiles([("1", base), ("2", ("",) + base[1:])])
        with self.assertRaisesRegex(ValueError, "Missing penA"):
            mod.load_alleles(self.alleles)

    def test_custom_cc_rules(self):
        a = {"1": ("1", "1", "1", "1", "1", "1", "1"),
             "2": ("1", "2", "1", "1", "1", "1", "1"),
             "3": ("1", "3", "1", "1", "1", "1", "1")}
        # Manually controlled MST: major founder 1 has 30 links.
        graph = {0: [(1, 1)], 1: [(0, 1), (2, 1)], 2: [(1, 1)]}
        degree = {0: 30, 1: 2, 2: 1}
        ccs, _ = mod.assign_ccs_from_adj(a, list(a), graph, degree,
                                          {0: "1", 1: "2", 2: "3"})
        self.assertEqual(ccs["1"], "CC1")
        self.assertEqual(ccs["2"], "CC1")
        self.assertEqual(ccs["3"], "CC1")  # 2-hop path reaches the major founder

    def test_phyloviz_json_mode(self):
        rows = [("1", ("1", "1", "1", "1", "1", "1", "1")),
                ("2", ("1", "2", "1", "1", "1", "1", "1"))]
        self.write_profiles(rows)
        tree = self.root / "tree.json"
        tree.write_text('{"nodes":[{"id":0,"profile":"1"},'
                        '{"id":1,"profile":"2"}],'
                        '"edges":[{"u":0,"v":1}],}')
        self.assertEqual(mod.main(["--alleles", str(self.alleles),
                                   "--json", str(tree), "--output", str(self.output)]), 0)
        self.assertTrue(self.output.exists())

    def test_compare_to_detects_changes(self):
        self.write_profiles([("1", ("1", "1", "1", "1", "1", "1", "1"))])
        old = self.root / "old.csv"
        old.write_text("ST,CC\n1,CC999\n")
        with self.assertRaisesRegex(ValueError, "disagree"):
            mod.main(["--alleles", str(self.alleles), "--output", str(self.output),
                      "--compare-to", str(old)])
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
