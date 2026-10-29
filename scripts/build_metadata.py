#!/usr/bin/env python3
"""Draft metadata/samples.tsv from the ENA run table of a BioProject / study.

Example:
    python scripts/build_metadata.py PRJNAxxxxxx --ko-pattern mig1

The script labels each run WT / KO by regex-matching the sample text. Heuristics
can be wrong, so ALWAYS open the output and verify it against the GEO sample page.
"""
import argparse
import io
import os
import re
import sys

import pandas as pd
import requests

ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = [
    "run_accession", "sample_accession", "sample_alias", "sample_title",
    "experiment_title", "library_layout", "library_strategy",
    "instrument_model", "read_count", "fastq_ftp", "fastq_md5",
]


def fetch_runs(accession: str) -> pd.DataFrame:
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(FIELDS),
        "format": "tsv",
    }
    r = requests.get(ENA_URL, params=params, timeout=60)
    if r.status_code != 200:
        sys.exit(
            f"ENA returned HTTP {r.status_code}: {r.text[:300]}\n"
            "If a field name was rejected, list valid ones with:\n"
            "  curl 'https://www.ebi.ac.uk/ena/portal/api/returnFields?result=read_run'"
        )
    df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str).fillna("")
    if df.empty:
        sys.exit(f"No runs found for {accession}. Check the accession.")
    return df


def classify(text: str, wt_re: re.Pattern, ko_re: re.Pattern) -> str:
    is_wt, is_ko = bool(wt_re.search(text)), bool(ko_re.search(text))
    if is_ko and not is_wt:
        return "KO"
    if is_wt and not is_ko:
        return "WT"
    return "CHECK"  # ambiguous or no match -> you decide


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accession", help="BioProject/study accession, e.g. PRJNA..., ERP..., SRP...")
    ap.add_argument("--ko-pattern", required=True, help="regex for knockout samples, e.g. 'mig1'")
    ap.add_argument("--wt-pattern", default=r"\bwt\b|wild[ _-]?type", help="regex for WT samples")
    ap.add_argument("--out", default="metadata/samples.tsv")
    ap.add_argument("--raw-out", default="metadata/ena_filereport_raw.tsv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df = fetch_runs(args.accession)
    df.to_csv(args.raw_out, sep="\t", index=False)
    print(f"[info] {len(df)} runs fetched -> raw table saved to {args.raw_out}")

    df = df[df["library_strategy"].str.upper() == "RNA-SEQ"].copy()
    text = (df["sample_title"] + " " + df["sample_alias"] + " " + df["experiment_title"]).str.lower()
    wt_re = re.compile(args.wt_pattern, re.I)
    ko_re = re.compile(args.ko_pattern, re.I)
    df["condition"] = [classify(t, wt_re, ko_re) for t in text]
    df["replicate"] = df.groupby("condition").cumcount() + 1

    keep = ["run_accession", "sample_alias", "sample_title", "condition", "replicate",
            "library_layout", "read_count", "fastq_ftp", "fastq_md5"]
    df[keep].to_csv(args.out, sep="\t", index=False)

    print(f"[info] draft sample sheet written to {args.out}\n")
    print(df[["run_accession", "sample_title", "condition", "library_layout"]].to_string(index=False))
    print("\nlayouts:", df["library_layout"].value_counts().to_dict())
    n_check = (df["condition"] == "CHECK").sum()
    if n_check:
        print(f"\n[ACTION] {n_check} runs are labelled CHECK. Edit {args.out}: set condition to WT/KO, "
              "or delete rows you do not want (other strains, other temperatures).")
    print("[ACTION] Set paired_end in config.yaml to match library_layout.")


if __name__ == "__main__":
    main()
