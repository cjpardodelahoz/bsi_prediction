#!/usr/bin/env python3
"""
growth_plot.py - Command-line tool to plot OD600 growth curves from plate reader .xlsx files.
"""
import argparse, os, sys
from typing import List, Dict
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def find_header_row(path: str):
    tmp = pd.read_excel(path, header=None, nrows=200)
    rows = tmp[tmp.iloc[:,0] == "Cycle Nr."].index.tolist()
    if not rows:
        raise ValueError('Could not find a row with "Cycle Nr." in column A.')
    return rows[0]

def read_plate_data(path: str):
    header_row = find_header_row(path)
    df = pd.read_excel(path, header=header_row)
    first_col = df.columns[0]
    end_idx = df[df[first_col].astype(str).str.contains("End Time", na=False)].index
    if len(end_idx)>0:
        df = df.loc[:end_idx[0]-1].copy()
    return df

def infer_time_column(df: pd.DataFrame):
    for c in df.columns:
        if "time" in str(c).lower():
            return c
    if len(df.columns) >= 2:
        return df.columns[1]
    raise ValueError("Could not infer time column.")

def parse_grouping_strings(grouping_strs: List[str], df_columns: List[str]) -> Dict[str, List[str]]:
    groups = {}
    for s in grouping_strs or []:
        if "=" not in s:
            raise argparse.ArgumentTypeError(f"Grouping must be NAME=Well1,Well2,... Got: {s}")
        name, wells = s.split("=",1)
        wells_list = parse_well_list(wells, df_columns)
        if wells_list:
            groups[name.strip()] = wells_list
    return groups

def parse_well_list(s: str, df_columns: List[str]) -> List[str]:
    """
    Accepts well IDs (e.g., A1,B2) and/or column numbers (e.g., col2).
    Returns list of well IDs present in df_columns.
    """
    wells = []
    items = [item.strip() for item in s.split(',') if item.strip()]
    for item in items:
        if item.lower().startswith('col'):
            try:
                col_num = int(item[3:])
                well_ids = [f"{row}{col_num}" for row in "ABCDEFGH"]
                wells.extend([w for w in well_ids if w in df_columns])
            except Exception as e:
                print(f"Could not parse column number from '{item}': {e}", file=sys.stderr)
        else:
            if item in df_columns:
                wells.append(item)
    return wells

def compute_mean_se(df: pd.DataFrame, cols: List[str]):
    subset = df[cols].astype(float)
    mean = subset.mean(axis=1)
    sd = subset.std(axis=1, ddof=1)
    n = subset.count(axis=1)
    se = sd / np.sqrt(n.replace(0, np.nan))
    return mean, se, n

def main():
    p = argparse.ArgumentParser(description="Plot OD600 growth curves from plate-reader .xlsx exports.")
    p.add_argument("xlsx", help="Path to Excel file.")
    p.add_argument("--blank-wells", default="", help="Comma-separated well IDs used as blank.")
    p.add_argument("--grouping", action="append", default=[], help='Grouping: "NAME=Well1,Well2,...". Can repeat.')
    p.add_argument("--separate-plots", action="store_true", help="Plot each grouping separately.")
    p.add_argument("--outdir", default=".", help="Output directory.")
    p.add_argument("--time-as-minutes", action="store_true", help="Label x-axis as minutes if in seconds.")
    p.add_argument("--suffix", default="", help="Suffix to append to output filenames.")
    args = p.parse_args()

    df = read_plate_data(args.xlsx)
    time_col = infer_time_column(df)
    potential_wells = [c for c in df.columns if isinstance(c, str) and len(c)>=2 and c[0].isalpha() and any(ch.isdigit() for ch in c)]
    if not potential_wells:
        potential_wells = list(df.columns[2:])
    df_wells = df[potential_wells].apply(pd.to_numeric, errors='coerce')
    time = pd.to_numeric(df[time_col], errors='coerce')
    # Convert time from seconds to hours
    time = time / 3600
    blank_wells = parse_well_list(args.blank_wells, df_wells.columns)
    if blank_wells:
        blanks_present = [w for w in blank_wells if w in df_wells.columns]
        if blanks_present:
            blank_avg = df_wells[blanks_present].mean(axis=1)
            df_wells = df_wells.subtract(blank_avg, axis=0)
    groups = parse_grouping_strings(args.grouping, df_wells.columns)
    if not groups:
        for w in df_wells.columns:
            groups[w] = [w]
    os.makedirs(args.outdir, exist_ok=True)
    if args.separate_plots:
        for name, wells in groups.items():
            mean, se, n = compute_mean_se(df_wells, wells)
            fig, ax = plt.subplots(figsize=(8,5))
            ax.plot(time, mean, label=name)
            ax.fill_between(time, mean-se, mean+se, alpha=0.25)
            ax.set_xlabel("Time (hours)")
            ax.set_ylabel("OD600 (blank subtracted)")
            ax.set_title(name)
            ax.legend(); ax.grid(True)
            ax.set_xticks(np.arange(0, max(time)+2, 2))
            fname = os.path.join(args.outdir, f"{name.replace(' ','_')}{('_'+args.suffix) if args.suffix else ''}.png")
            fig.savefig(fname, dpi=300); plt.close(fig)
    else:
        fig, ax = plt.subplots(figsize=(10,6))
        for name, wells in groups.items():
            mean, se, n = compute_mean_se(df_wells, wells)
            ax.plot(time, mean, label=name)
            ax.fill_between(time, mean-se, mean+se, alpha=0.25)
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("OD600 (blank subtracted)")
        ax.set_title(os.path.basename(args.xlsx))
        ax.legend(); ax.grid(True)
        ax.set_xticks(np.arange(0, max(time)+2, 2))
        fname = os.path.join(args.outdir, f"{os.path.splitext(os.path.basename(args.xlsx))[0]}{('_'+args.suffix) if args.suffix else ''}.png")
        fig.savefig(fname, dpi=300); plt.close(fig)

if __name__ == "__main__":
    main()
