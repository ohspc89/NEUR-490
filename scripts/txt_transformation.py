"""
ELAN output preprocessing -> time-binned SCVs per tier.

Features:
    - Configurable bin interval
    - Center-in-interval vs. any-overlap bin assignment
    - Explicit unknown-label handling
    - Optional "last-zero fix" toggle for LA/RA
"""
import argparse
from typing import Optional, Tuple
from pathlib import Path
import pandas as pd
import numpy as np

INFANT_BEHAV_CODE = {
    "MDG": 1,
    "MDT": 1,
    "MDX": 1,   # directed movement still
    "MO": 2,
    "N": 0,
    "Z": 4,
    "Q": 4
}

# R-tier label groupings
R_POSITIVE = {"P"}
R_NEG_LA = {"NPB", "NPL"}
R_NEG_RA = {"NPB", "NPR"}

INFANT_VALID_MOVES = ["MDG", "MDT", "MDX"]
PRE_BUFFER_ms = 500

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="ELAN output preprocessing")
    p.add_argument("input_dir", help="directory containing .txt files")
    p.add_argument("output_dir", help="directory to save .csv files")
    p.add_argument("--interval-ms", type=int, default=100,
                   help="bin width in milliseconds (default: 100)")
    p.add_argument("--unknown-value", type=int, default=7,
                   help="value to assign for unknown/unmapped labels (default: 7)")
    p.add_argument("--overlap-mode", choices=["center", "any"], default="center",
                   help="bin assignment rule: 'center' (default) or 'any' overlap")
    # When the last window is shorter than 50 ms, it will always be labeled
    # with the unknown-value. That should not change our CRQA result.
    # p.add_argument("--skip-last-zero-fix", action="store_true",
    #                help="disable the LA/RA last-bin 0->previous correction")
    p.add_argument("--tail-short-threshold-ms", type=int, default=None,
                   help="If the leftover tail at the end is shorter than this threshold, "
                   "copy the preceding bin value into the last bin. "
                   "Default: half of --interval-ms.")
    return p.parse_args()


def load_subset(filepath: Path, tier: str) -> pd.DataFrame:
    """Read whitespace-delimited ELAN export and filter by tier."""
    names = ["Tier", "code", "start_ms", "end_ms", "duration_ms", "label"]
    codes = pd.read_csv(filepath, sep=r"\s+", header=None, names=names,
                        engine="python")

    subset = (
        codes.loc[codes["Tier"] == tier, ["start_ms", "end_ms", "label"]]
        .dropna()
        .reset_index(drop=True)
    )
    # This is unnecessary, but just in case
    subset["label"] = subset["label"].astype(str).str.strip().str.upper()
    return subset


def build_bins(end_time_ms: int, interval_ms: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return (left_edges, centers) arrays for bins spanning [0, end_time]."""
    t = np.arange(0, end_time_ms + interval_ms, interval_ms, dtype=int)
    centers = t[:-1] + interval_ms // 2
    return t, centers


def assign_bins_center(centers: np.ndarray, start: int, end: int) -> np.ndarray:
    """Mask bins whose centers fall in [start, end)."""
    return (centers >= start) & (centers < end)


def assign_bins_any(left_edges: np.ndarray, interval_ms: int, start: int, end: int) -> np.ndarray:
    """Mask bins whose intervals [L, L+interval) overlap [start, end)."""
    right_edges = left_edges[:-1] + interval_ms
    return ~((right_edges <= start) | (left_edges[:-1] >= end))


def code_infant_behav(label: str, default: Optional[int] = None) -> Optional[int]:
    if label is None:
        return default
    return INFANT_BEHAV_CODE.get(label, default)


def count_episodes(sub_r):
    """
    Count episodes, which are defined as:
        P onset -> first subsequent NP* onset
    NP* includes NPR, NPL, NPB.
    This operationally marks the end of a presentation episode,
    regardless of which hand caused the termination.
    """
    sub_r = sub_r.sort_values("start_ms")
    category = ["NP", "NPB", "NPR", "NPL"]
    episodes = []

    for i, (_, rr) in enumerate(sub_r.iterrows()):
        if rr["label"] in category:
            if i == 0:
                continue
            else:
                if sub_r["label"].iloc[i-1] == "P":
                    episodes.append((sub_r["start_ms"].iloc[i-1],
                                     sub_r["end_ms"].iloc[i-1],
                                     rr["start_ms"]))
    return episodes


def truncate_LA_RA(sub: pd.DataFrame, episodes: list) -> pd.DataFrame:
    """
    Infant MD* intervals are truncated at NP onset
    to approximate the end of the reaching phase.
    The remaining segment after NP onset (INV) reflects
    holding/manipulation and is excluded from reach analysis.
    """
    out = []
    drop_idx = []
    for k, ir in sub.iterrows():
        if ir["label"] not in INFANT_VALID_MOVES:
            continue

        best = None
        for ep_idx, ep in enumerate(episodes):
        # compute overlap with P window
            ov = max(0, min(ir["end_ms"], ep[1]) - max(ir["start_ms"], ep[0]))
            if ov <= 0:
                continue

            # optional: onset must not precede P too much
            if ir["start_ms"] < ep[0] - PRE_BUFFER_ms:
                continue

            if best is None or ov > best["ov"]:
                best = {"idx": ep_idx, "ov": ov, "np_on": ep[2]}

        if best is None:
            continue

        new_off = min(ir["end_ms"], best["np_on"])
        if new_off > ir["start_ms"]:
            new_row1 = dict(zip(ir.index,
                            [ir["start_ms"], new_off, ir['label']]
                           )
                       )

            # INV intervals are kept as a distinct level.
            # In subsequent CRQA preprocessing, only 1-1 matches are counted,
            # so INV (coded as non-7) does not contribute to recurrence.
            new_row2 = dict(zip(ir.index,
                            [new_off, ir["end_ms"], "INV"]
                           )
                       )
            out.append(new_row1)
            if new_off != ir["end_ms"]:
                out.append(new_row2)

            drop_idx.append(k)

    to_add = pd.DataFrame(out)
    sub_filtered = sub.drop(index=drop_idx)
    final = pd.concat((sub_filtered, to_add))

    return final.sort_values(by="start_ms")


def process_LA_RA(subset: pd.DataFrame,
                  interval_ms: int,
                  overlap_mode: str,
                  unknown_value: Optional[int],
                  tail_short_threshold_ms: Optional[int]) -> pd.DataFrame:
    if subset.empty or subset["end_ms"].isna().all():
        return pd.DataFrame()

    end_time = int(np.ceil(subset['end_ms'].max() / interval_ms) * interval_ms)
    left_edges, centers = build_bins(end_time, interval_ms)
    series = np.full(len(centers), unknown_value if unknown_value is not None else np.nan, dtype=float)

    for _, row in subset.iterrows():
        if overlap_mode == "center":
            mask = assign_bins_center(centers, int(row.start_ms), int(row.end_ms))
        else:
            mask = assign_bins_any(left_edges, interval_ms, int(row.start_ms),
                                   int(row.end_ms))
        val = code_infant_behav(row.label, unknown_value)
        series[mask] = val

    # Tail-short fix:
    # If leftover time after the last annotation  is shorter than threshold,
    # copy the preceding bin's value into the last bin.
    if len(series) >= 2:
        last_annot_end = int(subset["end_ms"].max())
        tail = end_time - last_annot_end
        threshold = tail_short_threshold_ms if tail_short_threshold_ms is not None else (interval_ms // 2)
        if 0 < tail < threshold:
            series[-1] = series[-2]

    out = pd.DataFrame({
        "time_ms": left_edges[:-1],
        "center_ms": centers,
        "binary_value": series.astype(int, copy=False)
        })
    return out


def process_R(subset: pd.DataFrame,
              interval_ms: int,
              overlap_mode: str,
              unknown_value: int,
              tail_short_threshold_ms: Optional[int]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if subset.empty or subset["end_ms"].isna().all():
        return pd.DataFrame(), pd.DataFrame()

    end_time = int(np.ceil(subset["end_ms"].max() / interval_ms) * interval_ms)
    left_edges, centers = build_bins(end_time, interval_ms)
    series_L = np.full(len(centers), unknown_value, dtype=int)
    series_R = np.full(len(centers), unknown_value, dtype=int)

    for _, row in subset.iterrows():
        if overlap_mode == "center":
            mask = assign_bins_center(centers, int(row.start_ms),
                                      int(row.end_ms))
        else:
            mask = assign_bins_any(left_edges, interval_ms, int(row.start_ms),
                                   int(row.end_ms))
        lbl = row.label

        series_L[mask] = 1 if lbl in R_POSITIVE else (0 if lbl in R_NEG_LA else unknown_value)
        series_R[mask] = 1 if lbl in R_POSITIVE else (0 if lbl in R_NEG_RA else unknown_value)

    # Tail-short fix for both R outputs
    last_annot_end = int(subset["end_ms"].max())
    tail = end_time - last_annot_end
    threshold = tail_short_threshold_ms if tail_short_threshold_ms is not None else (interval_ms // 2)
    if 0 < tail < threshold and len(centers) >= 2:
        series_L[-1] = series_L[-2]
        series_R[-1] = series_R[-2]

    df_L = pd.DataFrame({"time_ms": left_edges[:-1], "center_ms": centers,
                         "binary_value": series_L})
    df_R = pd.DataFrame({"time_ms": left_edges[:-1], "center_ms": centers,
                         "binary_value": series_R})
    return df_L, df_R


def main():
    args = parse_args()
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    if not input_dir.is_dir():
        raise IOError(f"No directory exists: {input_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    txt_files = sorted(input_dir.glob("*.txt"))
    if not txt_files:
        print(f"[WARN] No .txt files in {input_dir}")
        return

    for file in txt_files:
        print(f"\n[FILE] {file.name}")
        # We first need to trim LA and RA tier
        sub_R  = load_subset(file, "R")
        sub_LA = load_subset(file, "LA")
        sub_RA = load_subset(file, "RA")

        episodes = count_episodes(sub_R)

        # Truncate and update the original DataFrames
        sub_LA = truncate_LA_RA(sub_LA, episodes)
        sub_RA = truncate_LA_RA(sub_RA, episodes)

        # sub_LA.to_csv(file.stem + "_LA_new.csv")
        # sub_RA.to_csv(file.stem + "_RA_new.csv")

        # R tier
        df_R_L, df_R_R = process_R(
            sub_R, args.interval_ms, args.overlap_mode, args.unknown_value,
            args.tail_short_threshold_ms
        )
        if not df_R_L.empty:
            df_R_L.to_csv(output_dir / f"{file.stem}_R_for_LA.csv", index=False)
            df_R_R.to_csv(output_dir / f"{file.stem}_R_for_RA.csv", index=False)
            print(f"  [OK] R tiers -> {file.stem}_R_for_LA.csv, {file.stem}_R_for_RA.csv")
        else:
            print("  [SKIP] No valid R tier rows or missing end_ms.")

        # LA tier
        df_LA = process_LA_RA(
            sub_LA, args.interval_ms, args.overlap_mode,
            args.unknown_value,
            tail_short_threshold_ms=args.tail_short_threshold_ms
        )
        if not df_LA.empty:
            df_LA.to_csv(output_dir / f"{file.stem}_LA.csv", index=False)
            print(f"  [OK] LA -> {file.stem}_LA.csv")
        else:
            print("  [SKIP] No valid LA tier rows or missing end_ms.")

        # RA tier
        df_RA = process_LA_RA(
            sub_RA, args.interval_ms, args.overlap_mode,
            args.unknown_value,
            tail_short_threshold_ms=args.tail_short_threshold_ms
        )
        if not df_RA.empty:
            df_RA.to_csv(output_dir / f"{file.stem}_RA.csv", index=False)
            print(f"  [OK] RA -> {file.stem}_RA.csv")
        else:
            print("  [SKIP] No valid RA tier rows or missing end_ms.")


if __name__ == "__main__":
    main()
