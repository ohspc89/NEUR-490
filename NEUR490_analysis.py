"""
Code to batch-process .txt files
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

INFANT_BEHAV_CODE = {
    "MDG": 1,
    "MDT": 1,
    "MDX": 2,
    "MO": 3,
    "N": 4,
    "Z": 4,
    "Q": 4
}

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="ELAN output preprocessing")
    parser.add_argument("input_dir", help="directory containing .txt files")
    parser.add_argument("output_dir", help="directory to save .csv files")

    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()

    if not input_dir.is_dir():
        raise IOError(f"No directory exists: {input_dirr}")


    if not output_dir.is_dir():
        output_dir.mkdir(parents=True, exist_ok=True)


    def make_subset(file, tier="R"):
        """Filter rows of a given tier"""
        names = ["Tier", "code", "start_ms", "end_ms", "duration_ms", "label"]
        codes = pd.read_csv(file, sep=r"\s+", header=None, names=names)

        subset = (
            codes.loc[codes["Tier"] == tier, ["start_ms", "end_ms", "label"]]
            .dropna()
            .reset_index(drop=True)
        )

        return subset


    def code_infant_behav(label: str, default: int | None = None) -> int | None:
        if label is None:
            return default
        return INFANT_BEHAV_CODE.get(label, default)

    txt_files = input_dir.glob("*.txt")

    # loop through each file
    for file in txt_files:
        for tier in ["R", "LA", "RA"]:
            subset = make_subset(file, tier)

            if subset.empty or subset['end_ms'].isna().all():
                print(f"Skipping {file.name} — no valid tier {tier} or missing end_ms.")
                continue

            # sample grid
            interval = 100  # ms
            end_time = int(np.ceil(subset['end_ms'].max() / interval) * interval)
            t = np.arange(0, end_time + interval, interval, dtype=int)
            centers = t[:-1] + interval // 2

            # series
            series = np.zeros(len(centers), dtype=int)
            for _, row in subset.iterrows():
                mask = (centers >= row.start_ms) & (centers < row.end_ms)
                label = str(row.label).strip().upper()
                print(label)
                if tier == "R":
                    series[mask] = 1 if label == 'P' else 0
                else:
                    series[mask] = code_infant_behav(label)

            # edge case?
            if tier is not "R" and series[-1] == 0:
                series[-1] = series[-2]

        # df
        binary_df = pd.DataFrame({
            'time_ms': t[:-1],
            'center_ms': centers,
            'binary_value': series
        })

        # save to CSV
        base_name = file.stem + f"_{tier}.csv"
        binary_df.to_csv(output_dir / base_name, index=False)

        print(f"{output_dir / base_name} Saved to {output_dir}")
