import pandas as pd
import numpy as np
from pathlib import Path
import argparse

REACH_SUCCESS = {"MDG", "MDT"}
REACH_ATTEMPT = {"MDX"}
REACH_ANY = REACH_SUCCESS | REACH_ATTEMPT


def read_elan_txt(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["tier", "file_id", "start_ms", "end_ms", "dur_ms", "label"]
    )

    df["start_ms"] = df["start_ms"].astype(int)
    df["end_ms"] = df["end_ms"].astype(int)

    return df


def get_reach_onsets(df):
    return df[df["label"].isin(REACH_ANY)].sort_values("start_ms")


def compute_hand_metrics(hand_df, p_on, p_off):

    reaches = hand_df[hand_df["label"].isin(REACH_ANY)]

    reach_onsets = reaches[
        (reaches["start_ms"] >= p_on) &
        (reaches["start_ms"] < p_off)
    ]

    reach_onsets_before_p = reaches[
        (reaches["start_ms"] < p_on) &
        (reaches["end_ms"] >= p_on)
    ]

    success_onsets = reaches[
        (reaches["label"].isin(REACH_SUCCESS)) &
        (reaches["start_ms"] >= p_on) &
        (reaches["start_ms"] < p_off)
    ]

    success_onsets_before_p = reaches[
        (reaches["label"].isin(REACH_SUCCESS)) &
        (reaches["start_ms"] < p_on) &
        (reaches["end_ms"] >= p_on)
    ]

    responded = int(len(reach_onsets) > 0)
    responded_before_p = int(len(reach_onsets_before_p) > 0)

    latency = np.nan
    if responded:
        latency = reach_onsets.iloc[0]["start_ms"] - p_on
    # if there's MDX that started earlier than p_on,
    # and it lasted past p_on, set the latency = 0
    elif responded_before_p:
        latency = 0

    latency_success = np.nan
    if len(success_onsets) > 0:
        latency_success = success_onsets.iloc[0]["start_ms"] - p_on
    elif len(success_onsets_before_p) > 0:
        latency_success = 0

    return {
        "responded": responded if responded else responded_before_p,
        "latency": latency,
        "latency_success": latency_success,
        "reach_count": len(reach_onsets) if responded else len(reach_onsets_before_p),
        "success_count": len(success_onsets) if len(success_onsets) > 0 else len(success_onsets_before_p),
        "attempt_count": len(reach_onsets[reach_onsets["label"].isin(REACH_ATTEMPT)])
    }


def process_file(path):

    df = read_elan_txt(path)

    r = df[df["tier"] == "R"].sort_values("start_ms")
    la = df[df["tier"] == "LA"]
    ra = df[df["tier"] == "RA"]

    episodes = r[r["label"] == "P"]

    rows = []

    for i, row in enumerate(episodes.itertuples(), start=1):

        p_on = row.start_ms
        p_off = row.end_ms

        la_metrics = compute_hand_metrics(la, p_on, p_off)
        ra_metrics = compute_hand_metrics(ra, p_on, p_off)

        any_responded = int(
            la_metrics["responded"] or ra_metrics["responded"]
        )

        latency_candidates = [
            la_metrics["latency"],
            ra_metrics["latency"]
        ]

        latency_candidates = [x for x in latency_candidates if not pd.isna(x)]

        any_latency = np.nan
        if latency_candidates:
            any_latency = min(latency_candidates)

        rows.append({

            "source_file": path.name,
            "episode_idx": i,

            "p_on_ms": p_on,
            "p_off_ms": p_off,
            "p_dur_ms": p_off - p_on,

            "any_responded": any_responded,
            "any_latency_to_first_reach_ms": any_latency,

            "la_responded": la_metrics["responded"],
            "la_latency_to_first_reach_ms": la_metrics["latency"],
            "la_latency_to_first_success_ms": la_metrics["latency_success"],
            "la_reach_onset_count": la_metrics["reach_count"],
            "la_success_onset_count": la_metrics["success_count"],

            "ra_responded": ra_metrics["responded"],
            "ra_latency_to_first_reach_ms": ra_metrics["latency"],
            "ra_latency_to_first_success_ms": ra_metrics["latency_success"],
            "ra_reach_onset_count": ra_metrics["reach_count"],
            "ra_success_onset_count": ra_metrics["success_count"],
        })

    return pd.DataFrame(rows)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
        help="folder with ELAN txt files"
    )

    parser.add_argument(
        "--pattern",
        default="*_CC_RE.txt"
    )

    parser.add_argument(
        "--out",
        default="episode_metrics_by_hand.csv"
    )

    args = parser.parse_args()
    path = Path(args.input)
    files = list(path.glob(args.pattern))

    all_ep = []

    for f in files:
        ep = process_file(f)
        all_ep.append(ep)

    df = pd.concat(all_ep)
    # Sort!
    df = df.sort_values(by=['source_file', 'episode_idx'])
    df.to_csv(args.out, index=False)
    print("Wrote", args.out)


if __name__ == "__main__":
    main()
