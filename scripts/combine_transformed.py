"""
Combine crqa results into a single file
"""
from pathlib import Path
import pandas as pd

copilot_folder = Path('../data/processed/crqa_results')
copilot_files = copilot_folder.glob('*_crqa_results.csv')

copilot_rows = []
for cf in copilot_files:
    crow = pd.read_csv(cf)
    crow["id"] = cf.stem
    copilot_rows.append(crow)

copilot_df = pd.concat(copilot_rows)
# Shuffling columns
# copilot_df = copilot_df[['id', 'RR', 'DET', 'L', 'maxL',
#                          'ENTR', 'LAM', 'TT', 'max_vertlength',
#                          'peak_lag_bins', 'peak_lag_ms',
#                          'peak_rr', 'pos_auc', 'neg_auc',
#                          'asymmetry']]

df = copilot_df.sort_values('id')
df = df.reset_index(drop=True)

df[['pid', 'temp']] = df['id'].str.split('-', expand=True)
df[['month', 'activity', 'drop1', 'drop2', 'hand', 'drop3', 'drop4']] = df['temp'].str.split(r'_+', expand=True)

df = df.drop(columns=['temp', 'drop1', 'drop2', 'drop3', 'drop4'])

# Reshuffle columns
colnames = df.columns
to_pull_forward = ['pid', 'month', 'activity', 'hand']
reordered = to_pull_forward + [x for x in colnames if x not in to_pull_forward]
df = df[reordered]
df = df.drop(columns='id')

# R's crqa() will return RR that's between 0-100. Make it a rate, not a percentage.
df = df.copy()
df['RR'] = df['RR']/100

df.to_csv(copilot_folder / 'combined/combined_maxdiag5.csv', index=False)
