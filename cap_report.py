#!/usr/bin/env python3

import pandas as pd
import sys

EXCEL_MAX = 32000
REPORT = sys.argv[1]

try:
	df = pd.read_csv(REPORT, sep='\t')
except Exception as e:
	df = pd.read_csv(REPORT, encoding="ISO-8859-1", sep='\t')

for col in df.columns:
	df[col] = df[col].apply(lambda x: str(x)[:EXCEL_MAX-1])

df.to_csv(REPORT, index=False, sep='\t')
