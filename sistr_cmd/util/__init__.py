# -*- coding: utf-8 -*-

import pandas as pd


def first_row(df: pd.DataFrame) -> pd.Series:
    for _, row in df.iterrows():
        return row