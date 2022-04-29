import pandas as pd


def filter_outlier(df: pd.DataFrame, column: str) -> pd.DataFrame:
    q25, q75 = df[column].quantile([0.25, 0.75])
    qrange = (q75 - q25) * 1.5
    upper = q75 + qrange
    lower = q25 - qrange
    # return filter_df, upperlimit, lowerlimit
    return df.loc[df[column].between(lower, upper), :]
