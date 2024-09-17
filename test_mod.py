import pandas as pd

df = pd.read_csv('data/test_data.csv')
df_html = df.to_html(classes="data")
col_vals = df.columns.values

test_variable = "new test"