import pandas as pd
import json

df = pd.read_csv('data/test_data_structure.csv', index_col=False)
df_html = df.to_html(classes="data")
df_json = json.loads(df.reset_index().to_json(orient="records"))
col_vals = df.columns.values

def match_on_id(dict_list, id):
    for value in dict_list:
        if str(value['index']) == id:
            match_dict = value
        else:
            continue
    return match_dict
