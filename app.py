from flask import Flask, render_template

import data_functions


app = Flask(__name__)

@app.route("/")
def home():
    return render_template('index.html')

@app.route("/ames")
def ames_data():
    return render_template('ames_index.html', json=[data_functions.df_json], tables=[data_functions.df_html], titles=data_functions.col_vals)

@app.route("/ames/<chem_id>")
def show_ames_data(chem_id):
    print(type(chem_id))
    chem_dict = data_functions.match_on_id(data_functions.df_json, chem_id)
    return render_template('ames_show.html', id=chem_id, dict=chem_dict)

if __name__ == "__main__":
    app.run()

