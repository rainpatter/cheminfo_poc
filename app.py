from flask import Flask, render_template

import test_mod

app = Flask(__name__)

@app.route("/")
def home():
    return render_template('index.html', name=test_mod.test_variable)

@app.route("/ames")
def ames_data():
    return render_template('ames_index.html', tables=[test_mod.df_html], titles=test_mod.col_vals)

print(test_mod.df_html)

if __name__ == "__main__":
    app.run()

