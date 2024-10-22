from flask import Flask, render_template, request, url_for, flash, redirect

import numpy as np

from ames.ames import AmesPredictor

from utils.data_handlers import generate_image
from utils.data_handlers import csv_to_dict
from targeted_risk_assessment.calculator import calculate_all

app = Flask(__name__)
app.config['SECRET_KEY'] = 'aicis'

np.set_printoptions(legacy='1.25')

all_results = []
all_images = []
hh_dict_result = []
calc_dict = {}

@app.route("/")
def home():
    return render_template('index.html')


@app.route('/ames_entry', methods=('GET', 'POST'))
def ames_entry():
    if request.method == 'POST':
        smiles = request.form['smiles']
        if not smiles:
            flash('SMILES are required')
        else:
            try:
                result = AmesPredictor(smiles)
                all_results.append(result)
                image = generate_image(smiles)
                all_images.append(image)

                return redirect(url_for('ames_data'))
            except:
                flash('Invalid SMILES string')
    return render_template('ames_form.html', list=all_results)


@app.route("/ames")
def ames_data():
    current_result = all_results[-1]
    current_image = all_images[-1]
    return render_template('ames_show.html', json=current_result, img_data=current_image.decode('utf-8'))


@app.route('/all_ames')
def all_ames_data():
    return render_template('ames_index.html', list=all_results)


@app.route('/about')
def about():
    return render_template('about.html')


@app.route('/model')
def model():
    return render_template('model.html')


@app.route('/high_hazard_chemicals')
def hh_chemicals():
    chem_dict = csv_to_dict('./data/high_hazard_chemicals.csv')
    return render_template('high_hazard.html', results=chem_dict)


@app.route('/high_hazard_chemical_result')
def hh_result():
    return render_template('hh_result.html', results=hh_dict_result)


@app.route('/high_hazard_search', methods=('GET', 'POST'))
def hh_search():
    if request.method == "POST":
        chem_dict = csv_to_dict('./data/high_hazard_chemicals.csv')
        cas = request.form['cas']
        chem_name = request.form['chem_name']
        if cas or chem_name:
            retrieved_dict = [dict for dict in chem_dict if dict['CAS RN'] == cas.strip(
            ) or dict['Chemical name'].upper() == chem_name.upper().strip()]
            if retrieved_dict:
                global hh_dict_result
                hh_dict_result = retrieved_dict
                return redirect(url_for('hh_result', results=hh_dict_result))
            else:
                flash('Search term not found')
        else:
            flash('A CAS number or chemical name is required')
    return render_template('hh_form.html')


@app.route('/exposure_calculator', methods=('GET', 'POST'))
def exposure_calc():
    if request.method == "POST":
        dict = request.form.to_dict()
        if any(x == '' for x in list(dict.values())):
            flash('Inputs are incomplete - please complete all entries')
        else:
            try:
                dict = {
                    'substance_name': request.form['substance_name'],
                    'cas_number': request.form['cas_no'],
                    'mol_weight': float(request.form['molecular_weight']),
                    'long_term_inhalation': float(request.form['lt_inhalation']),
                    'long_term_dermal': float(request.form['lt_dermal']),
                    'short_term_inhalation': float(request.form['st_inhalation']),
                    'local_dermal': float(request.form['local_dermal']),
                    'proc': request.form['proc'],
                    'ind_prof': request.form['ind_prof'],
                    'phys_state': request.form['phys_state'],
                    'fugacity': request.form['fugacity'],
                    'ventilation': request.form['ventilation'],
                    'duration': request.form['duration'],
                    'concentration': request.form['concentration'],
                    'lev': request.form['lev'],
                    'rpe_mask': request.form['rpe_mask'],
                    'ppe_gloves': request.form['ppe_gloves'],
                    'lev_dermal': request.form['lev_dermal'],
                }
                global calc_dict
                calc_dict = calculate_all(dict)
                print(calc_dict)
                return redirect(url_for('exposure_result', dict=calc_dict))
            except:
                flash('Inputs require changing for calculation')
    return render_template('exposure_form.html')

@app.route('/exposure_result')
def exposure_result():
    return render_template('exposure_result.html', dict=calc_dict)

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=5000, debug=True)

