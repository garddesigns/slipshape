from flask import Flask, render_template, request
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import main_app as shiphulldes

app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True


@app.route('/')
def index():
    return render_template('index.html')

@app.route("/result", methods = ['POST', 'GET'])
def result():
    output = request.form.to_dict()
    dims = request.form.getlist('dimensions')
    print(type(dims))
    speed = float(output["speed"])
    DWT = int(output["DWT"])
    ship_type = output['type']
    fuel = output['fuel']

    n = 7

    shiphulldes.shipdesign(speed, DWT, ship_type, fuel, {}, dims, [], n)

    return render_template('index.html', dims = dims, DWT = DWT, Vs = speed, n = n)

@app.route('/eedi')
def eedi():
    return render_template('eedi.html')

@app.route("/result_eedi", methods = ['POST', 'GET'])
def result_eedi():
    output = request.form.to_dict()
    print(output)
    speed = float(output["speed"])
    DWT = int(output["DWT"])
    ship_type = output['type']
    fuel = output['fuel']
    n = int(output['n'])
    ymin = int(output['year-min'])
    yman = int(output['year-max'])

    years = np.linspace(ymin, yman, n)
 
    result = shiphulldes.eedi_study(speed, DWT, ship_type, fuel, years)
    result = result.to_numpy()

    return render_template('eedi.html', speed=speed, result = result, output = output, DWT=DWT)

@app.route('/air-lube')
def airlube():
    return render_template('air-lube.html')

@app.route("/result_al", methods = ['POST', 'GET'])
def result_al():
    output = request.form.to_dict()
    speed = float(output["speed"])
    DWT = int(output["DWT"])
    ship_type = output['type']
    fuel = output['fuel']

    Cf_lube = float(output['Cf_lube'])/100
    airlubearea = float(output['airlubearea'])/100

    data = shiphulldes.air_lube_est(speed, DWT, ship_type, fuel, Cf_lube, airlubearea)

    return render_template('air-lube.html', speed=speed, DWT=DWT, airlubearea=airlubearea, Cf_lube=Cf_lube)

@app.route('/econ')
def econ():
    return render_template('econ.html')

