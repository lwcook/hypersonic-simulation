import json
import pdb

import numpy as hp
import matplotlib.pyplot as plt

import hypersonicsimulation.aerodynamics as aero
import hypersonicsimulation.geometry as geom
import hypersonicsimulation.vehicle as veh
import hypersonicsimulation.plotting as hsplot

def main():

    with open('validation/CCY_validation_data.txt', 'r') as f:
        lines = f.readlines()
        jstring = ''
        for line in lines:
            first_char = line.split()[0].lower()
            if first_char != '%':
                jstring += line.strip('\r').strip('\n').strip('\r')
        jdict = json.loads(jstring)

    ccy = veh.ConeCylinder(Nstrips=30, Npanels=30)

    fig, (ax1, ax2) = plt.subplots(1, 2)
    Ms = [6, 7, 8, 9]
    colors = [hsplot.blue, hsplot.red, hsplot.green, hsplot.grey]
    for iM, M in enumerate(Ms):
        print('M: ', M)
        lkey = 'M_' + str(M) + '_Lift'
        data = jdict[lkey]
        lalphas = [d[0] for d in data]
        Lifts_vlid = [d[1] for d in data]
        Lifts_pred = []
        for alpha in lalphas:
            adyn = aero.AeroModel(M=M, alpha=alpha, dynamic_pressure=50000)
            Cl, Cd = adyn.analyze_geometry(ccy.geometry, coeffs=True)
            Lifts_pred.append(Cl)

        dkey = 'M_' + str(M) + '_Drag'
        data = jdict[dkey]
        dalphas = [d[0] for d in data]
        Drags_vlid = [d[1] for d in data]
        Drags_pred = []
        for alpha in dalphas:
            adyn = aero.AeroModel(M=M, alpha=alpha, dynamic_pressure=50000)
            Cl, Cd = adyn.analyze_geometry(ccy.geometry, coeffs=True)
            Drags_pred.append(Cd)

        ax1.plot(lalphas, Lifts_vlid, c=colors[iM], linestyle='dashed',
            label='M'+str(M))
        ax1.plot(lalphas, Lifts_pred, c=colors[iM], linestyle='solid')
        ax2.plot(dalphas, Drags_vlid, c=colors[iM], linestyle='dashed',
            label='M'+str(M))
        ax2.plot(dalphas, Drags_pred, c=colors[iM], linestyle='solid')

    plt.show()


    hsplot.plot_geometry(ccy.geometry)
    plt.show()



if __name__ == "__main__":
    main()
