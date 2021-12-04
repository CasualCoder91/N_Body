import pandas as pd
import numpy as np
from uncertainties import ufloat

import config

df = pd.read_excel(config.output_base_path+r'\25_observations.xlsx','extinction', usecols = 'A:AV')

def output_mean_error(row, mass, angle):
    raw_data = df.query('Mass=='+str(mass)+'&'+'Angle=='+str(angle))[row]
    mean_error = ufloat(raw_data.mean(), raw_data.std())
    digits = np.abs(int(np.log10(abs(mean_error.s))-2))
    mean = np.round(mean_error.n,digits)
    error = np.round(mean_error.s,digits)
    format_string = '{:.'+str(digits)+'f}'
    print(format_string.format(mean),',',format_string.format(error))


def Vel2D_error():
    
    masses = [640] #df['Mass'].unique()
    angles = df['Angle'].unique()
    rows = ['avgVel2DCluster']#'avgVel2DCluster','avgVel2DFS','disp2DCluster','disp2DFS']

    for angle in angles:
        for mass in masses:
            for row in rows:
                output_mean_error(row,mass,angle)


def Precision_error(mass_range,simulated=False):

    masses = df['Mass'].unique()
    angles = df['Angle'].unique()

    for angle in angles:
        for mass in masses:
            query = 'Mass=='+str(mass)+'&'+'Angle=='+str(angle)
            if simulated:
                TPdf = df.query(query)['TP '+mass_range]
            else:
                TPdf = df.query(query)['CTP '+mass_range]
            TP = ufloat(TPdf.mean(), TPdf.std())

            if simulated:
                FPdf = df.query(query)['FP '+mass_range]
                FP = ufloat(FPdf.mean(), FPdf.std())
            else:
                CFPdf = df.query(query)['CFP '+mass_range]
                UPdf = df.query(query)['UP '+mass_range]
                FP = ufloat(CFPdf.mean(), CFPdf.std()) + ufloat(UPdf.mean(), UPdf.std())

            P = TP/(TP+FP)
            #print(uncertainty.s)
            if P.s == 0:
                digits = 10
            else:
                digits = np.abs(int(np.log10(abs(P.s))-2))

            mean = np.round(P.n,digits)
            error = np.round(P.s,digits)
            format_string = '{:.'+str(digits)+'f}'
            print(format_string.format(mean),',',format_string.format(error))

    return