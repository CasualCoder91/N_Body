import pandas as pd
import numpy as np
from uncertainties import ufloat

import config

df = pd.read_excel(config.output_base_path+r'\25_observations.xlsx','extinction', usecols = 'A:BM')


def output_mean_error(x,new_line=True):

    if x.s == 0:
        digits = 10
    else:
        digits = np.abs(int(np.log10(abs(x.s))-2))

    mean = np.round(x.n,digits)
    error = np.round(x.s,digits)
    format_string = '{:.'+str(digits)+'f}'
    if(new_line):
        print(format_string.format(mean),format_string.format(error))
    else:
        print(format_string.format(mean),format_string.format(error), end =" ")

def Vel2D_error():
    
    masses = df['Mass'].unique()
    angles = df['Angle'].unique()
    rows = ['avgVel2DFS']#'avgVel2DCluster','avgVel2DFS','disp2DCluster','disp2DFS']

    for angle in angles:
        for mass in masses:
            for row in rows:
                raw_data = df.query('Mass=='+str(mass)+'&'+'Angle=='+str(angle))[row]
                mean_error = ufloat(raw_data.mean(), raw_data.std())
                output_mean_error(mean_error)


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
            output_mean_error(P)

    return

def F1_error(mass_range):

    masses = df['Mass'].unique()
    angles = df['Angle'].unique()

    for angle in angles:
        for mass in masses:
            UPdf = df.query('Mass=='+str(mass)+'&'+'Angle=='+str(angle))['UP '+mass_range]
            UP = ufloat(UPdf.mean(), UPdf.std())
            CFPdf = df.query('Mass=='+str(mass)+'&'+'Angle=='+str(angle))['CFP '+mass_range]
            CFP = ufloat(CFPdf.mean(), CFPdf.std())
            CTPdf = df.query('Mass=='+str(mass)+'&'+'Angle=='+str(angle))['CTP '+mass_range]
            CTP = ufloat(CTPdf.mean(), CTPdf.std())
            CFNdf = df.query('Mass=='+str(mass)+'&'+'Angle=='+str(angle))['CFN '+mass_range]
            CFN = ufloat(CFNdf.mean(), CFNdf.std())
            F1 = CTP/(CTP+0.5*(CFP+UP+CFN))
            output_mean_error(F1)

def column_group_error(group):
    masses = df['Mass'].unique()
    angles = df['Angle'].unique()

    mass_ranges = ['Tot','> 2','2 - 0.5','0.5 - 0.08']

    for angle in angles:
        for mass in masses:
            query = 'Mass=='+str(mass)+'&'+'Angle=='+str(angle)
            for mass_range in mass_ranges:
                cdf = df.query(query)[group+mass_range]
                c = ufloat(cdf.mean(), cdf.std())
                output_mean_error(c,False)
            print('')
