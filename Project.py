import pandas as pd
import math
import numpy as np
from scipy.optimize import minimize

#lees data in
GRBdata = pd.read_csv('GRBdata.txt', sep=" ", header=None)
GRBdata.columns = ["name", "z", "Pbolo", "dPbolo", "Epeak", "dEpeak"]

Pantheon = pd.read_csv('GRBdata.txt', sep=" ", header=None)
Pantheon.columns = ["name", "zcmb", "zhel", "dz", "mb", "dmb"]

#enkele standaardwaardes
m_abs = -19.36
H_0 = 72 #(km/s)/Mpc
c=299792.458

#TODO: staat d_L in de juiste eenheden, nu in Mpc, moet dit pc zijn, km, lichtjaar??
#berekenen van L, hiervoor ook d_L nodig (2.141 in cursus)
L=[]
dL=[]
for index, row in GRBdata.iterrows():
    d_L =c*row['z']/H_0
    L_value = 4*math.pi*d_L**2*row['Pbolo']
    dL_value = 4*math.pi*d_L**2*row['dPbolo']
    L.append(L_value)
    dL.append(dL_value)
GRBdata['L'] = L
GRBdata['dL'] = dL


#testje: fit de data met standaard loss functie
# y = GRBdata['L'].to_numpy()
# x = []
# for index, row in GRBdata.iterrows():
#     x_value = row['Epeak']*(1+row['z'])
#     x.append(x_value)
# fit_values = np.polyfit(np.log10(x), np.log10(y), 1)
# a=fit_values[1]
# b = fit_values[0]
#print(a,b)

#TODO: moet np.log np.log10 worden? => gedaan maar is dit nu overal wel log10 => controleren
#fit de data rekening houdend met de gegeven minimalisatie, hierdoor vinden we een a en b die we kunnen gebruiken om in
# ons linair verband te steken.
y = GRBdata['L'].to_numpy()
y_error = GRBdata['dL'].to_numpy()
x = []
x_error = []
for index, row in GRBdata.iterrows():
    x_value = row['Epeak']*(1+row['z'])
    x_error_value = row['dEpeak']*(1+row['z'])
    x.append(x_value)
    x_error.append(x_error_value)
    #we nemen de log om zo lineaire regressie toe te kunnen passen, hierbij herbereken ik de error als error => log(f(x+x_error)-f(x))
#TODO is dit juist?
X = np.log10(x)
dX = np.subtract(np.add(x,x_error),np.log10(x))
Y = np.log10(y)
dY = np.subtract(np.log10(np.add(y,y_error)),np.log10(y))

def loss(X_0):
    a=X_0[0]
    b=X_0[1]
    return(sum((a+b*X[i]-Y[i])**2/(dY[i]**2+b**2*dX[i]**2) for i in range(len(X))))

X_0 = [1,1] #initiële gok, uitkomst kan hiervan afhankelijk zijn
X_0 = np.array(X_0)
X_min = minimize(loss, X_0).x
a=X_min[0]
b=X_min[1]

#test schatting
#TODO schatting van L geeft verkeerde waarde :'(
i=10
L_est = math.pow(10,(a+b*math.log10((1+GRBdata['z'][i])*(GRBdata['Epeak'][i]))))
print(GRBdata['L'][i], L_est)

