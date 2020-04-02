# -*- coding: utf-8 -*-
# 
# Práctica de laboratorio de infiltrómetro de doble anillo.
# ------------------------------------------------------------------------------
# Escrito por : Michael Heredia Pérez
# Fecha       : Nov/2019
# e-mail      : mherediap@unal.edu.co
# Universidad Nacional de Colombia, Manizales field.
# ------------------------------------------------------------------------------
# Lo siguiente es la optimización de la función de inflitración para el esanyo 
# realizado
# ------------------------------------------------------------------------------

# %% Importo librerías
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit     # Con esta realizo el ajuste
#from openpyxl import load_workbook      # Para sobreescribir el excel

# %% Defino variables

a_ani_int   = np.pi*30**2/4             # cm^2, área del anillo interior
a_esp_anu   = np.pi*60**2/4 - a_ani_int # cm^2, área del espacio anular
lam_ani_int = 31                        # cm, lámina en el anillo interior
lam_esp_anu = 33                        # cm, lámina en el espacio anular.

# importo series de tiempo.
df = pd.read_excel('info_de_campo.xlsx', sheet_name = None) # Datos de campo
tub_g = df['Hoja1']['tubo_grande'].to_numpy()               # Tubo grande
tub_p = df['Hoja1']['tubo_pequeno'].to_numpy()              # Tubo pequeño
tt    = df['Hoja1']['tiempo'].to_numpy()                    # Los minutos 
                                                            # transcurridos

# %% Se calcula la variación en la velocidad de infiltración para el anillo
# anillo interior medido por el tubo pequeño y para el espacio anular medido
# por el tubo grande.                                                            

vol_ani_int = []     # Almacena el volumen inflitrado en el anillo interior
                     # en función de tt
vol_esp_anu = []     # Almacena el volumen infiltrado en el espacio anular 
                     # en función de tt

mediciones = len(tt) # puntos que se tomaron 

for i in range(mediciones-1):
    int_i = tub_p[i]   # el dato inicial del anillo interior
    int_j = tub_p[i+1] # el dato final del anillo interior
    anu_i = tub_g[i]   # el dato inicial del espacio anular
    anu_j = tub_g[i+1] # el dato final del espacio anular
    
    # El volumen
    vol_int = np.abs(int_i - int_j)*a_ani_int # cm^3
    vol_anu = np.abs(anu_i - anu_j)*a_esp_anu # cm^3
    
    # Los añado a la lista
    vol_ani_int.append(vol_int)
    vol_esp_anu.append(vol_anu)

# paso a arrays de numpy
vol_ani_int = np.array(vol_ani_int)
vol_esp_anu = np.array(vol_esp_anu)

# Elimino el tiempo cero del ensayo
tt = np.delete(tt, 0)

# La velocidad de infiltración
infil_ani_int = vol_ani_int/(a_ani_int*tt)*10*60  # mm/h, velocidad de 
                                                  # infiltración en anillo                     
infil_esp_anu = vol_esp_anu/(a_esp_anu*tt)*10*60  # mm/h, velocidad de 
                                                  # infiltracion en espacio

# %% Grafico estos resultados

plt.figure(1)
plt.plot(tt, infil_ani_int, 'mo', label = 'Anillo interior')
plt.plot(tt, infil_esp_anu, 'bo', label = 'Espacio anular')
plt.legend(loc = 'best')
plt.xlabel(r'Tiempo transcurrido $[min]$')
plt.ylabel(r'Tasa de infiltración $[mm/h]$')
plt.title('Datos experimentales ensayo del infiltrómetro de doble anillo')
plt.savefig('datos_experimentales.jpg')
plt.show()

# %% Defino las ecuación a optimizar

# La ecuación de Kostiakov
#def kostiacov(t, a, b):
#    return (a/(1-b))*t**(1-b)

# La ecuación de Horton
def horton_int(t, g):       # Para el anillo interior
    return infil_ani_int[-1]+(infil_ani_int[0]-infil_ani_int[-1])*np.exp(-g*t)

def horton_esp(t, g):       # Para el espacio anular
    return infil_esp_anu[-1]+(infil_esp_anu[0]-infil_esp_anu[-1])*np.exp(-g*t)

# La ecuación de Mezencev
def mezencev_int(t, a, b):
    return infil_ani_int[-1] + a*t**(-b)

def mezencev_esp(t, a, b):
    return infil_esp_anu[-1] + a*t**(-b)

# La ecuación de Philip
def philip(t, A, B):
    return 0.5*A*t**(-0.5)+B

# La ecuación de Holtan no se tuvo en cuenta

# %% Realizo la optimización del ANILLO INTERIOR
    
params_hor_int, pcov_hor_int = curve_fit(horton_int  , tt, infil_ani_int)
params_mez_int, pcov_mez_int = curve_fit(mezencev_int, tt, infil_ani_int) 
params_phi_int, pcov_phi_int = curve_fit(philip      , tt, infil_ani_int) 
#params_kos_int, pcov_kos_int = curve_fit(kostiacov   , tt, infil_ani_int)

# se construyen los polinomios
fig_hor_int = horton_int(tt, *params_hor_int)
fig_mez_int = mezencev_int(tt, *params_mez_int)
fig_phi_int = philip(tt, *params_phi_int)

# se grafican los resultados
plt.figure(2)
plt.plot(tt, infil_ani_int, 'mo', label = 'Datos')
plt.plot(tt, fig_hor_int, 'y-', label = 'Horton')
plt.plot(tt, fig_mez_int, 'r-', label = 'Mezencev')
plt.plot(tt, fig_phi_int, 'g-', label = 'Philip')
plt.legend(loc = 'best')
plt.xlabel(r'Tiempo transcurrido $[min]$')
plt.ylabel(r'Tasa de infiltración $[mm/h]$')
plt.title('Ajuste a la tasa de infiltración del anillo interior')
plt.savefig('anillo_interior.jpg')
plt.show()
  
# %% Realizo la optimización del ESPACIO ANULAR
    
params_hor_esp, pcov_hor_esp = curve_fit(horton_esp  , tt, infil_esp_anu)
params_mez_esp, pcov_mez_esp = curve_fit(mezencev_esp, tt, infil_esp_anu) 
params_phi_esp, pcov_phi_esp = curve_fit(philip      , tt, infil_esp_anu) 
#params_kos_int, pcov_kos_int = curve_fit(kostiacov   , tt, infil_ani_int)

# se construyen los polinomios
fig_hor_esp = horton_esp(tt, *params_hor_esp)
fig_mez_esp = mezencev_int(tt, *params_mez_esp)
fig_phi_esp = philip(tt, *params_phi_esp)

# se grafican los resultados
plt.figure(3)
plt.plot(tt, infil_esp_anu, 'bo', label = 'Datos')
plt.plot(tt, fig_hor_esp, 'y-', label = 'Horton')
plt.plot(tt, fig_mez_esp, 'r-', label = 'Mezencev')
plt.plot(tt, fig_phi_esp, 'g-', label = 'Philip')
plt.legend(loc = 'best')
plt.xlabel(r'Tiempo transcurrido $[min]$')
plt.ylabel(r'Tasa de infiltración $[mm/h]$')
plt.title('Ajuste a la tasa de infiltración del espacio anular')
plt.savefig('espacio_anular.jpg')
plt.show()

# %% Calculo el Error cuadrático medio de cada ecuación (mean square error mse)

# Para el anillo interior
mse_hor_int = np.square(np.subtract(infil_ani_int, fig_hor_int)).mean()
mse_mez_int = np.square(np.subtract(infil_ani_int, fig_mez_int)).mean()
mse_phi_int = np.square(np.subtract(infil_ani_int, fig_phi_int)).mean()

# Para el espacio anular
mse_hor_esp = np.square(np.subtract(infil_esp_anu, fig_hor_esp)).mean()
mse_mez_esp = np.square(np.subtract(infil_esp_anu, fig_mez_esp)).mean()
mse_phi_esp = np.square(np.subtract(infil_esp_anu, fig_phi_esp)).mean()

# Selecciono la funciónq ue me da menor error
eq = ['Horton', 'Mezencev', 'Philip']    # Las ecuaciones empleadas

# Agrupo los errores en un diccionario
mse = {'Anillo interior': [mse_hor_int, mse_mez_int, mse_phi_int],
       'Espacio Anular' : [mse_hor_esp, mse_mez_esp, mse_phi_esp]}

mse_df = pd.DataFrame(data = mse,        # Se transforma todo en una tabla
                      index = eq)        # Se organiz el índice con las
                                         # funciones empleadas 

#Imprimo las desviaciones estandar
print('\nSegún el Error Cuadrático Medio (MSE por sus siglas en inglés): ')
print()
print(mse_df)

# %% Pensando que la mejor es la de Horton, se grafican estas dos para dejar
# registro

plt.figure()
plt.plot(tt, infil_ani_int, 'mo', label = 'Anillo interior')
plt.plot(tt, infil_esp_anu, 'bo', label = 'Espacio anular')
plt.plot(tt, fig_mez_int, '-', label = 'Mezencev anillo interior')
plt.plot(tt, fig_mez_esp, '-', label = 'Mezencev espacio anular')
plt.legend(loc = 'best')
plt.xlabel(r'Tiempo transcurrido $[min]$')
plt.ylabel(r'Tasa de infiltración $[mm/h]$')
plt.title('La ecuación de Mezencev para el anillo interno y el espacio anular')
plt.savefig('vistazo_mezencev.jpg')
plt.show()

# %% Evalúo el k (la constante de infiltración)
k_int = mezencev_int(1000, *params_mez_int)
k_esp = mezencev_esp(1000, *params_mez_esp)

print('La "k" según el anillo interior será: ', k_int.round(3))

print('La "k" según el espacio anular será: ', k_esp.round(3))

# Fin del código.