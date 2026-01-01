from ansys.aedt.core.modeler.cad.polylines import PolylineSegment
from ansys.aedt.core import Desktop
from ansys.aedt.core import Maxwell3d
from openpyxl import Workbook, load_workbook
import numpy as np
import os
import math
import time

#Guardar Resultados
archivo = "Resultados_Pruebas_Eddy_1.xlsx"
archivo_txt = "Objetivos.txt"

#desktop = Desktop(version="2025.1", new_desktop=False)
m3d = Maxwell3d(projectname="RectCoilProject4", designname="SqrCoilDesign3d", solution_type="EddyCurrent")
m3d.modeler.model_units = "mm"
m3d.modeler.delete("PL_line")
m3d.modeler.delete("PL_line2")
m3d.modeler.delete("AirBox1")
m3d.modeler.delete("AirBox2")
m3d.delete_setup("My_Setup")

################################## PARÁMETROS BOBINAS ####################################

dist = 50 # Distancia en el eje z entre bobinas
num_seg = 0 #Número de segmentos Polilinea
frequency = 500000 #Hz Frecuencia 

# Parámetros principales PRIMERA BOBINA
H = 50  # Altura (eje Y)
W = 50  # Ancho (eje X)
h = H/2 #Centro bobina
w = W/2 #Centro bobina
n = 1 # Número de capas
num_e = 4  #Número de espiras


wire_d = 0.81 # Calibre 20
extension_length = 5 
extension_gap = 8
z_offset = 1.4  # Separación entre capas
spc = 3        # Separación hacia adentro por cada espira


# Parámetros principales segunda bobina
H2 = 50  # Altura (eje Y)
W2 = 50  # Ancho (eje X)
h2 = H2/2
w2 = W2/2

n2 = 1 # Número de capas
num_e2 = 4 #Numero de espiras


wire_d2 = 0.81  # Calibre 22 (0.64)
extension_length2 = 5 
extension_gap2 = 8
z_offset2 = 1.4  # Separación entre capas
spc2 = 3        # Separación hacia adentro por cada espira

#################################################### TIEMPO #################################################

inicio = time.time() #Inicio tiempo

################################## PRIMERA  ##############################################

coil_points = []

# Terminal de entrada
coil_points.append([W + extension_length - w, H/2 - extension_gap/2 -h, 0])
coil_points.append([W - w, H/2 - extension_gap/2 - h, 0])


for capa in range(num_e):
    offset = capa * spc

    # Alterna entre orden ascendente y descendente para las espiras
    rango = range(n) if capa % 2 == 0 else range(n - 1, -1, -1)

    for i in rango:
        z = i * z_offset
        coil_points.append([W - capa * spc - w, 0 + offset -h, z])
        coil_points.append([0 + offset -w, 0 + offset -h, z])
        coil_points.append([0 + offset - w, H - offset -h, z])
        coil_points.append([W - capa * spc -w, H - offset -h, z])


# Último tramo hacia el terminal de salida

if num_e % 2 == 0:
    a = 1
else: a = -1

coil_points.append([W-capa*spc - w, H/2 + extension_gap/2 -h, z - a*z_offset ]) #a * (z + wire_d/2)
coil_points.append([W + extension_length - w, H/2 + extension_gap/2 - h, z - a*z_offset]) #a * (z + wire_d/2)


# Crear la polilínea en Maxwell
rectangular_coil = m3d.modeler.create_polyline(
    coil_points, 
    name="PL_line", 
    material="Copper", 
    xsection_type="Circle", 
    xsection_width=str(wire_d)+"mm",
    xsection_num_seg= num_seg
)

air_box = m3d.modeler.create_box(
    position=[-5 -w, -5 -h, -5],
    sizes=[W + extension_length + 5 , H + 10, (n*wire_d) + 5 + dist/2],
    name="AirBox1",
    material="Air"
)

######################################### SEGUNDA ##############################################

coil_points2 = []

# Terminal de entrada2
coil_points2.append([W2 + extension_length2 - w2, H2/2 - extension_gap2/2 - h2, dist])
coil_points2.append([W2 - w2, H2/2 - extension_gap2/2 - h2, dist])


for capa2 in range(num_e2):
    offset2 = capa2 * spc2

    # Alterna entre orden ascendente y descendente para las espiras
    rango2 = range(n2) if capa2 % 2 == 0 else range(n2 - 1, -1, -1)

    for i in rango2:
        z2 = i * z_offset2
        coil_points2.append([W2 - capa2 * spc2 - w2, 0 + offset2 -h2, z2+dist])
        coil_points2.append([0 + offset2 - w2, 0 + offset2 - h2, z2+dist])
        coil_points2.append([0 + offset2 - w2, H2 - offset2 - h2, z2+dist])
        coil_points2.append([W2 - capa2 * spc2 - w2, H2 - offset2 - h2, z2+dist])


# Último tramo hacia el terminal de salida

if num_e2 % 2 == 0:
    a2 = 1
else: a2 = -1

# Último tramo hacia el terminal de salida
coil_points2.append([W2-capa2*spc2 - w2, H2/2 + extension_gap2/2 - h2, z2 + dist - a2*z_offset2])
coil_points2.append([W2 + extension_length2 - w2, H2/2 + extension_gap2/2 - h2, z2 + dist - a2*z_offset2])



rectangular_coil2 = m3d.modeler.create_polyline(
    coil_points2,
    name="PL_line2",
    material="Copper",
    xsection_type="Circle",         # Cambia de "Circle" a "Polygon"
    xsection_width=str(wire_d2)+"mm",
    xsection_num_seg= num_seg
)


air_box = m3d.modeler.create_box(
    position=[-5 -w2, -5 - h2, (n*wire_d) + dist/2 ],
    sizes=[W2 + extension_length + 5 , H2 + 10, (n*wire_d) + dist/2 + 5],
    name="AirBox2",
    material="Air"
)


united_air = m3d.modeler.unite(
    ["AirBox1", "AirBox2"]
)

################################################## CORRIENTES #################################################

faces= rectangular_coil.faces

faces_ordenadas = sorted(faces, key=lambda f: f.area)

current1 = m3d.assign_current(faces_ordenadas[0], amplitude="3A", swap_direction=False, name= "term1", phase="0deg")
current2 = m3d.assign_current(faces_ordenadas[1], amplitude="3A", swap_direction=True, name= "term2", phase="0deg")

#---------------------------------------------------------------------------------------------------------------

faces2= rectangular_coil2.faces

faces_ordenadas2 = sorted(faces2, key=lambda f: f.area)

current3 = m3d.assign_current(faces_ordenadas2[0], amplitude="3A", swap_direction=False, name= "term3", phase="0deg" )
current4 = m3d.assign_current(faces_ordenadas2[1], amplitude="3A", swap_direction=True, name= "term4", phase="0deg")


######################################################## SIMULACION #################################################  
selection = ["term1", "term3"]
group_sources = {"Group1_Test": ["term1", "term3"]}

L = m3d.assign_matrix(selection, matrix_name="Matrix1", turns=None, return_path=None, group_sources=None, branches=None)


m3d.create_setup(name="My_Setup",
    setup_type="EddyCurrent",
    MaximumPasses=1,
    PercentError=2,
    Frequency=str(frequency)+"Hz",
)

#---------------------------------------------------------------------------------------------------------------

m3d.analyze_setup("My_Setup")

################################################ RESULTADOS ####################################################

inductance = ["Matrix1.L(term3, term1)","Matrix1.L(term1, term3)","Matrix1.L(term3, term3)","Matrix1.L(term1, term1)"]

inductance_report = m3d.post.create_report(
    expressions=inductance,
    setup_sweep_name="My_Setup",
    primary_sweep_variable="Matrix1.L(term1,term1)",
    report_category="L",
    plot_type="Data Table",
    plot_name="InductanceReport"
)

inductance_data1 = m3d.post.get_solution_data("Matrix1.L(term1, term1)", "My_Setup")
inductance_data2 = m3d.post.get_solution_data("Matrix1.L(term3, term3)", "My_Setup")

L_value1 = inductance_data1.data_real()[0]
L_value2 = inductance_data2.data_real()[0]  

print("Inductancia1:", L_value1, "nH", "---","Inductancia2:", L_value2, "nH")

##-----------------------------------------------------------------------------------------------------------------------------------------------------

coupling = ["Matrix1.CplCoef(term1, term1)","Matrix1.CplCoef(term3, term1)","Matrix1.CplCoef(term1, term3)","Matrix1.CplCoef(term3, term3)"]

coupling_coef_report = m3d.post.create_report(
    expressions=coupling,
    setup_sweep_name="My_Setup",
    primary_sweep_variable="Matrix1.CplCoef(term1,term1)",
    report_category="Coupling Coeff",
    plot_type="Data Table",
    plot_name="CouplingCoeffReport"
)

coupling_data = m3d.post.get_solution_data("Matrix1.CplCoef(term3, term1)", "My_Setup")

CouplingC = abs(coupling_data.data_real()[0])

print("Coeficiente de Acoplamiento:", CouplingC)
##----------------------------------------------------------------------------------------------------------------------------------------------------

resistance = ["Matrix1.R(term1, term1)","Matrix1.R(term3, term1)","Matrix1.R(term1, term3)","Matrix1.R(term3, term3)"]

resistance_report = m3d.post.create_report(
    expressions=resistance,
    setup_sweep_name="My_Setup",
    primary_sweep_variable="Matrix1.R(term1,term1)",
    report_category="R",
    plot_type="Data Table",
    plot_name="ResistanceReport"
)

resistance_data1 = m3d.post.get_solution_data("Matrix1.R(term1, term1)", "My_Setup")
resistance_data2 = m3d.post.get_solution_data("Matrix1.R(term3, term3)", "My_Setup")

R_value1 = resistance_data1.data_real()[0]
R_value2 = resistance_data2.data_real()[0]  

print("Resistencia1:", R_value1*1000, "mOhm", "---","Resistencia2:", R_value2*1000, "mOhm")

###----------------------------------------------------------------------------------------------------------------------------------------------------

impedance = ["Matrix1.Z(term1, term1)","Matrix1.Z(term3, term1)","Matrix1.Z(term1, term3)","Matrix1.Z(term3, term3)"]

impedance_report = m3d.post.create_report(
    expressions=impedance,
    setup_sweep_name="My_Setup",
    primary_sweep_variable="Matrix1.Z(term1,term1)",
    report_category="Z",
    plot_type="Data Table",
    plot_name="ImpedanceReport"
)

impedance_data1 = m3d.post.get_solution_data("Matrix1.Z(term1, term1)", "My_Setup")
impedance_data2 = m3d.post.get_solution_data("Matrix1.Z(term3, term3)", "My_Setup")

Z_value1 = impedance_data1.data_real()[0]
Z_value2 = impedance_data2.data_real()[0]
Z_value1_imag = impedance_data1.data_imag()[0]
Z_value2_imag = impedance_data2.data_imag()[0]  

print("Impedancia1: ", Z_value1*1000,"+", Z_value1_imag*1000, "i", "mOhm")
print("Impedancia2: ", Z_value2*1000,"+", Z_value2_imag*1000, "i", "mOhm")

###----------------------------------------------------------------------------------------------------------------------------------------------------

MagFlux = ["Matrix1.MagFlux(term1)","Matrix1.MagFlux(term3)"]

impedance_report = m3d.post.create_report(
    expressions=MagFlux,
    setup_sweep_name="My_Setup",
    primary_sweep_variable="Matrix1.MagFlux(term1)",
    report_category="MagFlux",
    plot_type="Data Table",
    plot_name="MagFluxReport"
)

MagFlux_data1 = m3d.post.get_solution_data("Matrix1.MagFlux(term1)", "My_Setup")
MagFlux_data2 = m3d.post.get_solution_data("Matrix1.MagFlux(term3)", "My_Setup")

MagFlux_value1 = MagFlux_data1.data_real()[0]
MagFlux_value2 = MagFlux_data2.data_real()[0]
MagFlux_value1_imag = MagFlux_data1.data_imag()[0]
MagFlux_value2_imag = MagFlux_data2.data_imag()[0]  

print("Flujo Magnético1: ", MagFlux_value1,"+", MagFlux_value1_imag, "i", "Wb")
print("Flujo Magnético2: ", MagFlux_value2,"+", MagFlux_value2_imag, "i", "Wb")

###----------------------------------------------------------------------------------------------------------------------------------------------------

M_calculated = CouplingC * math.sqrt(L_value1 * L_value2)

print("Inductancia Mutua:", M_calculated, "nH")

###----------------------------------------------------------------------------------------------------------------------------------------------------

#Calculo de Leq

omega = 2 * np.pi * frequency

X1 = Z_value1_imag
X2 = Z_value2_imag

Leq1 = X1 / omega
Leq2 = X2 / omega

print("InductanciaEq1:", Leq1 * 1e9, "nH", "---","InductanciaEq2:", Leq2 * 1e9, "nH")

###----------------------------------------------------------------------------------------------------------------------------------------------------

#Calculo de Qeq

Q1eq = Z_value1_imag/Z_value1

Q2eq = Z_value2_imag/Z_value2

print("FactorCalidad1:", Q1eq, "---","FactorCalidad2:", Q2eq)

###----------------------------------------------------------------------------------------------------------------------------------------------------

#Calculo de Eficiencia

eta = (CouplingC**2 * Q1eq * Q2eq) / (1 + CouplingC**2 * Q1eq * Q2eq)

print("Eficiencia:", eta * 100, "%")

################################################## TIEMPO ######################################################

fin = time.time()

tiempo = fin - inicio
print(f"Tiempo de ejecución: {tiempo:.4f} segundos")

################################################## LIBERAR DESKTOP #################################################

m3d.save_project()
m3d.release_desktop(close_projects=False, close_desktop=False)
