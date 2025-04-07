import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fluidproperties as fpo
import BeggsandBrill as BBe

# === Título
st.title("Curva VLP - Vertical Lift Performance")

# === Layout con dos columnas
col1, col2 = st.columns([1, 1])

with col1:
    st.header("Parámetros de entrada")
    oil_rate = st.number_input("Caudal de petróleo (m3/día)", value=3.0)
    water_rate = st.number_input("Caudal de agua (m3/día)", value=30.0)
    GOR = st.number_input("Razón gas/petróleo (m3/m3)", value=200.0)
    gas_grav = st.number_input("Densidad específica del gas", value=0.65)
    oil_grav = st.number_input("Gravedad API del petróleo", value=35.0)
    wtr_grav = st.number_input("Gravedad específica del agua", value=1.07)
    diameter = st.number_input("Diámetro interno de la tubería (pulgadas)", value=2.441)
    angle = st.number_input("Ángulo del pozo (grados)", value=90.0)
    thp = st.number_input("Presión en cabeza de pozo (psia)", value=150.0)
    tht = st.number_input("Temperatura en cabeza de pozo (°C)", value=20.0)
    twf = st.number_input("Temperatura en fondo de pozo (°C)", value=100.0)
    depth = st.number_input("Profundidad del pozo (metros)", value=2000.0)
    roughness = st.number_input("Rugosidad de la tubería (pulgadas)", value=0.0006)
    Psep = st.number_input("Presión del separador (psia)", value=114.7)
    Tsep = st.number_input("Temperatura del separador (°C)", value=50.0)

# Conversión de unidades
stb_to_m3 = 1 / 6.28981  # Factor de conversión

oil_rate_stb = oil_rate * 6.28981
water_rate_stb = water_rate * 6.28981
GOR = GOR * 5.61458
tht = tht * 9 / 5 + 32
twf = twf * 9 / 5 + 32
depth = depth * 3.28084
Tsep = Tsep * 9 / 5 + 32

# === Cálculos
WC = water_rate_stb / (oil_rate_stb + water_rate_stb)

def temperature_gradient(T1, T2, Depth):
    return abs(T1 - T2) / Depth if Depth else 0

t_grad = temperature_gradient(tht, twf, depth)
depths = np.linspace(0, depth, 50)
temps = tht + t_grad * depths

def pressure_traverse(oil_rate):
    pressure_list = []
    dpdz_list = []
    
    for i in range(len(depths)):
        if i == 0:
            pressure_list.append(thp)
        else:
            dz = (depths[i] - depths[i-1])
            pressure = pressure_list[i-1] + dz * dpdz_list[i-1]
            pressure_list.append(round(pressure, 3))

        dpdz_step = BBe.beggsandbrill(pressure_list[i], temps[i], oil_rate, WC, GOR, gas_grav, oil_grav, wtr_grav,
                                      diameter, angle, roughness, Psep, Tsep)
        dpdz_list.append(round(dpdz_step, 5))
    
    return pressure_list, dpdz_list

def VLP(rates):
    bhp_list = []
    for q in rates:
        p, _ = pressure_traverse(q)
        bhp_list.append(round(p[-1], 3))
    return bhp_list

# === Cálculo de la curva VLP
rate_array_stb = np.linspace(100, 5000, 50)
bhp_array = VLP(rate_array_stb)

# Conversión a m3/d para visualización
rate_array_m3 = rate_array_stb * stb_to_m3

with col2:
    st.header("Gráfico de la curva VLP")
    fig, ax = plt.subplots()
    ax.plot(rate_array_m3, bhp_array)
    ax.set_xlabel("Caudal total (m3/día)")
    ax.set_ylabel("Presión en fondo (psi)")
    ax.set_title("Curva VLP")
    st.pyplot(fig)

# === Mostrar tabla opcional debajo del gráfico
st.write("\n\n")  # Espaciado para mejorar la visualización
if st.checkbox("Mostrar tabla de resultados"):
    with col2:
        ratevspress = pd.DataFrame({"Caudal (m3/día)": rate_array_m3, "Presión en fondo (psi)": bhp_array})
        st.dataframe(ratevspress)
