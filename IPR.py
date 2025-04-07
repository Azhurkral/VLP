import matplotlib.pyplot as plt 
import streamlit as st
import math
def ipr(pr,pb,ptest,qtest,re=1000,rw=0.25,h=100,ko=100,mu=100, Bo=1):
    """
    Some values like well radius or formation oil factor are entered default.
    First, we need to find the pressure vs flow rate for the saturated reservoir conditions. 
    Then, these are found for the undersaturated reservoir conditions under 2 different scenarios.
      1) If test pressure is equal to or greater than bubble point
      2) If the test pressure is below the bubble point
    """
    if pb>=pr:
       list_q=[]
       list_pwf2=[]
       ylim=pr
       j_darcy=0.00708*ko*h/(mu*Bo*math.log(0.472*re/rw))
       print("Productivity index(J) is {} STB/D-psi".format(round(j_darcy,4)))
       qob=0
       for i in range(round(pr,-2),-100,-100):
           qomax=qtest/(1-0.2*ptest/pr-0.8*(ptest/pr)**2)
           qo=qomax*(1-0.2*(i/pr)-0.8*(i/pr)**2)
           list_pwf2.append(i)
           list_q.append(round(qo))
       print(list_q,list_pwf2)
       
    else:   
        #The case for undersaturated reservoirs where Pb is less than Pr
        ylim=pr    
        list_pwf2=[]
        list_q=[]
        for i in range(round(pr,-1),-100,-100):
            list_pwf2.append(i)
            
            if ptest>=pb:
                #When the entered test pressure is below bubble point
                j1=qtest/(pr-ptest)       
                qob=j1*(pr-pb)
                
                if i>=pb:
                    #When the loop pressure value of i is greater than bubble point                    
                    qo=j1*(pr-i)
                    list_q.append(round(qo))
                else:
                    #When the loop pressure value of i is less than bubble point
                    qo=qob+(j1*pb*(1-0.2*i/pb-0.8*(i/pb)**2)/1.8)
                    list_q.append(round(qo))
                    
            else:
                #When the entered test pressure is above bubble point
                j1=((qtest/((pr-pb)+pb/1.8*(1-0.2*ptest/pr-0.8*(ptest/pr)**2))))  
                qob=((qtest/((pr-pb)+pb/1.8*(1-0.2*ptest/pr-0.8*(ptest/pr)**2))*(pr-pb)))
                if i>=pb:
                    qo=j1*(pr-i)
                    list_q.append(round(qo))
                else:
                    qo=qob+(j1*pb/1.8*(1-0.2*i/pb-0.8*(i/pb)**2))
                    list_q.append(round(qo))  
    
    plt.plot(list_q, list_pwf2, marker="x",ms=5)
    plt.ylim(0,ylim)
    plt.xlabel("Flow Rate, STB/D")  
    plt.ylabel("Flowing Bottomhole Pressure, psi")
    plt.scatter(qob, pb, marker="o",color="r")
    plt.annotate("Bubble Point", (qob, pb))
    plt.plot(0, pr, marker="o",color="b")
    plt.title(label="Inflow Performance Relation Curve") 

# Título de la app
st.title("Curva de IPR - Análisis de Productividad")

# Entradas interactivas
pr = st.slider("Presión del reservorio (psi)", 1000, 5000, 3500)
pb = st.slider("Presión de burbuja (psi)", 1000, pr, 2000)
ptest = st.slider("Presión de prueba (psi)", 1000, 4000, 1500)
qtest = st.slider("Caudal de prueba (STB/d)", 100, 5000, 500)

ipr(pr,pb,ptest,qtest) #Presión de reservorio, Presión de burbuja, Presión de ensayo, Caudal de ensayo
st.pyplot(plt.gcf())
