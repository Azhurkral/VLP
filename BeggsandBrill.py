# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 20:35:36 2023

@author: emink
"""

import math
import fluidproperties as pvt
import random
import matplotlib.pyplot as plt


def beggsandbrill(P,T,liquid_rate, WC, GOR, gas_grav, oil_grav, wtr_grav, diameter, angle, roughness, Psep,Tsep):
    #P              Pressure, psia
    #T              Temperature, °F
    #liquid_rate    Oil flowrate, stb/D
    #WC             Water Cut, %
    #GOR            Producing gas-oil ratio, scf/stb
    #gas_grav       Gas specific gravity
    #oil_grav       Oil gravity, API
    #wtr_grav       Water specific gravity
    #diameter       Inner diameter of pipe, in.
    #angle          angle of pipe inclination in degrees
    #roughness      Roughness of the pipe or tubing
    #Psep           Seperator pressure, psia
    #Tsep           Seperator temperature, °F    
    
    angle_pi = angle * math.pi / 180                       #angle_pi is converted angle in terms of pi
    area = math.pi / 4 * (diameter / 12)**2  #X-sectional area of pipe, ft^2
    gas_rate=liquid_rate*GOR                            #Gas rate, SCF/D 
    liquid_rate1=liquid_rate* 0.000065                  #Unit conversion to ft^3/s from stb/D
    
    #Velocity under zero slippage
    oil_rate=liquid_rate*(1-WC)
    water_rate=liquid_rate*WC
    wor=water_rate/oil_rate
    # usl = liquid_rate1 / area                            #Liquid superficial velocity, ft/s
    # usg = gas_rate / area /86400                         #Gas superficial velocity, ft/s
    
    
    #Importing the properties of fluids
    Z = pvt.zfactor(P,T,gas_grav)                                           #Gas compressibility factor
    TDS = pvt.salinity(wtr_grav)                                            #Water salinity, wt% total dissolved solids
    Pb = pvt.bubble_point2(T, Tsep, Psep, gas_grav, oil_grav, GOR)          #Bubble point pressure, psia
    Rso = pvt.sol_gor(T, P, Tsep, Psep, Pb, gas_grav, oil_grav)             #Solution gas-oil ratio, scf/stb
    Rsw = pvt.sol_gwr(P, T, TDS)                                            #Solution gas_water ratio, scf/stb
    Bo = pvt.oil_fvf(T, P, Tsep, Psep, Pb, Rso, gas_grav, oil_grav)                           #Oil formation volume factor, rb/stb
    Bw = pvt.water_fvf(P, T, TDS)                                           #Water formation volume factor, rb/stb
    Bg = pvt.gas_fvf(P, T, gas_grav)                                        #Gas formation volume factor, ft_/scf
    muo = pvt.oil_visc(T, P, Pb, Rso, gas_grav, oil_grav)                   #Oil viscosity, cp
    muw = pvt.wtr_visc(P, T, TDS)                                           #Water viscosity, cp
    mug = pvt.gvisc(P, T, Z, gas_grav)                                        #Gas viscosity, cp
    rhoo = pvt.oil_dens(P, T, Tsep, Psep, Pb, Bo, Rso, gas_grav, oil_grav)  #Oil density, lb/ft
    rhow = 62.368 * wtr_grav / Bw                                      #Water density, lb/ft^3
    rhog = 2.699 * gas_grav * P / (T + 460) / Z                         #Gas density, lb/ft^3
    sigo = pvt.oil_tens(P, T, oil_grav)                                     #Gas-oil interfacial tension, dynes/cm
    sigw = pvt.wtr_tens(P, T)                                               #Gas-water interfacial tension, dynes/cm
    
    #Volume fraction weighted liquid properties
    rhol = (Bw * wor * rhow + Bo * rhoo) / (Bw * wor + Bo)  
    #rhol     Liquid density, lb/ft                                                  
    mul = (Bw * wor * rhow) / (Bw * wor * rhow + Bo * rhoo) * muw + (Bo * rhoo) / (Bw * wor * rhow + Bo * rhoo) * muo              
    #mul      Liquid viscosity, cp
    sigl = (Bw * wor * rhow) / (Bw * wor * rhow + Bo * rhoo) * sigw + (Bo * rhoo) / (Bw * wor * rhow + Bo * rhoo) * sigo           
    #sigl     Gas-liquid interfacial tension, dynes/cm
    #Calculate bottomhole fluid velocity in ft_/s
    qo = Bo * oil_rate / 15387                                          #Oil flowrate
    qw = Bw * wor * oil_rate / 15387                                    #Water flowrate
    ql = qo + qw                                                        #Liquid flowrate
    if ((GOR - Rso - Rsw*wor) <= 0):                           #If gas flowrate is negative, set to zero
        qg = 0
    else:
        qg = Bg * (GOR - Rso - Rsw * wor) * oil_rate / 86400
       
    usl = ql / area                                                      #Liquid superficial velocity
    usg = qg / area
    um = usl + usg                                       #Mixture superficial velocity, ft/s


    #um       mixture velocity
    #diameter pipe inside diameter
    NFr = um**2/32.174/(diameter/12)  
    CL = ql/ (ql+ qg)
    CG = 1- CL
    NLV = 1.938 * usl * (rhol / sigl) ** 0.25 
    """
      usl is no slip liquid velocity, 
      ρL, liquid density, 
      g, gravitational constant and 
      σ,  is surface tension
    """
    #The transition lines for correlation    
    L1 = 316*CL**0.302
    L2 = 0.0009252*CL**(-2.4684)
    L3 = 0.1*CL**(-1.4516)
    L4 = 0.5*CL**(-6.738) 
                             
    """"
    EL0 = a*CL**b / Frm**c
    """
    #For uphill&downhill flow flow
    #β    the inclination correction factor coefficient
    #B    the liquid hold-up inclination correction factor
    #EL0  the horizontal hold-up 
    #CL   the non slip holdup
    #EL   final hold-up after inclination and slip correlation
    
    
    

    def Flow_type(NFr, CL, L1, L2, L3, L4):
        """Function to Determine the Flow regime by the Method of Beggs and Brill"""
        #The function returns a number indicating the flow flow_type
        #   Segregated flow
        #   Transition flow
        #   Intermittent flow
        #   Distributed flow    
        #NFr        Froude Number
        #CL       Input liquid fraction
        #L1,2,3,4   Dimensionless constants
    
        #flow_type 1 - Segregated flow
        if (((CL < 0.01) and (NFr < L1)) or ((CL >= 0.01) and (NFr < L2))):
            flow_type = "Segregated Flow"
        
            
        #flow_type 2 - Transition flow
        if ((CL >= 0.01) and (L2 < NFr) and (NFr <= L3)):
            flow_type = "Transition Flow"
        
            
        #flow_type 3 - Intermittent flow
        if ((((0.01 <= CL) and (CL < 0.4)) and ((L3 < NFr) and (NFr < L1))) or ((CL >= 0.4) and (L3 < NFr) and (NFr <= L4))):
            flow_type = "Intermittent Flow"
        
            
        #flow_type 4 - Distributed flow
        if (((CL < 0.4) and (NFr >= L1)) or ((CL >= 0.4) and (NFr > L4))):
            flow_type = "Distributed Flow"
        
        
        return flow_type


    def holdup(flow_type, NFr, NLV, CL, angle):
        angle_pi = angle * math.pi / 180      
        if flow_type=="Segregated Flow":
            a,b,c=0.98,0.4846,0.0868
            if (angle >= 0):        #Uphill Flow
                d,e,f,g=0.011,-3.768,3.539,-1.614
            else:                   #Downhill Flow
                d,e,f,g=4.7,-0.3692,-0.3692,-0.5056
                
        elif flow_type=="Intermittent Flow":
            a,b,c=0.845,0.5351,0.0173
            if (angle >= 0):        #Uphill Flow
                d,e,f,g=2.96,0.305,-0.4473,0.0978
            else:                   #Downhill Flow
                d,e,f,g=4.7,-0.3692,-0.3692,-0.5056
        elif flow_type=="Distributed Flow":
            a,b,c=1.065,0.5824,0.0609
            if (angle >= 0):        #Uphill Flow
                d,e,f,g=1,0,0,0
            else:                   #Downhill Flow
                d,e,f,g=4.7,-0.3692,-0.3692,-0.5056
    
        β = (1 - CL)*math.log( d*CL**e*(NLV**f*(NFr**g )))
        if β<0:
            β=0
            
        B = 1 + β * (math.sin(1.8 * angle_pi) - (math.sin(1.8 * angle_pi)) ** 3 / 3)
        EL0 = a*CL**b / NFr**c
        EL= EL0* B
        return EL
    
    flow_type= Flow_type(NFr, CL, L1, L2, L3, L4)
    
    if flow_type=="Transition Flow":
        EL_segregated = holdup("Segregated Flow", NFr, NLV, CL, angle)
        EL_intermittent = holdup("Intermittent Flow", NFr, NLV, CL, angle)
        A = ( L3 - NFr)/(L3 - L2)
        B = 1- A
        
        EL = A*EL_segregated + B*EL_intermittent
        
        
    else:
        EL = holdup(flow_type, NFr, NLV, CL, angle)
        

    y = CL / EL**2
    if y>1 and y<1.2:
        S = math.log(2.2*y - 1.2)
    else:
        S = math.log(y)/(-0.0523 + 3.182*math.log(y) - 0.8725*(math.log(y))**2 + 0.01853*(math.log(y))**4 )


    
    #Calculate non slip fluid mixture properties
    rhom = rhol * CL+ rhog * CG      #Input fraction weighted density, lb/ft_
    mum = mul * CL + mug * CG         #Input fraction weighted viscosity, cp
    EG = 1-EL
    rhobar = rhol * EL + rhog * EG

    Nre = 1488 * (rhom) * um * (diameter/12) / mum    #Non slip Reynolds number
    
    if Nre<2000:
            fNS=64/Nre
    else:
            #Friction factor estimation
            f0=random.random()  #Randomly selecting a friction factor value between 0 to 1
            constant=(roughness/(3.7*diameter))+(2.51/(Nre*math.sqrt(f0)))
            f=(1/(-2*math.log(constant)))   #Calculating the first f value with Re and rougness
            
            def error(a1,a2):   #Function that finds the error between 2 numbers
                return abs(a1-a2)/a1*100
            
            liste=[]   
            liste.append(round(f0,5))  #Adding the first random f value
            liste.append(round(f,5))    
            for i in range(10):
                
                if error(f0,f)<0.05:
                    break
                else:
                    f0=f
                    constant=(roughness/(3.7*diameter))+(2.51/(Nre*math.sqrt(f0)))
                    f=(1/(-2*math.log(constant,10)))**2
                    f=round(f,5)
                    liste.append(f)
            
            fNS=f                          #fNS is then calculated using Colebrook-White equation
        
    fTP = (math.e**S)*fNS                       #Two phase friction factor

    elevation_loss = rhobar * math.sin(angle_pi)/144  #Potential energy pressure gradient, psi/ft
    friction_loss = 2*fTP*um**2*rhom /32.17/(diameter/12)/144  #Frictional pressure gradient, psi/ft
    Ek = um * usg * rhobar / 32.17 / P / 144
    total_pres_loss_grad= (friction_loss +elevation_loss)/(1-Ek)     #Overall pressure gradient, psi/ft   

    
    return total_pres_loss_grad

    
P,T,liquid_rate, WC, GOR, gas_grav, oil_grav, wtr_grav, diameter, angle, roughness, Psep,Tsep=3000,150,100,0.1,100,0.75,35,1.121,3,90,0.005, 114.7,50
beggsandbrill(P,T,liquid_rate, WC, GOR, gas_grav, oil_grav, wtr_grav, diameter, angle, roughness, Psep,Tsep)