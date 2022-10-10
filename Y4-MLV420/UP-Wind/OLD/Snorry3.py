import numpy as np
import matplotlib.pyplot as plt


### Sorty:
Laps = np.linspace(2,5,100)                   #
g = 9.81
rho = 1.0687                # kg/m3
Pitstop = 45                #sec 40 s minimum
Time = 5*60                 #sec
TurnRadius = 50             #m  Used in Speed approximation
TurnRadiusMin = 20          #m  Used in power requirement
TakeOffDistance = 25        #m  25 m maximum
RequiredClimbRate = 1.4    #m/s  1m/s minimum

### Estimate performance figures
C_L_TO = 1.3
C_L_Max_Flaps = 2.2         #Langing Configuration

C_Dmin = 0.028
C_D_TO = 0.04

Power_engine = 155  ## Estimate of what will be in with prop
eta_engine = 0.65
eta_prop = 0.7
Power_engine_shaft = Power_engine*eta_engine

    #Speed estimates relative to Cruise:
K_Stall = 0.5
K_TakeOff = 0.7
K_Climb = 0.85

AR = 10
e = 1.78*(1-0.045*(AR**0.68)) - 0.64    ### Oswalds
k = 1/(np.pi * AR * e)                  ### Oswalds
mu = 0.05                              ### Rolling resistance

def CalcTotalMassMoved(Laps):
    TotAns = []
    TotSpan = []
        
    for Lap in Laps:
        LapDistance = 2*150 + np.pi*TurnRadius
        FlightTimePerLap = (Time - (Lap-1)*Pitstop)/Lap
        V_min_ave = LapDistance/FlightTimePerLap
        V_min_Cruise = V_min_ave*1.25
        LoadFactor = (((V_min_Cruise**4)/((TurnRadiusMin**2)*(g**2))) + 1)**0.5  ## eq 19-37

        ### Calcs:


        ### Thrust To weight requirement:
        q_Cruise = 0.5*rho*(V_min_Cruise**2)
        q_TakeOff = 0.5*rho*((V_min_Cruise*K_TakeOff/(2**0.5))**2)
        q_Climb = 0.5*rho*((V_min_Cruise*K_Climb)**2)
        q_Stall = 0.5*rho*((V_min_Cruise*K_Stall)**2)

        Max_W_to_S = C_L_Max_Flaps*q_Stall

        xaxismax = 1.5*round(Max_W_to_S)

        WtoS = np.linspace(0.1*round(Max_W_to_S),xaxismax,1000)           ###N/m2


        TtoW_b = q_Cruise*((C_Dmin/WtoS)+(k*((LoadFactor/q_Cruise)**2)*WtoS))
        TtoW_c = q_Cruise*C_Dmin*(1/WtoS)+(k*(1/q_Cruise)* WtoS)
        TtoW_to = (((V_min_Cruise*K_TakeOff)**2)/(2*g*TakeOffDistance))+((C_D_TO*q_TakeOff)/WtoS) + mu*(1-((C_L_TO*q_TakeOff)/WtoS))
        TtoW_climb = (RequiredClimbRate/(K_Climb*V_min_Cruise))+((C_Dmin*q_Climb)/WtoS)+(k*(1/q_Climb)*WtoS)

        PtoW_Bank = (V_min_Cruise/(eta_engine*eta_prop))*TtoW_b  ### W per N
        PtoW_Cruise = (V_min_Cruise/(eta_engine*eta_prop))*TtoW_c  ### W per N
        PtoW_TO = ((V_min_Cruise*K_TakeOff)/(eta_engine*eta_prop))*TtoW_to  ### W per N
        PtoW_Climb = ((V_min_Cruise*K_Climb)/(eta_engine*eta_prop))*TtoW_climb  ### W per N


        ylim = np.max([PtoW_Bank,PtoW_Cruise,PtoW_TO,PtoW_Climb])

        ### finding optimal
        minimum = 1000
        index = 0
        for idx,i in enumerate(WtoS):
            if i <= Max_W_to_S:
                ans = np.max([PtoW_Bank[idx],PtoW_Cruise[idx],PtoW_TO[idx],PtoW_Climb[idx]])
                if ans < minimum:
                    minimum = ans
                    index = idx



        ChoosenPtoW = minimum*1.08     ## Safety factor
        ChoosenWtoS = WtoS[index]

        Weight = Power_engine_shaft/ChoosenPtoW
        Mass = Weight/g
        S_ref = Weight/ChoosenWtoS

        Span = (AR*S_ref)**0.5
        Chord = S_ref/Span

        
        TotAns.append(Mass*Lap)
        TotSpan.append(Span)

    return TotAns, TotSpan

TotalMass,TotalSpan = CalcTotalMassMoved(Laps)
        


plt.figure(1)
plt.title('Total Mass')

color_r = [0.3,0.3,0.3]
color_b = [0,0,1]
plt.plot(Laps,TotalMass)
plt.xlabel('Laps compleated')
plt.ylabel('Mass Moved [kg]')
plt.grid()

plt.figure(2)
plt.title('Total Span')

color_r = [0.3,0.3,0.3]
color_b = [0,0,1]
plt.plot(Laps,TotalSpan)
plt.xlabel('Laps compleated')
plt.ylabel('Span [kg]')
plt.grid()
plt.show()