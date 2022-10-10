

from tkinter.font import nametofont


W_To_S = 50  #50#N/m2
P_To_W = 7.2  #W/N
AR = 10

MotorShaftPower = 175*0.65

Weight = MotorShaftPower/P_To_W
Mass = Weight/10

S_ref = Weight/W_To_S

Span = (AR*S_ref)**0.5
Chord = S_ref/Span

print('Weight:', Weight)
print('Mass:', Mass)
print('S_ref:', S_ref)
print('Span:', Span)
print('Chord:', Chord)