#Function to generate duty cycle of PWM voltage and also internal step sizes
import numpy as np
import math as math
#Inputs: current from control law, Resistance
def duty_generate(i_control,R):
    duty=(i_control*R)/3.3
    return duty
def internal_step_for_torquer(duty,freq):
    dt_h = 0.01*duty/freq 	#step size for high cycle
    dt_l = 0.01*(1.0-duty)/freq 	#step size for low cycle
    dt_l=dt_l.min();
    return np.hstack((dt_h,dt_l))





# Test Code
#i_control=0.003
#R=100.0
#duty= duty_generate(i_control,R)   
#print duty