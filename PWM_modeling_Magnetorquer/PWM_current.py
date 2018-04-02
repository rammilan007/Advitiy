#Writer of the code ---> Sanket
#Function to calculate the current from t=0 to t=tf  
import numpy as np
import math as math
#Inputs: inductance, resistance, final time, frequency, duty cycle(from 0.0 to 1.0) [all floats]
def current_LR_PWM(L,R,tf,freq,duty,dt_h,dt_l):	#Returns the tuple with [t,i] where t is from 0 to tf.
	#dt = 0.001*duty/freq
	V_m = 3.3
	N = int(duty*tf/dt_h)+int((1.-duty)*tf/dt_l) + 2
	print N
	time = np.zeros(N)
	i = np.zeros(N)
	v = np.zeros(N)
	T = 1.0/freq
	tau = L/R
	i0_h = 0.0	#current at start of high cycle
	i0_l = 0.0	#current at the start of low cycle
	
	t = 0.0
	k = 0 
	t_curr = 0.0
	while (t<=tf):

		 
		if math.fmod(t,T) < duty*T :
			t_curr = math.fmod(t,T)
			i[k] = (V_m/R) + (i0_h - V_m/R)*np.exp(-t_curr/tau)
			v[k] = V_m
			i0_l = i[k]
			t = t + dt_h
		else:
			t_curr = math.fmod(t,T) - duty*T
			i[k] = i0_l*np.exp(-t_curr/tau)
			v[k] = 0.0
			i0_h = i[k]
			t = t +dt_l
		
		k = k+1
		time[k] = t
		#print t
	print k
	return time, i, v



#Test code

import matplotlib.pyplot as plt
tf = 10.0
freq = 100.0
duty = 0.5
dt_h = 0.01*duty/freq 	#step size for high cycle
dt_l = 0.01*(1.0-duty)/freq 	#step size for low cycle
x = current_LR_PWM(1.,100.0,tf,freq,duty,dt_h,dt_l)

plt.plot(x[0][0:1000],x[1][0:1000])
plt.show()
plt.plot(x[0][0:1000],x[2][0:1000])
plt.show()

