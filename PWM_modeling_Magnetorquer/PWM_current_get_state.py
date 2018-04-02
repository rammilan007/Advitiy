def R_B_O_q(q):  #q should be 4x1 vecotr with first scalar element
    a=q[0];b=q[1];c=q[2];d=q[3];
    m11=a*a+b*b-c*c-d*d; m12=2.0*(b*c-a*d); m13=2.0*(b*d+a*c);
    m21=2.0*(b*c+a*d); m22=a*a-b*b+c*c-d*d; m23=2.0*(c*d-a*b);
    m31=2.0*(b*d-a*c); m32=2.0*(c*d+a*b); m33=a*a-b*b-c*c+d*d;
    rot_mat=matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])
    return rot_mat
#---------------------------------------------------------------------------------------------------------------------
#code using 7 Dimensional vector
#if something is wrong please check for 2.0(float) to be written instead of 2(int)
#Aim of this code is to generated angular rate and quaternions data for a orbit
import numpy as np                                     
import matplotlib.pyplot as plt
from numpy import matrix
from random import *
import csv
import math
import matplotlib.pyplot as plt
#-------------Array of Current from PWM_current.py____Written by Sanket's code----------
from PWM_current import current_LR_PWM
#Inputs to above function: inductance, resistance, final time, frequency, duty cycle(from 0.0 to 1.0) [all floats]

tf = 10.0; freq = 100.0; duty = 0.5
dt_h = 0.01*duty/freq 	#step size for high cycle
dt_l = 0.01*(1.0-duty)/freq 	#step size for low cycle
time_i_v = current_LR_PWM(1.0,100.0,tf,freq,duty,dt_h,dt_l)
I=time_i_v[1][0:1000];  #array of current of size 1000
time_arr= time_i_v[0][0:1000]; #array of time as per dt_h and dt_l(variable) step sizes
#-------------------------------------------------------------------------------------------

time_interval=50 # time in seconds for which code should run
h=0.1    #time step =h =0.1 second
time_arr =np.linspace(0,time_interval,1+time_interval/h) 

        

Area= 0.01 # 10cmX10cm =0.01 sq. metre-->Area of coil

mag_moment=I*Area             #  array of magnetic moment across one coil[In SI unit]

Inertia_matrix=matrix( [[1.,2,0.5],[1,2.,5.0],[1.2,2.9,3.5]])   # Creates a 3x3 matrix
def x_dot(x,Torque):    # x should be a matrix of dimension 7X1 with single column first 4 quaternions and next 3 omega
                        #first element of quaternion is scalar
    q=x[0:4]; w=x[4:7];
    omega_dot=(Inertia_matrix.I)*(Torque - matrix(np.cross(w.T, (Inertia_matrix*w).T)).T)
    W=matrix([[0,w[2],-w[1],w[0]],[-w[2],0,w[0],w[1]],[w[1],-w[0],0,w[2]],[-w[0],-w[1],w[2],0]])  #notice 'W' and 'w' here 
    q_dot=0.5*W*q;    #quaternion kinematics differential equation 
    return np.vstack((q_dot,omega_dot))   #dimension 7x1

x_initial=matrix( [[0.84313],[0.30189],[0.44275],[0.0443],[5.0],[6.0],[7.0]] )   #initial omega 7X1
x=x_initial
x_store=x_initial.T

#read data from Bi_120K.csv --> magnetic field in ECIF
b1=[];b2=[];b3=[];mag_file_address='C:/Users/150010037/Desktop/Advitiy/PWM_modeling_Torquer/Bi_120k.csv'
with open(mag_file_address) as csvfile:
    data = csv.reader(csvfile)
    for j in range(int(1+time_interval/h)):     # for running upto the time i want because SGP file may contain data of larger time
        row=data.next()
        b1.append(float(row[1]));b2.append(float(row[2]));b3.append(float(row[3]));
#b1,b2 and b3 are list of components of Magnetic field in ECIF at successive time

#read data from SGP and construct orbit frame
sgp_file_address='C:/Users/150010037/Desktop/Advitiy/PWM_modeling_Torquer/SGP_120k.csv'
with open(sgp_file_address) as csvfile:
    data = csv.reader(csvfile)
    for i in range(int(1+time_interval/h)):     # for running upto the time i want because SGP file may contain data of larger time
        row=data.next()
        
        r1=float(row[1]); r2=float(row[2]); r3=float(row[3]);  # components of r vector in ECIF
        z_O_unit=np.array([r1,r2,r3])/math.sqrt(r1*r1+r2*r2+r3*r3)
        v1=float(row[4]); v2=float(row[5]); v3=float(row[6]);  # components of v vector in ECIF
        r_cross_v=np.cross([v1,v2,v3],[r1,r2,r3])
        y_O_unit=r_cross_v/math.sqrt(r_cross_v[0]**2 + r_cross_v[1]**2 + r_cross_v[2]**2)
        x_O_unit = np.cross(y_O_unit,z_O_unit)
        
        R_O_I = matrix( [x_O_unit,y_O_unit,-z_O_unit])  # rotation matrix for ECIF to Body 
        q=x[0:4] ;q=np.squeeze(np.asarray(q.T))  #to squeeze matrix to array since R_B_O_q function takes array
        R_B_O=R_B_O_q(q)   #look one cell up for definition
        R_B_I=R_B_O*R_O_I   #rotation matrix for Inertial to body 
        B_vector_I= matrix([[b1[i]],[b2[i]],[b3[i]]])
        B_vector_B=R_B_I*B_vector_I
        
        mag_moment_vector=matrix([[mag_moment[i]],[mag_moment[i]],[mag_moment[i]]])
        Torque=np.cross(mag_moment_vector.T,B_vector_B.T).T
        #Torque=np.cross(np.squeeze(np.asarray(mag_moment_vector.T)),np.squeeze(np.asarray(B.T)))
        #Torque=matrix(Torque).T
        #RK4 Solver for x
        a=x_dot(x,Torque)         
        b=x_dot(x + a*h/2.0  ,Torque)  
        c=x_dot(x + b*h/2.0 ,Torque)
        d=x_dot(x + c*h ,Torque)
        x=x + (h/6.0)*(a + 2*b + 2*c + d)     # dimension 7x1
        x_store=np.concatenate((x_store,x.T), axis=0); #nx7 matrix-> first 4 column quaternion then next 3 is omega
    print x_store


    
