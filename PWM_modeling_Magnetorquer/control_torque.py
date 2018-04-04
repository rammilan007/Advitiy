#control torque calculation
from frames import*
import numpy as np
#for Area of the coil(A), import constants
#control_i should  come from PWM_current.py file
def control_torque(control_i,A,v_mag_field_I,q,):
        mag_moment=control_i*A;
        v_mag_field_o=ecif2orbit(r,v,v_mag_field_I);
        v_mag_field_b=np.dot(R_B_O_q(q),v_mag_field_o)
        Torque=np.cross(mag_moment,v_mag_field_b);
        return Torque