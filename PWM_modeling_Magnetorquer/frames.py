#frames 
import numpy as np                                     
import numpy.linalg as la
#------------------ECIF to Orbit frame--------------
def ecif2orbit(r,v,vector_ecif):
    z=-r/la.norm(r);
    y= np.cross(v,r)/la.norm(np.cross(v,r));
    x=np.cross(y,z);
    Rot_mat=np.array([[x],[y],[z]]);
    return np.dot(Rot_mat,vector_ecif)
#---------------------------------------------------------------------------
#Rotation matrix for Orbit to body from  # input is quaternion
def R_B_O_q(q):  #q should be 4x1 vecotr with first scalar element
    a=q[0];b=q[1];c=q[2];d=q[3];
    m11=a*a+b*b-c*c-d*d; m12=2.0*(b*c-a*d); m13=2.0*(b*d+a*c);
    m21=2.0*(b*c+a*d); m22=a*a-b*b+c*c-d*d; m23=2.0*(c*d-a*b);
    m31=2.0*(b*d-a*c); m32=2.0*(c*d+a*b); m33=a*a-b*b-c*c+d*d;
    rot_mat=np.array([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])
    return rot_mat   #3x3 matrix