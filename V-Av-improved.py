import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy.optimize import minimize

def jack(x,j):
    r=0
    for i in range(len(x)):
        if i!=j:
            r=r+x[i]
    return 1/(len(x)-1)*r        

def exp_val(data):
    return sum(data)/len(data)

def var(data):
    res=0
    for i in range(len(data)):
        res=res+(jack(data,i)-exp_val(data))**2
    return np.sqrt((len(data)-1)/len(data)*res)



#########decide here which nsq
nsq=5
##########


f = h5py.File("BsDsStar.h5", "r")

mom=[[''],
     ['final_state_GX/operator_GammaY/n2_1/0_0_1','final_state_GY/operator_GammaZ/n2_1/1_0_0', 'final_state_GZ/operator_GammaX/n2_1/0_1_0','final_state_GY/operator_GammaX/n2_1/0_0_1','final_state_GZ/operator_GammaY/n2_1/1_0_0', 'final_state_GX/operator_GammaZ/n2_1/0_1_0'],
     ['final_state_GY/operator_GammaZ/n2_2/1_1_0','final_state_GX/operator_GammaY/n2_2/0_1_1','final_state_GZ/operator_GammaX/n2_2/1_1_0','final_state_GX/operator_GammaZ/n2_2/1_1_0','final_state_GZ/operator_GammaY/n2_2/1_1_0','final_state_GY/operator_GammaX/n2_2/1_0_1','final_state_GY/operator_GammaZ/n2_2/1_0_1','final_state_GX/operator_GammaY/n2_2/1_0_1','final_state_GZ/operator_GammaX/n2_2/0_1_1','final_state_GX/operator_GammaZ/n2_2/0_1_1','final_state_GZ/operator_GammaY/n2_2/1_0_1','final_state_GY/operator_GammaX/n2_2/0_1_1'],
     ['final_state_GX/operator_GammaY/n2_3/1_1_1','final_state_GY/operator_GammaZ/n2_3/1_1_1','final_state_GZ/operator_GammaX/n2_3/1_1_1'],
     ['final_state_GX/operator_GammaY/n2_4/0_0_2','final_state_GY/operator_GammaZ/n2_4/2_0_0', 'final_state_GZ/operator_GammaX/n2_4/0_2_0','final_state_GY/operator_GammaX/n2_4/0_0_2','final_state_GZ/operator_GammaY/n2_4/2_0_0', 'final_state_GX/operator_GammaZ/n2_4/0_2_0'],
     ['final_state_GX/operator_GammaY/n2_5/0_2_1','final_state_GY/operator_GammaZ/n2_5/1_2_0','final_state_GZ/operator_GammaX/n2_5/2_1_0','final_state_GX/operator_GammaZ/n2_5/2_1_0','final_state_GY/operator_GammaX/n2_5/2_0_1','final_state_GZ/operator_GammaY/n2_5/1_2_0','final_state_GX/operator_GammaY/n2_5/2_0_1','final_state_GY/operator_GammaZ/n2_5/1_0_2','final_state_GZ/operator_GammaX/n2_5/0_1_2','final_state_GX/operator_GammaZ/n2_5/0_1_2','final_state_GY/operator_GammaX/n2_5/0_2_1','final_state_GZ/operator_GammaY/n2_5/1_0_2',
     'final_state_GX/operator_GammaY/n2_5/0_1_2','final_state_GY/operator_GammaZ/n2_5/2_1_0','final_state_GZ/operator_GammaX/n2_5/1_2_0','final_state_GX/operator_GammaZ/n2_5/1_2_0','final_state_GY/operator_GammaX/n2_5/1_0_2','final_state_GZ/operator_GammaY/n2_5/2_1_0','final_state_GX/operator_GammaY/n2_5/1_0_2','final_state_GY/operator_GammaZ/n2_5/2_0_1','final_state_GZ/operator_GammaX/n2_5/0_2_1','final_state_GX/operator_GammaZ/n2_5/0_2_1','final_state_GY/operator_GammaX/n2_5/0_1_2','final_state_GZ/operator_GammaY/n2_5/2_0_1']]
     

ptmom=['0_0_0','1_0_0','1_1_0','1_1_1', '2_0_0','2_1_0']


if len(mom[nsq])==3:
    dset1xn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][0])]
    dset1xn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][0])]
    dset1yn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][1])]
    dset1yn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][1])]
    dset1zn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][2])]
    dset1zn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][2])]
    nmom=3
elif len(mom[nsq])==6:
    dset1xn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][0])]
    dset1xn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][0])]
    dset1yn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][1])]
    dset1yn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][1])]
    dset1zn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][2])]
    dset1zn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][2])]   
    dset1xn0m=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][3])]
    dset1xn0bm=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][3])]
    dset1yn0m=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][4])]
    dset1yn0bm=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][4])]
    dset1zn0m=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][5])]
    dset1zn0bm=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][5])]   
    nmom=6
elif len(mom[nsq])==12:
    dset1xn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][0])]
    dset1xn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][0])]
    dset1yn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][1])]
    dset1yn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][1])]
    dset1zn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][2])]
    dset1zn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][2])]   
    dset1xn0m=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][3])]
    dset1xn0bm=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][3])]
    dset1yn0m=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][4])]
    dset1yn0bm=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][4])]
    dset1zn0m=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][5])]
    dset1zn0bm=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][5])]   
    
    dset1xn0n=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][6])]
    dset1xn0bn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][6])]
    dset1yn0n=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][7])]
    dset1yn0bn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][7])]
    dset1zn0n=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][8])]
    dset1zn0bn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][8])]   
    dset1xn0mn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][9])]
    dset1xn0bmn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][9])]
    dset1yn0mn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][10])]
    dset1yn0bmn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][10])]
    dset1zn0mn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][11])]
    dset1zn0bmn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][11])]   
    nmom=12  
    
elif len(mom[nsq])==24:
    dset1xn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][0])]
    dset1xn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][0])]
    dset1yn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][1])]
    dset1yn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][1])]
    dset1zn0=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][2])]
    dset1zn0b=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][2])]   
    dset1xn0m=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][3])]
    dset1xn0bm=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][3])]
    dset1yn0m=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][4])]
    dset1yn0bm=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][4])]
    dset1zn0m=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][5])]
    dset1zn0bm=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][5])]   
    
    dset1xn0n=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][6])]
    dset1xn0bn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][6])]
    dset1yn0n=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][7])]
    dset1yn0bn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][7])]
    dset1zn0n=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][8])]
    dset1zn0bn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][8])]   
    dset1xn0mn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][9])]
    dset1xn0bmn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][9])]
    dset1yn0mn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][10])]
    dset1yn0bmn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][10])]
    dset1zn0mn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][11])]
    dset1zn0bmn=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][11])]   
    
    dset1xn0p=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][12])]
    dset1xn0bp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][12])]
    dset1yn0p=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][13])]
    dset1yn0bp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][13])]
    dset1zn0p=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][14])]
    dset1zn0bp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][14])]   
    dset1xn0mp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][15])]
    dset1xn0bmp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][15])]
    dset1yn0mp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][16])]
    dset1yn0bmp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][16])]
    dset1zn0mp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][17])]
    dset1zn0bmp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][17])]   
    
    dset1xn0np=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][18])]
    dset1xn0bnp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][18])]
    dset1yn0np=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][19])]
    dset1yn0bnp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][19])]
    dset1zn0np=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][20])]
    dset1zn0bnp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][20])]   
    dset1xn0mnp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][21])]
    dset1xn0bmnp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][21])]
    dset1yn0mnp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][22])]
    dset1yn0bmnp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][22])]
    dset1zn0mnp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][23])]
    dset1zn0bmnp=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][23])]   
    nmom=24      


dsxn0=f["/CHARM_SM12.14_SM12.14_s0.02144/c0.248/operator_GammaX/{}/data".format(ptmom[nsq])]
dsyn0=f["/CHARM_SM12.14_SM12.14_s0.02144/c0.248/operator_GammaY/{}/data".format(ptmom[nsq])]
dszn0=f["/CHARM_SM12.14_SM12.14_s0.02144/c0.248/operator_GammaZ/{}/data".format(ptmom[nsq])]
bsn0=f["/rhq_m2.42_csw2.68_zeta1.52_SM12.14_SM12.14_s0.02144/operator_Gamma5/0_0_0/data"]



#####Eff Masses
mdlist=[0.7348303222656272,0.7458868408203149,0.7567413330078149,0.7673809814453149,0.77745605,0.787656860351565]
md=mdlist[nsq]

# Constants

mb = 1.9257122802734448
md0 = 0.73483032
pre=-(mb+md0)/(2*mb)

if nsq ==4: pre=pre*1/2

nconf=98
dt=30
ts=96

# Initialize arrays
av1n0x = np.zeros((dt, nconf))
av1n0y = np.zeros((dt, nconf))
av1n0z = np.zeros((dt, nconf))
av1n0xm = np.zeros((dt, nconf))
av1n0ym = np.zeros((dt, nconf))
av1n0zm = np.zeros((dt, nconf))

av1n0xn = np.zeros((dt, nconf))
av1n0yn = np.zeros((dt, nconf))
av1n0zn = np.zeros((dt, nconf))
av1n0xmn = np.zeros((dt, nconf))
av1n0ymn = np.zeros((dt, nconf))
av1n0zmn = np.zeros((dt, nconf))

av1n0xp = np.zeros((dt, nconf))
av1n0yp = np.zeros((dt, nconf))
av1n0zp = np.zeros((dt, nconf))
av1n0xmp = np.zeros((dt, nconf))
av1n0ymp = np.zeros((dt, nconf))
av1n0zmp = np.zeros((dt, nconf))

av1n0xnp = np.zeros((dt, nconf))
av1n0ynp = np.zeros((dt, nconf))
av1n0znp = np.zeros((dt, nconf))
av1n0xmnp = np.zeros((dt, nconf))
av1n0ymnp = np.zeros((dt, nconf))
av1n0zmnp = np.zeros((dt, nconf))

avdx = np.zeros((dt, nconf))
avdy = np.zeros((dt, nconf))
avdz = np.zeros((dt, nconf))
avb = np.zeros((dt, nconf))
avn0 = np.zeros(dt)


# Loop optsmizatsons
for j in range(dt):
    temp1x = 0
    temp1y = 0
    temp1z = 0
    temp1xm = 0
    temp1ym = 0
    temp1zm = 0
    
    temp1xn = 0
    temp1yn = 0
    temp1zn = 0
    temp1xmn = 0
    temp1ymn = 0
    temp1zmn = 0
    
    temp1xp = 0
    temp1yp = 0
    temp1zp = 0
    temp1xmp = 0
    temp1ymp = 0
    temp1zmp = 0
    
    temp1xnp = 0
    temp1ynp = 0
    temp1znp = 0
    temp1xmnp = 0
    temp1ymnp = 0
    temp1zmnp = 0
    
    tempdx = 0
    tempdy = 0
    tempdz = 0
    tempb = 0
    
    for k in range(nconf):
        tmp1x = np.mean((np.imag(dset1xn0[k, :, j]) - np.imag(dset1xn0b[k, :, dt-j]))) / 2
        tmp1y = np.mean((np.imag(dset1yn0[k, :, j]) - np.imag(dset1yn0b[k, :, dt-j]))) / 2
        tmp1z = np.mean((np.imag(dset1zn0[k, :, j]) - np.imag(dset1zn0b[k, :, dt-j]))) / 2
        if nmom==6:
            tmp1xm = np.mean((np.imag(dset1xn0m[k, :, j]) - np.imag(dset1xn0bm[k, :, dt-j]))) / 2
            tmp1ym = np.mean((np.imag(dset1yn0m[k, :, j]) - np.imag(dset1yn0bm[k, :, dt-j]))) / 2
            tmp1zm = np.mean((np.imag(dset1zn0m[k, :, j]) - np.imag(dset1zn0bm[k, :, dt-j]))) / 2
            tmp1xn=0
            tmp1yn=0
            tmp1zn=0
            
            tmp1xmn=0
            tmp1ymn=0
            tmp1zmn=0
            
            tmp1xp=0
            tmp1yp=0
            tmp1zp=0
            
            tmp1xmp=0
            tmp1ymp=0
            tmp1zmp=0
            
            tmp1xnp=0
            tmp1ynp=0
            tmp1znp=0
            
            tmp1xmnp=0
            tmp1ymnp=0
            tmp1zmnp=0
            
        elif nmom==12:
            tmp1xm = np.mean((np.imag(dset1xn0m[k, :, j]) - np.imag(dset1xn0bm[k, :, dt-j]))) / 2
            tmp1ym = np.mean((np.imag(dset1yn0m[k, :, j]) - np.imag(dset1yn0bm[k, :, dt-j]))) / 2
            tmp1zm = np.mean((np.imag(dset1zn0m[k, :, j]) - np.imag(dset1zn0bm[k, :, dt-j]))) / 2
            
            tmp1xn = np.mean((np.imag(dset1xn0n[k, :, j]) - np.imag(dset1xn0bn[k, :, dt-j]))) / 2
            tmp1yn = np.mean((np.imag(dset1yn0n[k, :, j]) - np.imag(dset1yn0bn[k, :, dt-j]))) / 2
            tmp1zn = np.mean((np.imag(dset1zn0n[k, :, j]) - np.imag(dset1zn0bn[k, :, dt-j]))) / 2
            tmp1xmn = np.mean((np.imag(dset1xn0mn[k, :, j]) - np.imag(dset1xn0bmn[k, :, dt-j]))) / 2
            tmp1ymn = np.mean((np.imag(dset1yn0mn[k, :, j]) - np.imag(dset1yn0bmn[k, :, dt-j]))) / 2
            tmp1zmn = np.mean((np.imag(dset1zn0mn[k, :, j]) - np.imag(dset1zn0bmn[k, :, dt-j]))) / 2
            
            tmp1xp=0
            tmp1yp=0
            tmp1zp=0
            
            tmp1xmp=0
            tmp1ymp=0
            tmp1zmp=0
            
            tmp1xnp=0
            tmp1ynp=0
            tmp1znp=0
            
            tmp1xmnp=0
            tmp1ymnp=0
            tmp1zmnp=0
        
        elif nmom==24:
             tmp1xm = np.mean((np.imag(dset1xn0m[k, :, j]) - np.imag(dset1xn0bm[k, :, dt-j]))) / 2
             tmp1ym = np.mean((np.imag(dset1yn0m[k, :, j]) - np.imag(dset1yn0bm[k, :, dt-j]))) / 2
             tmp1zm = np.mean((np.imag(dset1zn0m[k, :, j]) - np.imag(dset1zn0bm[k, :, dt-j]))) / 2
             
             tmp1xn = np.mean((np.imag(dset1xn0n[k, :, j]) - np.imag(dset1xn0bn[k, :, dt-j]))) / 2
             tmp1yn = np.mean((np.imag(dset1yn0n[k, :, j]) - np.imag(dset1yn0bn[k, :, dt-j]))) / 2
             tmp1zn = np.mean((np.imag(dset1zn0n[k, :, j]) - np.imag(dset1zn0bn[k, :, dt-j]))) / 2
             tmp1xmn = np.mean((np.imag(dset1xn0mn[k, :, j]) - np.imag(dset1xn0bmn[k, :, dt-j]))) / 2
             tmp1ymn = np.mean((np.imag(dset1yn0mn[k, :, j]) - np.imag(dset1yn0bmn[k, :, dt-j]))) / 2
             tmp1zmn = np.mean((np.imag(dset1zn0mn[k, :, j]) - np.imag(dset1zn0bmn[k, :, dt-j]))) / 2
             
             tmp1xp = np.mean((np.imag(dset1xn0p[k, :, j]) - np.imag(dset1xn0bp[k, :, dt-j]))) / 2
             tmp1yp = np.mean((np.imag(dset1yn0p[k, :, j]) - np.imag(dset1yn0bp[k, :, dt-j]))) / 2
             tmp1zp = np.mean((np.imag(dset1zn0p[k, :, j]) - np.imag(dset1zn0bp[k, :, dt-j]))) / 2
             tmp1xmp = np.mean((np.imag(dset1xn0mp[k, :, j]) - np.imag(dset1xn0bmp[k, :, dt-j]))) / 2
             tmp1ymp = np.mean((np.imag(dset1yn0mp[k, :, j]) - np.imag(dset1yn0bmp[k, :, dt-j]))) / 2
             tmp1zmp = np.mean((np.imag(dset1zn0mp[k, :, j]) - np.imag(dset1zn0bmp[k, :, dt-j]))) / 2
             
             tmp1xnp = np.mean((np.imag(dset1xn0np[k, :, j]) - np.imag(dset1xn0bnp[k, :, dt-j]))) / 2
             tmp1ynp = np.mean((np.imag(dset1yn0np[k, :, j]) - np.imag(dset1yn0bnp[k, :, dt-j]))) / 2
             tmp1znp = np.mean((np.imag(dset1zn0np[k, :, j]) - np.imag(dset1zn0bnp[k, :, dt-j]))) / 2
             tmp1xmnp = np.mean((np.imag(dset1xn0mnp[k, :, j]) - np.imag(dset1xn0bmnp[k, :, dt-j]))) / 2
             tmp1ymnp = np.mean((np.imag(dset1yn0mnp[k, :, j]) - np.imag(dset1yn0bmnp[k, :, dt-j]))) / 2
             tmp1zmnp = np.mean((np.imag(dset1zn0mnp[k, :, j]) - np.imag(dset1zn0bmnp[k, :, dt-j]))) / 2
             
        
        else:
            tmp1xm=0
            tmp1ym=0
            tmp1zm=0
            
            tmp1xn=0
            tmp1yn=0
            tmp1zn=0
            
            tmp1xmn=0
            tmp1ymn=0
            tmp1zmn=0
            
            
        tmpdx = np.mean(np.real(dsxn0[k, :, j+1]) + np.real(dsxn0[k, :, ts-1-j])) / 2 if j != 0 else np.mean(np.real(dsxn0[k, :, 0]))
        tmpdy = np.mean(np.real(dsyn0[k, :, j+1]) + np.real(dsyn0[k, :, ts-1-j])) / 2 if j != 0 else np.mean(np.real(dsyn0[k, :, 0]))
        tmpdz = np.mean(np.real(dszn0[k, :, j+1]) + np.real(dszn0[k, :, ts-1-j])) / 2 if j != 0 else np.mean(np.real(dszn0[k, :, 0]))
        tmpb = np.mean(np.real(bsn0[k, :, dt-1-j]) + np.real(bsn0[k, :, ts-dt+1+j])) / 2 if j != dt-1 else np.mean(np.real(bsn0[k, :, 0]))


        av1n0x[j, k] = tmp1x
        av1n0y[j, k] = tmp1y
        av1n0z[j, k] = tmp1z
        av1n0xm[j, k] = tmp1xm
        av1n0ym[j, k] = tmp1ym
        av1n0zm[j, k] = tmp1zm
        
        av1n0xn[j, k] = tmp1xn
        av1n0yn[j, k] = tmp1yn
        av1n0zn[j, k] = tmp1zn
        av1n0xmn[j, k] = tmp1xmn
        av1n0ymn[j, k] = tmp1ymn
        av1n0zmn[j, k] = tmp1zmn
        
        av1n0xp[j, k] = tmp1xp
        av1n0yp[j, k] = tmp1yp
        av1n0zp[j, k] = tmp1zp
        av1n0xmp[j, k] = tmp1xmp
        av1n0ymp[j, k] = tmp1ymp
        av1n0zmp[j, k] = tmp1zmp
        
        av1n0xnp[j, k] = tmp1xnp
        av1n0ynp[j, k] = tmp1ynp
        av1n0znp[j, k] = tmp1znp
        av1n0xmnp[j, k] = tmp1xmnp
        av1n0ymnp[j, k] = tmp1ymnp
        av1n0zmnp[j, k] = tmp1zmnp
        
        avdx[j, k] = tmpdx
        avdy[j, k] = tmpdy
        avdz[j, k] = tmpdz
        avb[j, k] = tmpb

        temp1x += tmp1x
        temp1y += tmp1y
        temp1z += tmp1z
        temp1xm += tmp1xm
        temp1ym += tmp1ym
        temp1zm += tmp1zm
        
        temp1xn += tmp1xn
        temp1yn += tmp1yn
        temp1zn += tmp1zn
        temp1xmn += tmp1xmn
        temp1ymn += tmp1ymn
        temp1zmn += tmp1zmn
        
        temp1xp += tmp1xp
        temp1yp += tmp1yp
        temp1zp += tmp1zp
        temp1xmp += tmp1xmp
        temp1ymp += tmp1ymp
        temp1zmp += tmp1zmp
        
        temp1xnp += tmp1xnp
        temp1ynp += tmp1ynp
        temp1znp += tmp1znp
        temp1xmnp += tmp1xmnp
        temp1ymnp += tmp1ymnp
        temp1zmnp += tmp1zmnp
        
        tempdx += tmpdx
        tempdy += tmpdy
        tempdz += tmpdz
        tempb += tmpb
    avn0[j] = pre * (((temp1x + temp1y + temp1z-temp1xm - temp1ym - temp1zm +temp1xn +temp1yn + temp1zn-temp1xmn - temp1ymn - temp1zmn+(temp1xp + temp1yp + temp1zp-temp1xmp - temp1ymp - temp1zmp +temp1xnp +temp1ynp + temp1znp-temp1xmnp - temp1ymnp - temp1zmnp)) / nmom) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    #+temp1xn + temp1yn + temp1zn+temp1xmn + temp1ymn + temp1zmn

errn0=np.zeros(shape=(dt))
#errnn0=np.zeros(shape=(96))

for j in range(dt):
    x=0
    for i in range(nconf):
        x=x+((((jack(av1n0x[j],i)+jack(av1n0y[j],i)+jack(av1n0z[j],i)-jack(av1n0xm[j],i)-jack(av1n0ym[j],i)-jack(av1n0zm[j],i)+jack(av1n0xn[j],i)+jack(av1n0yn[j],i)+jack(av1n0zn[j],i)-jack(av1n0xmn[j],i)-jack(av1n0ymn[j],i)-jack(av1n0zmn[j],i)+(jack(av1n0xp[j],i)+jack(av1n0yp[j],i)+jack(av1n0zp[j],i)-jack(av1n0xmp[j],i)-jack(av1n0ymp[j],i)-jack(av1n0zmp[j],i)+jack(av1n0xnp[j],i)+jack(av1n0ynp[j],i)+jack(av1n0znp[j],i)-jack(av1n0xmnp[j],i)-jack(av1n0ymnp[j],i)-jack(av1n0zmnp[j],i)))/nmom)/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(30-j))))*pre-avn0[j])**2
    errn0[j]=np.sqrt((nconf-1)/nconf*x) 
    #errnn0[j]=np.sqrt((98-1)/98*x)*10**(85)

avn0[np.isnan(avn0)] = 0
errn0[np.isnan(errn0)] = 0

#plt.plot(list(range(96)), np.absolute(avn0)[0:96],'ro')
plt.xlabel('time')
plt.ylabel(r'$\widetilde{V}$')
plt.errorbar(list(range(dt)), np.absolute(avn0)[0:dt], yerr=errn0[0:dt],ls='none',fmt='x',label='nsq={}'.format(nsq))
#plt.yscale('log')
plt.legend()
plt.savefig('./V/V-nsq{}.png'.format(nsq))

np.savetxt('./V/V-nsq{}.txt'.format(nsq), np.c_[np.absolute(avn0), errn0])

###############################################################################

reg_low=18
reg_up=25

#Covarianze matrix (without prefactor, not squarrooted)
cut=ts/2-1-reg_up
cut1=ts/2+1-reg_up
covmat=np.zeros(shape=(int(ts/2-1-reg_low-cut),int(ts/2-1-reg_low-cut)))
for t1 in range(int(ts/2-1-reg_low-cut)):
    for t2 in range(int(ts/2-1-reg_low-cut)):
        x=0
        for i in range(nconf):  
            #+jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)+jack(av1n0zm[t1+reg_low],i)
            #+jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)+jack(av1n0zm[t2+reg_low],i)
            if nmom==3:
                x=x+((((jack(av1n0x[t1+reg_low],i)+jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i))/3)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)+jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i))/3)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(30-(t2+reg_low)))))*pre-avn0[t2])            
            elif nmom==6:
                x=x+((((jack(av1n0x[t1+reg_low],i)+jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)-jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)-jack(av1n0zm[t1+reg_low],i))/6)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)+jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)-jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)-jack(av1n0zm[t2+reg_low],i))/6)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
            elif nmom==12:
                x=x+((((jack(av1n0x[t1+reg_low],i)+jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)-jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)-jack(av1n0zm[t1+reg_low],i)+jack(av1n0xn[t1+reg_low],i)+jack(av1n0yn[t1+reg_low],i)+jack(av1n0zn[t1+reg_low],i)-jack(av1n0xmn[t1+reg_low],i)-jack(av1n0ymn[t1+reg_low],i)-jack(av1n0zmn[t1+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)+jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)-jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)-jack(av1n0zm[t2+reg_low],i)+jack(av1n0xn[t2+reg_low],i)+jack(av1n0yn[t2+reg_low],i)+jack(av1n0zn[t2+reg_low],i)-jack(av1n0xmn[t2+reg_low],i)-
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  jack(av1n0ymn[t2+reg_low],i)-jack(av1n0zmn[t2+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
            elif nmom==24:
                x=x+((((jack(av1n0x[t1+reg_low],i)+jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)-jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)-jack(av1n0zm[t1+reg_low],i)+jack(av1n0xn[t1+reg_low],i)+jack(av1n0yn[t1+reg_low],i)+jack(av1n0zn[t1+reg_low],i)-jack(av1n0xmn[t1+reg_low],i)-jack(av1n0ymn[t1+reg_low],i)-jack(av1n0zmn[t1+reg_low],i)+jack(av1n0xp[t1+reg_low],i)+jack(av1n0yp[t1+reg_low],i)+jack(av1n0zp[t1+reg_low],i)-jack(av1n0xmp[t1+reg_low],i)-jack(av1n0ymp[t1+reg_low],i)-jack(av1n0zmp[t1+reg_low],i)+jack(av1n0xnp[t1+reg_low],i)+jack(av1n0ynp[t1+reg_low],i)+jack(av1n0znp[t1+reg_low],i)-jack(av1n0xmnp[t1+reg_low],i)-jack(av1n0ymnp[t1+reg_low],i)-jack(av1n0zmnp[t1+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)+jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)-jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)-jack(av1n0zm[t2+reg_low],i)+jack(av1n0xn[t2+reg_low],i)+jack(av1n0yn[t2+reg_low],i)+jack(av1n0zn[t2+reg_low],i)-jack(av1n0xmn[t2+reg_low],i)-
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  jack(av1n0ymn[t2+reg_low],i)-jack(av1n0zmn[t2+reg_low],i)+jack(av1n0xp[t2+reg_low],i)+jack(av1n0yp[t2+reg_low],i)+jack(av1n0zp[t2+reg_low],i)-jack(av1n0xmp[t2+reg_low],i)-jack(av1n0ymp[t2+reg_low],i)-jack(av1n0zmp[t2+reg_low],i)+jack(av1n0xnp[t2+reg_low],i)+jack(av1n0ynp[t2+reg_low],i)+jack(av1n0znp[t2+reg_low],i)-jack(av1n0xmnp[t2+reg_low],i)-
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  jack(av1n0ymnp[t2+reg_low],i)-jack(av1n0zmnp[t2+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
   
        covmat[t1][t2]=x
        covmat[t2][t1]=x  
        
        
def chi(a):
    return (nconf-1-reg_low-cut)/(nconf-reg_low-cut)*np.dot(np.transpose([i-a for i in avn0[reg_low:reg_up]]),np.matmul(np.linalg.inv(covmat),[i-a for i in avn0[reg_low:reg_up]]))

mbar=minimize(chi,0.1,method='Nelder-Mead', tol=1e-6)


def jackmass(t1,i):
    #+jack(av1n0xm[t1],i)-jack(av1n0ym[t1],i)+jack(av1n0zm[t1],i)
    return (((jack(av1n0x[t1],i)+jack(av1n0y[t1],i)+jack(av1n0z[t1],i)-jack(av1n0xm[t1],i)-jack(av1n0ym[t1],i)-jack(av1n0zm[t1],i)+jack(av1n0xn[t1],i)+jack(av1n0yn[t1],i)+jack(av1n0zn[t1],i)-jack(av1n0xmn[t1],i)-jack(av1n0ymn[t1],i)-jack(av1n0zmn[t1],i)+jack(av1n0xp[t1],i)+jack(av1n0yp[t1],i)+jack(av1n0zp[t1],i)-jack(av1n0xmp[t1],i)-jack(av1n0ymp[t1],i)-jack(av1n0zmp[t1],i)+jack(av1n0xnp[t1],i)+jack(av1n0ynp[t1],i)+jack(av1n0znp[t1],i)-jack(av1n0xmnp[t1],i)-jack(av1n0ymnp[t1],i)-jack(av1n0zmnp[t1],i))/nmom)/(1/3*np.sqrt((jack(avdx[t1],i)+jack(avdy[t1],i)+jack(avdz[t1],i))*jack(avb[t1],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1))*np.exp(-mb*(30-(t1)))))*pre

def chijack(a,k):
    return np.dot(np.transpose([jackmass(i+reg_low,k)-a for i in range(int(ts/2-1-reg_low-cut))]),np.matmul(np.linalg.inv(covmat),[jackmass(i+reg_low,k)-a for i in range(int(ts/2-1-reg_low-cut))]))


#Std Deviatson for all jakcknife blocks
h=0
for i in range(nconf):
    h=h+(minimize(chijack,0.1,args=(i),method='Nelder-Mead', tol=1e-6).x[0]-mbar.x[0])**2
sigma=np.sqrt((nconf-1-reg_low-cut)/(nconf-reg_low-cut)*h)

print(mbar,sigma)

f = plt.figure(figsize=(10,3))
ax = f.add_subplot(121)

ax.set_xlabel('time')
ax.set_ylabel('V')
ax.errorbar(list(range(dt))[1:dt], avn0[1:dt], yerr=errn0[1:dt],fmt='x', label='nsq={}'.format(nsq))
plt.axhline(y = mbar.x[0], color = 'r', linestyle = 'dashed')
plt.fill_between(list(range(dt))[reg_low:reg_up], mbar.x[0]+sigma, mbar.x[0]-sigma, color='r',alpha=0.2)
#ax.set_yscale('log')
plt.savefig('./Fits/V-Av-nsq{}-Fit.png'.format(nsq))

df3 = pd.DataFrame(columns=['EffectiveMass','Error','RegUp','RegLow'])
df3['EffectiveMass']=mbar.x
df3['Error']=sigma  
df3['RegUp']=reg_up
df3['RegLow']=reg_low    
df3.to_csv('./Fits/V-Av-nsq{}-Fit.csv'.format(nsq), sep='\t')
