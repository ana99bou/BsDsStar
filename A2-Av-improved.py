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

def jack_mean(dat, nmom,j,i):
    temp=0
    for k in range(nmom):
        temp += jack(dat[k][j],i)
    return temp/nmom

def sum_with_exceptions(lst):
    total = 0
    for i, num in enumerate(lst):
        if (i + 1) % 3 == 2:
            total -= num  # Subtract every third element starting from the second
        else:
            total += num  # Add other elements
    return total/len(lst)


def sum_with_exceptions_jack(lst,j,i):
    total = 0
    for i, num in enumerate(lst):
        if (i + 1) % 3 == 2:
            total -= jack(num[j],i)  # Subtract every third element starting from the second
        else:
            total += jack(num[j],i)  # Add other elements
    return total/len(lst)


#########decide here which nsq
nsq=3
##########


f = h5py.File("BsDsStar.h5", "r")

mom=[[''],
     ['final_state_GX/operator_GammaXGamma5/n2_1/1_0_0','final_state_GY/operator_GammaYGamma5/n2_1/0_1_0', 'final_state_GZ/operator_GammaZGamma5/n2_1/0_0_1','final_state_GX/operator_GammaXGamma5/n2_1/-1_0_0','final_state_GY/operator_GammaYGamma5/n2_1/0_-1_0', 'final_state_GZ/operator_GammaZGamma5/n2_1/0_0_-1'],
     ['final_state_GX/operator_GammaXGamma5/n2_2/1_1_0','final_state_GY/operator_GammaYGamma5/n2_2/1_1_0','final_state_GZ/operator_GammaZGamma5/n2_2/0_1_1','final_state_GX/operator_GammaXGamma5/n2_2/1_0_1', 'final_state_GY/operator_GammaYGamma5/n2_2/0_1_1','final_state_GZ/operator_GammaZGamma5/n2_2/1_1_0',
      'final_state_GX/operator_GammaXGamma5/n2_2/-1_1_0','final_state_GY/operator_GammaYGamma5/n2_2/-1_1_0','final_state_GZ/operator_GammaZGamma5/n2_2/0_-1_1','final_state_GX/operator_GammaXGamma5/n2_2/-1_0_1', 'final_state_GY/operator_GammaYGamma5/n2_2/0_-1_1','final_state_GZ/operator_GammaZGamma5/n2_2/-1_1_0',
      'final_state_GX/operator_GammaXGamma5/n2_2/1_-1_0','final_state_GY/operator_GammaYGamma5/n2_2/1_-1_0','final_state_GZ/operator_GammaZGamma5/n2_2/0_1_-1','final_state_GX/operator_GammaXGamma5/n2_2/1_0_-1', 'final_state_GY/operator_GammaYGamma5/n2_2/0_1_-1','final_state_GZ/operator_GammaZGamma5/n2_2/1_-1_0',
       'final_state_GX/operator_GammaXGamma5/n2_2/-1_-1_0','final_state_GY/operator_GammaYGamma5/n2_2/-1_-1_0','final_state_GZ/operator_GammaZGamma5/n2_2/0_-1_-1','final_state_GX/operator_GammaXGamma5/n2_2/-1_0_-1', 'final_state_GY/operator_GammaYGamma5/n2_2/0_-1_-1','final_state_GZ/operator_GammaZGamma5/n2_2/-1_-1_0'],
     ['final_state_GX/operator_GammaXGamma5/n2_3/1_1_1','final_state_GY/operator_GammaYGamma5/n2_3/1_1_1','final_state_GZ/operator_GammaZGamma5/n2_3/1_1_1','final_state_GX/operator_GammaXGamma5/n2_3/-1_1_1','final_state_GY/operator_GammaYGamma5/n2_3/-1_1_1','final_state_GZ/operator_GammaZGamma5/n2_3/-1_1_1',
      'final_state_GX/operator_GammaXGamma5/n2_3/1_-1_1','final_state_GY/operator_GammaYGamma5/n2_3/1_-1_1','final_state_GZ/operator_GammaZGamma5/n2_3/1_-1_1','final_state_GX/operator_GammaXGamma5/n2_3/1_1_-1','final_state_GY/operator_GammaYGamma5/n2_3/1_1_-1','final_state_GZ/operator_GammaZGamma5/n2_3/1_1_-1',
      'final_state_GX/operator_GammaXGamma5/n2_3/-1_-1_1','final_state_GY/operator_GammaYGamma5/n2_3/-1_-1_1','final_state_GZ/operator_GammaZGamma5/n2_3/-1_-1_1','final_state_GX/operator_GammaXGamma5/n2_3/-1_1_-1','final_state_GY/operator_GammaYGamma5/n2_3/-1_1_-1','final_state_GZ/operator_GammaZGamma5/n2_3/-1_1_-1',
      'final_state_GX/operator_GammaXGamma5/n2_3/1_-1_-1','final_state_GY/operator_GammaYGamma5/n2_3/1_-1_-1','final_state_GZ/operator_GammaZGamma5/n2_3/1_-1_-1','final_state_GX/operator_GammaXGamma5/n2_3/-1_-1_-1','final_state_GY/operator_GammaYGamma5/n2_3/-1_-1_-1','final_state_GZ/operator_GammaZGamma5/n2_3/-1_-1_-1'],
     ['final_state_GX/operator_GammaXGamma5/n2_4/2_0_0','final_state_GY/operator_GammaYGamma5/n2_4/0_2_0', 'final_state_GZ/operator_GammaZGamma5/n2_4/0_0_2','final_state_GX/operator_GammaXGamma5/n2_4/-2_0_0','final_state_GY/operator_GammaYGamma5/n2_4/0_-2_0', 'final_state_GZ/operator_GammaZGamma5/n2_4/0_0_-2'],
     ['final_state_GX/operator_GammaXGamma5/n2_5/2_1_0','final_state_GY/operator_GammaYGamma5/n2_5/2_1_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_2_1','final_state_GX/operator_GammaXGamma5/n2_5/2_0_1', 'final_state_GY/operator_GammaYGamma5/n2_5/0_2_1','final_state_GZ/operator_GammaZGamma5/n2_5/2_0_1',
      'final_state_GX/operator_GammaXGamma5/n2_5/1_2_0','final_state_GY/operator_GammaYGamma5/n2_5/1_2_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_1_2','final_state_GX/operator_GammaXGamma5/n2_5/1_0_2', 'final_state_GY/operator_GammaYGamma5/n2_5/0_1_2','final_state_GZ/operator_GammaZGamma5/n2_5/1_0_2',
      'final_state_GX/operator_GammaXGamma5/n2_5/-2_1_0','final_state_GY/operator_GammaYGamma5/n2_5/-2_1_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_-2_1','final_state_GX/operator_GammaXGamma5/n2_5/-2_0_1', 'final_state_GY/operator_GammaYGamma5/n2_5/0_-2_1','final_state_GZ/operator_GammaZGamma5/n2_5/-2_0_1',
       'final_state_GX/operator_GammaXGamma5/n2_5/-1_2_0','final_state_GY/operator_GammaYGamma5/n2_5/-1_2_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_-1_2','final_state_GX/operator_GammaXGamma5/n2_5/-1_0_2', 'final_state_GY/operator_GammaYGamma5/n2_5/0_-1_2','final_state_GZ/operator_GammaZGamma5/n2_5/-1_0_2',
      'final_state_GX/operator_GammaXGamma5/n2_5/2_-1_0','final_state_GY/operator_GammaYGamma5/n2_5/2_-1_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_2_-1','final_state_GX/operator_GammaXGamma5/n2_5/2_0_-1', 'final_state_GY/operator_GammaYGamma5/n2_5/0_2_-1','final_state_GZ/operator_GammaZGamma5/n2_5/2_0_-1',
        'final_state_GX/operator_GammaXGamma5/n2_5/1_-2_0','final_state_GY/operator_GammaYGamma5/n2_5/1_-2_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_1_-2','final_state_GX/operator_GammaXGamma5/n2_5/1_0_-2', 'final_state_GY/operator_GammaYGamma5/n2_5/0_1_-2','final_state_GZ/operator_GammaZGamma5/n2_5/1_0_-2',
        'final_state_GX/operator_GammaXGamma5/n2_5/-2_-1_0','final_state_GY/operator_GammaYGamma5/n2_5/-2_-1_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_-2_-1','final_state_GX/operator_GammaXGamma5/n2_5/-2_0_-1', 'final_state_GY/operator_GammaYGamma5/n2_5/0_-2_-1','final_state_GZ/operator_GammaZGamma5/n2_5/-2_0_-1',
         'final_state_GX/operator_GammaXGamma5/n2_5/-1_-2_0','final_state_GY/operator_GammaYGamma5/n2_5/-1_-2_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_-1_-2','final_state_GX/operator_GammaXGamma5/n2_5/-1_0_-2', 'final_state_GY/operator_GammaYGamma5/n2_5/0_-1_-2','final_state_GZ/operator_GammaZGamma5/n2_5/-1_0_-2']]

mom2=[[''],
     ['final_state_GX/operator_GammaXGamma5/n2_1/0_1_0','final_state_GY/operator_GammaYGamma5/n2_1/1_0_0', 'final_state_GZ/operator_GammaZGamma5/n2_1/1_0_0','final_state_GX/operator_GammaXGamma5/n2_1/0_0_1','final_state_GY/operator_GammaYGamma5/n2_1/0_0_1', 'final_state_GZ/operator_GammaZGamma5/n2_1/0_1_0',
      'final_state_GX/operator_GammaXGamma5/n2_1/0_-1_0','final_state_GY/operator_GammaYGamma5/n2_1/-1_0_0', 'final_state_GZ/operator_GammaZGamma5/n2_1/-1_0_0','final_state_GX/operator_GammaXGamma5/n2_1/0_0_-1','final_state_GY/operator_GammaYGamma5/n2_1/0_0_-1', 'final_state_GZ/operator_GammaZGamma5/n2_1/0_-1_0'],
     ['final_state_GX/operator_GammaXGamma5/n2_2/0_1_1','final_state_GY/operator_GammaYGamma5/n2_2/1_0_1','final_state_GZ/operator_GammaZGamma5/n2_2/1_1_0',
      'final_state_GX/operator_GammaXGamma5/n2_2/0_-1_1','final_state_GY/operator_GammaYGamma5/n2_2/-1_0_1','final_state_GZ/operator_GammaZGamma5/n2_2/-1_1_0',
      'final_state_GX/operator_GammaXGamma5/n2_2/0_1_-1','final_state_GY/operator_GammaYGamma5/n2_2/1_0_-1','final_state_GZ/operator_GammaZGamma5/n2_2/1_-1_0',
      'final_state_GX/operator_GammaXGamma5/n2_2/0_-1_-1','final_state_GY/operator_GammaYGamma5/n2_2/-1_0_-1','final_state_GZ/operator_GammaZGamma5/n2_2/-1_-1_0'],
     ['final_state_GX/operator_GammaXGamma5/n2_3/1_1_1','final_state_GY/operator_GammaYGamma5/n2_3/1_1_1','final_state_GZ/operator_GammaZGamma5/n2_3/1_1_1','final_state_GX/operator_GammaXGamma5/n2_3/1_-1_1','final_state_GY/operator_GammaYGamma5/n2_3/1_-1_1','final_state_GZ/operator_GammaZGamma5/n2_3/1_-1_1',
      'final_state_GX/operator_GammaXGamma5/n2_3/-1_1_1','final_state_GY/operator_GammaYGamma5/n2_3/-1_1_1','final_state_GZ/operator_GammaZGamma5/n2_3/-1_1_1','final_state_GX/operator_GammaXGamma5/n2_3/1_1_-1','final_state_GY/operator_GammaYGamma5/n2_3/1_1_-1','final_state_GZ/operator_GammaZGamma5/n2_3/1_1_-1',
      'final_state_GX/operator_GammaXGamma5/n2_3/-1_-1_1','final_state_GY/operator_GammaYGamma5/n2_3/-1_-1_1','final_state_GZ/operator_GammaZGamma5/n2_3/-1_-1_1','final_state_GX/operator_GammaXGamma5/n2_3/-1_1_-1','final_state_GY/operator_GammaYGamma5/n2_3/-1_1_-1','final_state_GZ/operator_GammaZGamma5/n2_3/-1_1_-1',
     'final_state_GX/operator_GammaXGamma5/n2_3/-1_1_-1','final_state_GY/operator_GammaYGamma5/n2_3/-1_1_-1','final_state_GZ/operator_GammaZGamma5/n2_3/-1_1_-1','final_state_GX/operator_GammaXGamma5/n2_3/-1_-1_-1','final_state_GY/operator_GammaYGamma5/n2_3/-1_-1_-1','final_state_GZ/operator_GammaZGamma5/n2_3/-1_-1_-1'],
     ['final_state_GX/operator_GammaXGamma5/n2_4/0_2_0','final_state_GY/operator_GammaYGamma5/n2_4/2_0_0', 'final_state_GZ/operator_GammaZGamma5/n2_4/2_0_0','final_state_GX/operator_GammaXGamma5/n2_4/0_0_2','final_state_GY/operator_GammaYGamma5/n2_4/0_0_2', 'final_state_GZ/operator_GammaZGamma5/n2_4/0_2_0',
      'final_state_GX/operator_GammaXGamma5/n2_4/0_-2_0','final_state_GY/operator_GammaYGamma5/n2_4/-2_0_0', 'final_state_GZ/operator_GammaZGamma5/n2_4/-2_0_0','final_state_GX/operator_GammaXGamma5/n2_4/0_0_-2','final_state_GY/operator_GammaYGamma5/n2_4/0_0_-2', 'final_state_GZ/operator_GammaZGamma5/n2_4/0_-2_0'],
     ['final_state_GX/operator_GammaXGamma5/n2_5/0_2_1','final_state_GY/operator_GammaYGamma5/n2_5/2_0_1','final_state_GZ/operator_GammaZGamma5/n2_5/2_1_0','final_state_GX/operator_GammaXGamma5/n2_5/0_1_2','final_state_GY/operator_GammaYGamma5/n2_5/1_0_2','final_state_GZ/operator_GammaZGamma5/n2_5/1_2_0',
      'final_state_GX/operator_GammaXGamma5/n2_5/0_-2_1','final_state_GY/operator_GammaYGamma5/n2_5/-2_0_1','final_state_GZ/operator_GammaZGamma5/n2_5/-2_1_0','final_state_GX/operator_GammaXGamma5/n2_5/0_-1_2','final_state_GY/operator_GammaYGamma5/n2_5/-1_0_2','final_state_GZ/operator_GammaZGamma5/n2_5/-1_2_0',
      'final_state_GX/operator_GammaXGamma5/n2_5/0_2_-1','final_state_GY/operator_GammaYGamma5/n2_5/2_0_-1','final_state_GZ/operator_GammaZGamma5/n2_5/2_-1_0','final_state_GX/operator_GammaXGamma5/n2_5/0_1_-2','final_state_GY/operator_GammaYGamma5/n2_5/1_0_-2','final_state_GZ/operator_GammaZGamma5/n2_5/1_-2_0']]




ptmom=['0_0_0','1_0_0','1_1_0','1_1_1', '2_0_0','2_1_0']

dsets=[f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]
dsetsb=[f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]

dsets2=[f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom2[nsq][i])] for i in range(len(mom2[nsq]))]
dsetsb2=[f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom2[nsq][i])] for i in range(len(mom2[nsq]))]


nmom=len(mom[nsq])
nmom2=len(mom2[nsq])


dsxn0=f["/CHARM_SM12.14_SM12.14_s0.02144/c0.248/operator_GammaX/{}/data".format(ptmom[nsq])]
dsyn0=f["/CHARM_SM12.14_SM12.14_s0.02144/c0.248/operator_GammaY/{}/data".format(ptmom[nsq])]
dszn0=f["/CHARM_SM12.14_SM12.14_s0.02144/c0.248/operator_GammaZ/{}/data".format(ptmom[nsq])]
bsn0=f["/rhq_m2.42_csw2.68_zeta1.52_SM12.14_SM12.14_s0.02144/operator_Gamma5/0_0_0/data"]

bsfit=pd.read_csv('./2pt/Bs-blocks.csv',sep='\s')
dsfit=pd.read_csv('./2pt/Ds248-nsq{}-blocks.csv'.format(nsq),sep='\s')
ds0fit=pd.read_csv('./2pt/Ds248-nsq0-blocks.csv',sep='\s')


#####Eff Masses
mdlist=[0.7348303222656272,0.7458868408203149,0.7567413330078149,0.7673809814453149,0.77745605,0.787656860351565]
md=mdlist[nsq]

# Constants

mb = 1.92498840332032
md0 = 0.73483032
pre=-md0**2*(mb-md)/(2*(mb**2*md))
#pre2=-(mb+md)
pre2=1/md0**2

nconf=98
dt=30
ts=96

# Initialize arrays

av1n0=np.zeros((nmom, dt+1,nconf))
for i in range(nmom):
    av1n0[i]=np.zeros((dt+1, nconf))

av1n02=np.zeros((nmom2, dt+1,nconf))
for i in range(nmom2):
    av1n02[i]=np.zeros((dt+1, nconf))


avdx = np.zeros((dt+1, nconf))
avdy = np.zeros((dt+1, nconf))
avdz = np.zeros((dt+1, nconf))
avb = np.zeros((dt+1, nconf))
avn0 = np.zeros(dt+1)

findx = np.zeros(dt+1)
findy = np.zeros(dt+1)
findz = np.zeros(dt+1)
finb = np.zeros(dt+1)
fin3pt=np.zeros((dt+1,nmom))

fin3pt2 =np.zeros((dt+1,nmom2))


# Loop optsmizatsons
for j in range(dt+1):
    
    temp=np.zeros(nmom)
    temp2=np.zeros(nmom2)  
    
    tempdx = 0
    tempdy = 0
    tempdz = 0
    tempb = 0
    
    for k in range(nconf):
        
        
        tmp=[np.mean((np.real(dsets[i][k, :, j]) + np.real(dsetsb[i][k, :, dt-j]))) / 2 for i in range(nmom)]
        tmp2=[np.mean((np.real(dsets2[i][k, :, j]) + np.real(dsetsb2[i][k, :, dt-j]))) / 2 for i in range(nmom2)]
        
            
        tmpdx = np.mean(np.real(dsxn0[k, :, j]) + np.real(dsxn0[k, :, ts-j])) / 2 if j != 0 else np.mean(np.real(dsxn0[k, :, 0]))
        tmpdy = np.mean(np.real(dsyn0[k, :, j]) + np.real(dsyn0[k, :, ts-j])) / 2 if j != 0 else np.mean(np.real(dsyn0[k, :, 0]))
        tmpdz = np.mean(np.real(dszn0[k, :, j]) + np.real(dszn0[k, :, ts-j])) / 2 if j != 0 else np.mean(np.real(dszn0[k, :, 0]))
        tmpb = np.mean(np.real(bsn0[k, :, j]) + np.real(bsn0[k, :, ts-j])) / 2 if j != 0 else np.mean(np.real(bsn0[k, :, 0]))

        for l in range(nmom):
            av1n0[l][j,k]=tmp[l]
            temp[l] +=tmp[l]
        for l in range(nmom2):
            av1n02[l][j,k]=tmp2[l]
            temp2[l] +=tmp2[l]
        
        avdx[j, k] = tmpdx
        avdy[j, k] = tmpdy
        avdz[j, k] = tmpdz
        avb[j, k] = tmpb

        tempdx += tmpdx
        tempdy += tmpdy
        tempdz += tmpdz
        tempb += tmpb
    #avn0[j]=pre*(pre2*(temp1xa-temp1ya+temp1za+temp1xma-temp1yma+temp1zma)/nmom2*(1/md-1)*(mb+md0)-((((temp1x-temp1y+temp1z+temp1xm-temp1ym+temp1zm+temp1xn-temp1yn+temp1zn+temp1xmn-temp1ymn+temp1zmn)/nmom)/(1/3*np.sqrt((tempdx+tempdy+tempdz)*tempb))))*np.sqrt((4*mb*md)/(np.exp(-md*j)*np.exp(-mb*(dt-j)))))
    findx[j]=tempdx
    findy[j]=tempdy
    findz[j]=tempdz
    finb[j]=tempb
    fin3pt[j]=temp
    fin3pt2[j]=temp2
    
for j in range(dt+1):    
    avn0[j]=pre*(pre2*sum_with_exceptions(fin3pt2[j])+sum_with_exceptions(fin3pt[j]))/(np.sqrt(1/3*(findx[j]+findy[j]+findz[j])*finb[dt-j]))*np.sqrt((4*mb*md)/(np.exp(-md*j)*np.exp(-mb*(dt-j))))
    

errn0=np.zeros(shape=(dt+1))
#errnn0=np.zeros(shape=(96))

for j in range(dt+1):
    x=0
    for i in range(nconf):
        x=x+(pre*(pre2*sum_with_exceptions_jack(av1n02, j, i)+sum_with_exceptions_jack(av1n0, j, i))/(np.sqrt(1/3*(jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[dt-j],i)))*np.sqrt((4*dsfit['EffectiveMass'][i]*bsfit['EffectiveMass'][i])/(np.exp(-dsfit['EffectiveMass'][i]*j)*np.exp(-bsfit['EffectiveMass'][i]*(dt-j))))-avn0[j])**2
    errn0[j]=np.sqrt((nconf-1)/nconf*x) 
    #errnn0[j]=np.sqrt((98-1)/98*x)*10**(85)

avn0[np.isnan(avn0)] = 0
errn0[np.isnan(errn0)] = 0

#plt.plot(list(range(96)), np.absolute(avn0)[0:96],'ro')
plt.xlabel('time')
plt.ylabel(r'$\widetilde{A}_2$')
plt.errorbar(list(range(dt)), np.absolute(avn0)[0:dt], yerr=errn0[0:dt],ls='none',fmt='x',label='nsq={}'.format(nsq))
#plt.yscale('log')
plt.legend()
plt.savefig('./A2/A2-nsq{}.png'.format(nsq))

np.savetxt('./A2/A2-nsq{}.txt'.format(nsq), np.c_[np.absolute(avn0), errn0])

###############################################################################

reg_low=19
reg_up=26

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
            x=x+(((pre2*sum_with_exceptions_jack(av1n02, t1+reg_low, i)+sum_with_exceptions_jack(av1n0, t1+reg_low, i))/(np.sqrt(1/3*(jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[dt-(t1+reg_low)],i))))*np.sqrt((4*dsfit['EffectiveMass'][i]*bsfit['EffectiveMass'][i])/(np.exp(-dsfit['EffectiveMass'][i]*(t1+reg_low))*np.exp(-bsfit['EffectiveMass'][i]*(dt-(t1+reg_low)))))*pre-avn0[t1])*(((pre2*sum_with_exceptions_jack(av1n02, t2+reg_low, i)+sum_with_exceptions_jack(av1n0, t2+reg_low, i))/(np.sqrt(1/3*(jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[dt-(t2+reg_low)],i))))*np.sqrt((4*dsfit['EffectiveMass'][i]*bsfit['EffectiveMass'][i])/(np.exp(-dsfit['EffectiveMass'][i]*(t2+reg_low))*np.exp(-bsfit['EffectiveMass'][i]*(dt-(t2+reg_low)))))*pre-avn0[t2])            

            '''
            if nmom==3 and nmom2==3:
                x=x+(((pre2*(jack(av1n0xa[t1+reg_low],i)-jack(av1n0ya[t1+reg_low],i)+jack(av1n0za[t1+reg_low],i))/3*(1/md-1)*(mb+md0))-(((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i))/3)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(30-(t1+reg_low))))))*pre-avn0[t1])*(((pre2*(jack(av1n0xa[t2+reg_low],i)-jack(av1n0ya[t2+reg_low],i)+jack(av1n0za[t2+reg_low],i))/3*(1/md-1)*(mb+md0))-(((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i))/3)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low))))))*pre-avn0[t2])            
            elif nmom==3 and nmom2==6:
                x=x+(((pre2*(jack(av1n0xa[t1+reg_low],i)-jack(av1n0ya[t1+reg_low],i)+jack(av1n0za[t1+reg_low],i)+jack(av1n0xma[t1+reg_low],i)-jack(av1n0yma[t1+reg_low],i)+jack(av1n0zma[t1+reg_low],i))/nmom2*(1/md-1)*(mb+md0))-(((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low))))))*pre-avn0[t1])*(((pre2*(jack(av1n0xa[t2+reg_low],i)-jack(av1n0ya[t2+reg_low],i)+jack(av1n0za[t2+reg_low],i)+jack(av1n0xma[t2+reg_low],i)-jack(av1n0yma[t2+reg_low],i)+jack(av1n0zma[t2+reg_low],i))/6*(1/md-1)*(mb+md0))-(((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i))/3)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low))))))*pre-avn0[t2])            
            elif nmom==6 and nmom2==3:
                x=x+(((pre2*(jack(av1n0xa[t1+reg_low],i)-jack(av1n0ya[t1+reg_low],i)+jack(av1n0za[t1+reg_low],i))/3*(1/md-1)*(mb+md0))-(((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)+jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)+jack(av1n0zm[t1+reg_low],i))/6)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low))))))*pre-avn0[t1])*(((pre2*(jack(av1n0xa[t2+reg_low],i)-jack(av1n0ya[t2+reg_low],i)+jack(av1n0za[t2+reg_low],i))/3*(1/md-1)*(mb+md0))-(((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)+jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)+jack(av1n0zm[t2+reg_low],i))/6)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low))))))*pre-avn0[t2])            
            elif nmom==12 and nmom2==6:
                x=x+(((pre2*(jack(av1n0xa[t1+reg_low],i)-jack(av1n0ya[t1+reg_low],i)+jack(av1n0za[t1+reg_low],i)+jack(av1n0xma[t1+reg_low],i)-jack(av1n0yma[t1+reg_low],i)+jack(av1n0zma[t1+reg_low],i))/nmom2*(1/md-1)*(mb+md0))-(((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)+jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)+jack(av1n0zm[t1+reg_low],i)+jack(av1n0xn[t1+reg_low],i)-jack(av1n0yn[t1+reg_low],i)+jack(av1n0zn[t1+reg_low],i)+jack(av1n0xmn[t1+reg_low],i)-jack(av1n0ymn[t1+reg_low],i)+jack(av1n0zmn[t1+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low))))))*pre-avn0[t1])*(((pre2*(jack(av1n0xa[t2+reg_low],i)-jack(av1n0ya[t2+reg_low],i)+jack(av1n0za[t2+reg_low],i)+jack(av1n0xma[t2+reg_low],i)-jack(av1n0yma[t2+reg_low],i)+jack(av1n0zma[t2+reg_low],i))/6*(1/md-1)*(mb+md0))-(((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)+jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)+jack(av1n0zm[t2+reg_low],i)+jack(av1n0xn[t2+reg_low],i)-jack(av1n0yn[t2+reg_low],i)+jack(av1n0zn[t2+reg_low],i)+jack(av1n0xmn[t2+reg_low],i)-jack(av1n0ymn[t2+reg_low],i)+jack(av1n0zmn[t2+reg_low],i))/3)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low))))))*pre-avn0[t2])            
            '''
        covmat[t1][t2]=x
        covmat[t2][t1]=x  
        
        
def chi(a):
    return (nconf-1-reg_low-cut)/(nconf-reg_low-cut)*np.dot(np.transpose([i-a for i in avn0[reg_low:reg_up]]),np.matmul(np.linalg.inv(covmat),[i-a for i in avn0[reg_low:reg_up]]))

mbar=minimize(chi,0.1,method='Nelder-Mead', tol=1e-6)


def jackmass(t1,i):
    #+jack(av1n0xm[t1],i)-jack(av1n0ym[t1],i)+jack(av1n0zm[t1],i)
    return pre*((pre2*sum_with_exceptions_jack(av1n02, t1, i)+sum_with_exceptions_jack(av1n0, t1, i))/(np.sqrt(1/3*(jack(avdx[t1],i)+jack(avdy[t1],i)+jack(avdz[t1],i))*jack(avb[dt-t1],i))))*np.sqrt((4*dsfit['EffectiveMass'][i]*bsfit['EffectiveMass'][i])/(np.exp(-dsfit['EffectiveMass'][i]*(t1))*np.exp(-bsfit['EffectiveMass'][i]*(dt-(t1)))))

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

ax.set_xlabel('Time')
ax.set_ylabel('A2')
ax.errorbar(list(range(dt))[1:dt], avn0[1:dt], yerr=errn0[1:dt],fmt='x', label='nsq={}'.format(nsq))
plt.axhline(y = mbar.x[0], color = 'r', linestyle = 'dashed')
plt.fill_between(list(range(dt))[reg_low+1:reg_up+1], mbar.x[0]+sigma, mbar.x[0]-sigma, color='r',alpha=0.2)
#ax.set_yscale('log')
plt.savefig('./Fits/A2-Av-nsq{}-Fit.png'.format(nsq))

df3 = pd.DataFrame(columns=['EffectiveMass','Error','RegUp','RegLow'])
df3['EffectiveMass']=mbar.x
df3['Error']=sigma  
df3['RegUp']=reg_up
df3['RegLow']=reg_low    
df3.to_csv('./Fits/A2-Av-nsq{}-Fit.csv'.format(nsq), sep='\t')






'''
# nsq=1
dset1xn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GX/operator_GammaXGamma5/n2_1/1_0_0/forward/data"]
dset1xn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GX/operator_GammaXGamma5/n2_1/1_0_0/backward/data"]

dset1yn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GY/operator_GammaYGamma5/n2_1/0_1_0/forward/data"]
dset1yn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GY/operator_GammaYGamma5/n2_1/0_1_0/backward/data"]

dset1zn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GZ/operator_GammaZGamma5/n2_1/0_0_1/forward/data"]
dset1zn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GZ/operator_GammaZGamma5/n2_1/0_0_1/backward/data"]

md=0.7458868408203149
'''

'''
#nsq2
dset1xn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GX/operator_GammaXGamma5/n2_2/1_1_0/forward/data"]
dset1xn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GX/operator_GammaXGamma5/n2_2/1_1_0/backward/data"]

dset1yn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GY/operator_GammaYGamma5/n2_2/1_1_0/forward/data"]
dset1yn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GY/operator_GammaYGamma5/n2_2/1_1_0/backward/data"]

dset1zn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GZ/operator_GammaZGamma5/n2_2/0_1_1/forward/data"]
dset1zn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GZ/operator_GammaZGamma5/n2_2/0_1_1/backward/data"]

dset1xn0m = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GX/operator_GammaXGamma5/n2_2/1_0_1/forward/data"]
dset1xn0bm = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GX/operator_GammaXGamma5/n2_2/1_0_1/backward/data"]

dset1yn0m = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GY/operator_GammaYGamma5/n2_2/0_1_1/forward/data"]
dset1yn0bm = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GY/operator_GammaYGamma5/n2_2/0_1_1/backward/data"]

dset1zn0m = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GZ/operator_GammaZGamma5/n2_2/1_1_0/forward/data"]
dset1zn0bm = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GZ/operator_GammaZGamma5/n2_2/1_1_0/backward/data"]

md=0.7567413330078149
'''

'''
#nsq3
dset1xn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GX/operator_GammaXGamma5/n2_3/1_1_1/forward/data"]
dset1xn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GX/operator_GammaXGamma5/n2_3/1_1_1/backward/data"]

dset1yn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GY/operator_GammaYGamma5/n2_3/1_1_1/forward/data"]
dset1yn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GY/operator_GammaYGamma5/n2_3/1_1_1/backward/data"]

dset1zn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GZ/operator_GammaZGamma5/n2_3/1_1_1/forward/data"]
dset1zn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GZ/operator_GammaZGamma5/n2_3/1_1_1/backward/data"]

md=0.7673809814453149
'''
'''
#nsq4
dset1xn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GX/operator_GammaXGamma5/n2_4/2_0_0/forward/data"]
dset1xn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GX/operator_GammaXGamma5/n2_4/2_0_0/backward/data"]

dset1yn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GY/operator_GammaYGamma5/n2_4/0_2_0/forward/data"]
dset1yn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GY/operator_GammaYGamma5/n2_4/0_2_0/backward/data"]

dset1zn0 = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GZ/operator_GammaZGamma5/n2_4/0_0_2/forward/data"]
dset1zn0b = f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/final_state_GZ/operator_GammaZGamma5/n2_4/0_0_2/backward/data"]

md=0.77745605

mom2=[[''],
     ['final_state_GX/operator_GammaXGamma5/n2_1/0_1_0','final_state_GY/operator_GammaYGamma5/n2_1/1_0_0', 'final_state_GZ/operator_GammaZGamma5/n2_1/1_0_0','final_state_GX/operator_GammaXGamma5/n2_1/0_0_1','final_state_GY/operator_GammaYGamma5/n2_1/0_0_1', 'final_state_GZ/operator_GammaZGamma5/n2_1/0_1_0'],
     ['final_state_GX/operator_GammaXGamma5/n2_2/0_1_1','final_state_GY/operator_GammaYGamma5/n2_2/1_0_1','final_state_GZ/operator_GammaZGamma5/n2_2/1_1_0'],
     ['final_state_GX/operator_GammaXGamma5/n2_3/1_1_1','final_state_GY/operator_GammaYGamma5/n2_3/1_1_1','final_state_GZ/operator_GammaZGamma5/n2_3/1_1_1'],
     ['final_state_GX/operator_GammaXGamma5/n2_4/0_2_0','final_state_GY/operator_GammaYGamma5/n2_4/2_0_0', 'final_state_GZ/operator_GammaZGamma5/n2_4/2_0_0','final_state_GX/operator_GammaXGamma5/n2_4/0_0_2','final_state_GY/operator_GammaYGamma5/n2_4/0_0_2', 'final_state_GZ/operator_GammaZGamma5/n2_4/0_2_0'],
     ['final_state_GX/operator_GammaXGamma5/n2_5/0_2_1','final_state_GY/operator_GammaYGamma5/n2_5/2_0_1','final_state_GZ/operator_GammaZGamma5/n2_5/2_1_0','final_state_GX/operator_GammaXGamma5/n2_5/0_1_2','final_state_GY/operator_GammaYGamma5/n2_5/1_0_2','final_state_GZ/operator_GammaZGamma5/n2_5/1_2_0']]

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


if len(mom2[nsq])==3:
    dset1xn0a=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom2[nsq][0])]
    dset1xn0ba=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom2[nsq][0])]
    dset1yn0a=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom2[nsq][1])]
    dset1yn0ba=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom2[nsq][1])]
    dset1zn0a=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom2[nsq][2])]
    dset1zn0ba=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom2[nsq][2])]
    nmom2=3
elif len(mom2[nsq])==6:
    dset1xn0a=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom2[nsq][0])]
    dset1xn0ba=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom2[nsq][0])]
    dset1yn0a=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom2[nsq][1])]
    dset1yn0ba=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom2[nsq][1])]
    dset1zn0a=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom2[nsq][2])]
    dset1zn0ba=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom2[nsq][2])]   
    dset1xn0ma=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom2[nsq][3])]
    dset1xn0bma=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom2[nsq][3])]
    dset1yn0ma=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom2[nsq][4])]
    dset1yn0bma=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom2[nsq][4])]
    dset1zn0ma=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom2[nsq][5])]
    dset1zn0bma=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom2[nsq][5])]   
    nmom2=6


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

av1n0xa = np.zeros((dt, nconf))
av1n0ya = np.zeros((dt, nconf))
av1n0za = np.zeros((dt, nconf))
av1n0xma = np.zeros((dt, nconf))
av1n0yma = np.zeros((dt, nconf))
av1n0zma = np.zeros((dt, nconf))

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

temp1xa = 0
temp1ya = 0
temp1za = 0
temp1xma = 0
temp1yma = 0
temp1zma = 0

tmp1x = np.mean((np.real(dset1xn0[k, :, j]) + np.real(dset1xn0b[k, :, dt-j]))) / 2
tmp1y = np.mean((np.real(dset1yn0[k, :, j]) + np.real(dset1yn0b[k, :, dt-j]))) / 2
tmp1z = np.mean((np.real(dset1zn0[k, :, j]) + np.real(dset1zn0b[k, :, dt-j]))) / 2

tmp1xa = np.mean((np.real(dset1xn0a[k, :, j]) + np.real(dset1xn0ba[k, :, dt-j]))) / 2
tmp1ya = np.mean((np.real(dset1yn0a[k, :, j]) + np.real(dset1yn0ba[k, :, dt-j]))) / 2
tmp1za = np.mean((np.real(dset1zn0a[k, :, j]) + np.real(dset1zn0ba[k, :, dt-j]))) / 2

tmp1xm=0
tmp1ym=0
tmp1zm=0
tmp1xma=0
tmp1yma=0
tmp1zma=0

tmp1xn=0
tmp1yn=0
tmp1zn=0
tmp1xmn=0
tmp1ymn=0
tmp1zmn=0

if nmom==6:
    tmp1xm = np.mean((np.real(dset1xn0m[k, :, j]) + np.real(dset1xn0bm[k, :, dt-j]))) / 2
    tmp1ym = np.mean((np.real(dset1yn0m[k, :, j]) + np.real(dset1yn0bm[k, :, dt-j]))) / 2
    tmp1zm = np.mean((np.real(dset1zn0m[k, :, j]) + np.real(dset1zn0bm[k, :, dt-j]))) / 2
elif nmom==12:
    tmp1xm = np.mean((np.real(dset1xn0m[k, :, j]) + np.real(dset1xn0bm[k, :, dt-j]))) / 2
    tmp1ym = np.mean((np.real(dset1yn0m[k, :, j]) + np.real(dset1yn0bm[k, :, dt-j]))) / 2
    tmp1zm = np.mean((np.real(dset1zn0m[k, :, j]) + np.real(dset1zn0bm[k, :, dt-j]))) / 2
    
    tmp1xn = np.mean((np.real(dset1xn0n[k, :, j]) + np.real(dset1xn0bn[k, :, dt-j]))) / 2
    tmp1yn = np.mean((np.real(dset1yn0n[k, :, j]) + np.real(dset1yn0bn[k, :, dt-j]))) / 2
    tmp1zn = np.mean((np.real(dset1zn0n[k, :, j]) + np.real(dset1zn0bn[k, :, dt-j]))) / 2
    tmp1xmn = np.mean((np.real(dset1xn0mn[k, :, j]) + np.real(dset1xn0bmn[k, :, dt-j]))) / 2
    tmp1ymn = np.mean((np.real(dset1yn0mn[k, :, j]) + np.real(dset1yn0bmn[k, :, dt-j]))) / 2
    tmp1zmn = np.mean((np.real(dset1zn0mn[k, :, j]) + np.real(dset1zn0bmn[k, :, dt-j]))) / 2 
    
if nmom2==6:
    tmp1xma = np.mean((np.real(dset1xn0ma[k, :, j]) + np.real(dset1xn0bma[k, :, dt-j]))) / 2
    tmp1yma = np.mean((np.real(dset1yn0ma[k, :, j]) + np.real(dset1yn0bma[k, :, dt-j]))) / 2
    tmp1zma = np.mean((np.real(dset1zn0ma[k, :, j]) + np.real(dset1zn0bma[k, :, dt-j]))) / 2
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
    


'''
