import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy.optimize import minimize
import scipy


def pvalue(chi2, dof):
    r"""Compute the $p$-value corresponding to a $\chi^2$ with `dof` degrees
    of freedom."""
    return 1 - scipy.stats.chi2.cdf(chi2, dof)

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

def parse_string(input_string):
    stmp=input_string.split('/')
    fs=stmp[0].split('_G')[1]
    op=stmp[1].split('_Gamma')[1]
    mom=stmp[3]
    return fs,op,mom
    

def eps_over_k(fs,op,mom):
    cycl=['X','Y','Z']
    pos=['012','120','201']
    neg=['021','210','102']
    iV=cycl.index(fs)
    iJ=cycl.index(op)
    if(iV==iJ): return 0.
    iMom=3-iV-iJ
    kComp=(float)(mom.split('_')[iMom])
    if(kComp==0.): return 0.
    epsStr=str(iV)+str(iJ)+str(iMom)
    if(epsStr in pos): sign=1.
    elif(epsStr in neg): sign=-1.
    else: sign=0. 
    return sign/kComp

def sum_comps(lst):
    total=0
    for i in range(len(lst)):
        a,b,c=parse_string(mom[nsq][i])
        total += eps_over_k(a,b,c)*lst[i]
    return total/len(lst)




#########decide here which nsq
nsq=5
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
      'final_state_GX/operator_GammaXGamma5/n2_3/1_-1_-1','final_state_GY/operator_GammaYGamma5/n2_3/1_-1_-1','final_state_GZ/operator_GammaZGamma5/n2_3/1_-1_-1','final_state_GX/operator_GammaXGamma5/n2_3/-1_-1_-1'],
     ['final_state_GX/operator_GammaXGamma5/n2_4/2_0_0','final_state_GY/operator_GammaYGamma5/n2_4/0_2_0', 'final_state_GZ/operator_GammaZGamma5/n2_4/0_0_2','final_state_GX/operator_GammaXGamma5/n2_4/-2_0_0','final_state_GY/operator_GammaYGamma5/n2_4/0_-2_0', 'final_state_GZ/operator_GammaZGamma5/n2_4/0_0_-2'],
     ['final_state_GX/operator_GammaXGamma5/n2_5/2_1_0','final_state_GY/operator_GammaYGamma5/n2_5/2_1_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_2_1','final_state_GX/operator_GammaXGamma5/n2_5/2_0_1', 'final_state_GY/operator_GammaYGamma5/n2_5/0_2_1','final_state_GZ/operator_GammaZGamma5/n2_5/2_0_1',
      'final_state_GX/operator_GammaXGamma5/n2_5/1_2_0','final_state_GY/operator_GammaYGamma5/n2_5/1_2_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_1_2','final_state_GX/operator_GammaXGamma5/n2_5/1_0_2', 'final_state_GY/operator_GammaYGamma5/n2_5/0_1_2','final_state_GZ/operator_GammaZGamma5/n2_5/1_0_2',
      'final_state_GX/operator_GammaXGamma5/n2_5/-2_1_0','final_state_GY/operator_GammaYGamma5/n2_5/-2_1_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_-2_1','final_state_GX/operator_GammaXGamma5/n2_5/-2_0_1', 'final_state_GY/operator_GammaYGamma5/n2_5/0_-2_1','final_state_GZ/operator_GammaZGamma5/n2_5/-2_0_1',
       'final_state_GX/operator_GammaXGamma5/n2_5/-1_2_0','final_state_GY/operator_GammaYGamma5/n2_5/-1_2_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_-1_2','final_state_GX/operator_GammaXGamma5/n2_5/-1_0_2', 'final_state_GY/operator_GammaYGamma5/n2_5/0_-1_2','final_state_GZ/operator_GammaZGamma5/n2_5/-1_0_2',
      'final_state_GX/operator_GammaXGamma5/n2_5/2_-1_0','final_state_GY/operator_GammaYGamma5/n2_5/2_-1_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_2_-1','final_state_GX/operator_GammaXGamma5/n2_5/2_0_-1', 'final_state_GY/operator_GammaYGamma5/n2_5/0_2_-1','final_state_GZ/operator_GammaZGamma5/n2_5/2_0_-1',
        'final_state_GX/operator_GammaXGamma5/n2_5/1_-2_0','final_state_GY/operator_GammaYGamma5/n2_5/1_-2_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_1_-2','final_state_GX/operator_GammaXGamma5/n2_5/1_0_-2', 'final_state_GY/operator_GammaYGamma5/n2_5/0_1_-2','final_state_GZ/operator_GammaZGamma5/n2_5/1_0_-2',
        'final_state_GX/operator_GammaXGamma5/n2_5/-2_-1_0','final_state_GY/operator_GammaYGamma5/n2_5/-2_-1_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_-2_-1','final_state_GX/operator_GammaXGamma5/n2_5/-2_0_-1', 'final_state_GY/operator_GammaYGamma5/n2_5/0_-2_-1','final_state_GZ/operator_GammaZGamma5/n2_5/-2_0_-1',
         'final_state_GX/operator_GammaXGamma5/n2_5/-1_-2_0','final_state_GY/operator_GammaYGamma5/n2_5/-1_-2_0','final_state_GZ/operator_GammaZGamma5/n2_5/0_-1_-2','final_state_GX/operator_GammaXGamma5/n2_5/-1_0_-2', 'final_state_GY/operator_GammaYGamma5/n2_5/0_-1_-2','final_state_GZ/operator_GammaZGamma5/n2_5/-1_0_-2']]

ptmom=['0_0_0','1_0_0','1_1_0','1_1_1', '2_0_0','2_1_0']


dsets=[f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]
dsetsb=[f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]

nmom=len(mom[nsq])

dsxn0=f["/CHARM_SM12.14_SM12.14_s0.02144/c0.248/operator_GammaX/{}/data".format(ptmom[nsq])]
dsyn0=f["/CHARM_SM12.14_SM12.14_s0.02144/c0.248/operator_GammaY/{}/data".format(ptmom[nsq])]
dszn0=f["/CHARM_SM12.14_SM12.14_s0.02144/c0.248/operator_GammaZ/{}/data".format(ptmom[nsq])]
bsn0=f["/rhq_m2.42_csw2.68_zeta1.52_SM12.14_SM12.14_s0.02144/operator_Gamma5/0_0_0/data"]

bsfit=pd.read_csv('./2pt/Bs-blocks.csv',sep='\s')
dsfit=pd.read_csv('./2pt/Ds248-nsq{}-blocks.csv'.format(nsq),sep='\s')


#####Eff Masses
mdlist=['',0.7458868408203149,0.7567413330078149,0.7673809814453149,0.77745605,0.787656860351565]
md=mdlist[nsq]

# Constants

mb = 1.92498840332032
md0 = 0.73483032
pre = md0 / (2 * md * mb)

nconf=98
dt=30
ts=96

# Initialize arrays

av1n0=np.zeros((nmom, dt,nconf))
for i in range(nmom):
    av1n0[i]=np.zeros((dt, nconf))

# Initialize arrays
av1n01=np.zeros((nmom, dt,nconf))
av1n01back=np.zeros((nmom, dt,nconf))
for i in range(nmom):
    av1n01[i]=np.zeros((dt, nconf))
    av1n01back[i]=np.zeros((dt, nconf))
    

avdx = np.zeros((dt, nconf))
avdy = np.zeros((dt, nconf))
avdz = np.zeros((dt, nconf))
avb = np.zeros((dt, nconf))
avn01 = np.zeros((nmom,dt))
avn01back = np.zeros((nmom,dt))


# Loop optsmizatsons
for j in range(dt):
    
    temp1=np.zeros(nmom)
    temp1back=np.zeros(nmom)
    
    tempdx = 0
    tempdy = 0
    tempdz = 0
    tempb = 0
    
    for k in range(nconf):
        
        #tmp=[np.mean((np.imag(dsets[i][k, :, j]) - np.imag(dsetsb[i][k, :, dt-j]))) / 2 for i in range(nmom)]
    
        tmp1=[np.mean((np.real(dsets[i][k, :, j])))  for i in range(nmom)]
        tmp1back=[np.mean((np.real(dsetsb[i][k, :, dt-j])))  for i in range(nmom)]
    

        for l in range(3):
            av1n01[l][j,k]=tmp1[l]
            av1n01back[l][j,k]=tmp1back[l]
            temp1[l] +=tmp1[l]
            temp1back[l] +=tmp1back[l]
        
    #avn0[j] = pre * (sum_with_exceptions(temp, nsq) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    for l in range(5):
        avn01[l][j]=temp1[l]/nconf
        avn01back[l][j]=temp1back[l]/nconf
    #print(sum_with_exceptions(temp1, nsq))
    print(temp1[1])
    #avn0[j] = pre * ((temp[0]+temp[1]+temp[2]-temp[3]-temp[4]-temp[5])/6 / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    
    #+temp1xn + temp1yn + temp1zn+temp1xmn + temp1ymn + temp1zmn


errn01=np.zeros(shape=(nmom,dt))
errn01back=np.zeros(shape=(nmom,dt))

for l in range(5):
    for j in range(dt):
        x=0
        xb=0
        for i in range(nconf):
            x=x+(jack(av1n01[l][j],i)-avn01[l][j])**2
            xb=xb+(jack(av1n01back[l][j],i)-avn01back[l][j])**2
        errn01[l,j]=np.sqrt((nconf-1)/nconf*x) 
        errn01back[l,j]=np.sqrt((nconf-1)/nconf*xb) 

for l in range(5):
    plt.errorbar(list(range(dt)),-avn01[:][l],yerr=errn01[l],fmt='x')
    #plt.errorbar(list(range(dt)),avn01back[:][l],yerr=errn01back[l],fmt='x')
plt.yscale('log')

