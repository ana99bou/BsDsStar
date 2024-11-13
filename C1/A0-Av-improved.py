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
#########decide here which nsq
nsq=1
##########


f = h5py.File("BsDsStar_C1.h5", "r")
f2p = h5py.File("BsDsStar_2ptData_C1.h5", "r")
f2pB = h5py.File("BsDsStar-2ptBs_C1.h5", "r")

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


dsets=[f["/CHARM_PT_SEQ_SM7.86_s0.03224/c0.400/dT20/{}/forward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]
dsetsb=[f["/CHARM_PT_SEQ_SM7.86_s0.03224/c0.400/dT20/{}/backward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]

nmom=len(mom[nsq])

dsxn0=f2p["/cl_SM7.86_SM7.86_0.03224/c0.400/operator_GammaX/n2_{}/data".format(nsq)]
dsyn0=f2p["/cl_SM7.86_SM7.86_0.03224/c0.400/operator_GammaY/n2_{}/data".format(nsq)]
dszn0=f2p["/cl_SM7.86_SM7.86_0.03224/c0.400/operator_GammaZ/n2_{}/data".format(nsq)]
bsn0=f2pB["/hl_SM7.86_SM7.86_0.03224_m7.47_csw4.92_zeta2.93/operator_Gamma5/n2_0/data"]



#####Eff Masses
mdlist=[1.1056561279296913,1.1354833984375037,1.165108642578129,1.194660644531254,1.2244628906250041,1.254152832031254]
md=mdlist[nsq]

# Constants

mb = 3.009638061523448
md0 = 1.1056561279296913
pre = md0 / (2 * md * mb)

nconf=1636
dt=20
ts=64

# Initialize arrays

av1n0=np.zeros((nmom, dt,nconf))
for i in range(nmom):
    av1n0[i]=np.zeros((dt, nconf))


avdx = np.zeros((dt, nconf))
avdy = np.zeros((dt, nconf))
avdz = np.zeros((dt, nconf))
avb = np.zeros((dt, nconf))
avn0 = np.zeros(dt)


# Loop optsmizatsons
for j in range(dt):
    
    temp=np.zeros(nmom)
    
    tempdx = 0
    tempdy = 0
    tempdz = 0
    tempb = 0
    
    for k in range(nconf):
        
        tmp=[np.mean((np.real(dsets[i][k, :, j]) + np.real(dsetsb[i][k, :, dt-j]))) / 2 for i in range(nmom)]
            
            
        tmpdx = np.mean(np.real(dsxn0[k, :, j+1]) + np.real(dsxn0[k, :, ts-1-j])) / 2 if j != 0 else np.mean(np.real(dsxn0[k, :, 0]))
        tmpdy = np.mean(np.real(dsyn0[k, :, j+1]) + np.real(dsyn0[k, :, ts-1-j])) / 2 if j != 0 else np.mean(np.real(dsyn0[k, :, 0]))
        tmpdz = np.mean(np.real(dszn0[k, :, j+1]) + np.real(dszn0[k, :, ts-1-j])) / 2 if j != 0 else np.mean(np.real(dszn0[k, :, 0]))
        tmpb = np.real(bsn0[k, dt-1-j]) + np.real(bsn0[k, ts-dt+1+j]) / 2 if j != dt-1 else np.real(bsn0[k,  0])

        for l in range(nmom):
            av1n0[l][j,k]=tmp[l]
            temp[l] +=tmp[l]
        
        avdx[j, k] = tmpdx
        avdy[j, k] = tmpdy
        avdz[j, k] = tmpdz
        avb[j, k] = tmpb

        tempdx += tmpdx
        tempdy += tmpdy
        tempdz += tmpdz
        tempb += tmpb
    #avn0[j] = pre * (((temp1x - temp1y + temp1z+temp1xm - temp1ym + temp1zm+temp1xn - temp1yn + temp1zn+temp1xmn - temp1ymn + temp1zmn+temp1xp - temp1yp + temp1zp+temp1xmp - temp1ymp + temp1zmp+temp1xnp - temp1ynp + temp1znp+temp1xmnp - temp1ymnp + temp1zmnp) / nmom) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    avn0[j] = pre * (((sum_with_exceptions(temp))) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    #avn0[j] = pre * (((temp[0]+temp[1]+temp[2]+temp[3]+temp[4]+temp[5]) / nmom) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    

errn0=np.zeros(shape=(dt))
#errnn0=np.zeros(shape=(96))

for j in range(dt):
    x=0
    for i in range(nconf):
#        x=x+((((jack(av1n0x[j],i)-jack(av1n0y[j],i)+jack(av1n0z[j],i)+jack(av1n0xm[j],i)-jack(av1n0ym[j],i)+jack(av1n0zm[j],i)+jack(av1n0xn[j],i)-jack(av1n0yn[j],i)+jack(av1n0zn[j],i)+jack(av1n0xmn[j],i)-jack(av1n0ymn[j],i)+jack(av1n0zmn[j],i)+jack(av1n0xp[j],i)-jack(av1n0yp[j],i)+jack(av1n0zp[j],i)+jack(av1n0xmp[j],i)-jack(av1n0ymp[j],i)+jack(av1n0zmp[j],i)+jack(av1n0xnp[j],i)-jack(av1n0ynp[j],i)+jack(av1n0znp[j],i)+jack(av1n0xmnp[j],i)-jack(av1n0ymnp[j],i)+jack(av1n0zmnp[j],i))/nmom)/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(30-j))))*pre-avn0[j])**2
        #x=x+((((sum_with_exceptions_jack(av1n0, j, i))/nmom)/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(30-j))))*pre-avn0[j])**2
        #x=x+((((jack(av1n0[0][j],i)+jack(av1n0[1][j],i)+jack(av1n0[2][j],i)+jack(av1n0[3][j],i)+jack(av1n0[4][j],i)+jack(av1n0[5][j],i))/nmom)/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(30-j))))*pre-avn0[j])**2
        x=x+((((sum_with_exceptions_jack(av1n0, j, i)))/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(30-j))))*pre-avn0[j])**2

    errn0[j]=np.sqrt((nconf-1)/nconf*x) 
    #errnn0[j]=np.sqrt((98-1)/98*x)*10**(85)

avn0[np.isnan(avn0)] = 0
errn0[np.isnan(errn0)] = 0

#plt.plot(list(range(96)), np.absolute(avn0)[0:96],'ro')
plt.xlabel('time')
plt.ylabel(r'$\widetilde{A}_0$')
plt.errorbar(list(range(dt)), np.absolute(avn0)[0:dt], yerr=errn0[0:dt],ls='none',fmt='x',label='nsq={}'.format(nsq))
#plt.yscale('log')
plt.legend()
plt.savefig('./A0/A0-nsq{}.png'.format(nsq))

np.savetxt('./A0/A0-nsq{}.txt'.format(nsq), np.c_[np.absolute(avn0), errn0])

###############################################################################

reg_low=13
reg_up=16

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
            x=x+(((sum_with_exceptions_jack(av1n0, t1+reg_low, i))/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*(((sum_with_exceptions_jack(av1n0, t2+reg_low, i))/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
            '''                                                                                                                                                                                                                                                              
            if nmom==3:
                x=x+((((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i))/3)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i))/3)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
            elif nmom==6:
                x=x+((((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)+jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)+jack(av1n0zm[t1+reg_low],i))/6)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)+jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)+jack(av1n0zm[t2+reg_low],i))/6)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
            elif nmom==16:
                x=x+((((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)+jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)+jack(av1n0zm[t1+reg_low],i)+jack(av1n0xn[t1+reg_low],i)-jack(av1n0yn[t1+reg_low],i)+jack(av1n0zn[t1+reg_low],i)+jack(av1n0xmn[t1+reg_low],i)-jack(av1n0ymn[t1+reg_low],i)+jack(av1n0zmn[t1+reg_low],i)+jack(av1n0xp[t1+reg_low],i)-jack(av1n0yp[t1+reg_low],i)+jack(av1n0zp[t1+reg_low],i)+jack(av1n0xmp[t1+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)+jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)+jack(av1n0zm[t2+reg_low],i)+jack(av1n0xn[t2+reg_low],i)-jack(av1n0yn[t2+reg_low],i)+jack(av1n0zn[t2+reg_low],i)+jack(av1n0xmn[t2+reg_low],i)-jack(av1n0ymn[t2+reg_low],i)+jack(av1n0zmn[t2+reg_low],i)+jack(av1n0xp[t2+reg_low],i)-jack(av1n0yp[t2+reg_low],i)+jack(av1n0zp[t2+reg_low],i)+jack(av1n0xmp[t2+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
            elif nmom==12:
                x=x+((((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)+jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)+jack(av1n0zm[t1+reg_low],i)+jack(av1n0xn[t1+reg_low],i)-jack(av1n0yn[t1+reg_low],i)+jack(av1n0zn[t1+reg_low],i)+jack(av1n0xmn[t1+reg_low],i)-jack(av1n0ymn[t1+reg_low],i)+jack(av1n0zmn[t1+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)+jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)+jack(av1n0zm[t2+reg_low],i)+jack(av1n0xn[t2+reg_low],i)-jack(av1n0yn[t2+reg_low],i)+jack(av1n0zn[t2+reg_low],i)+jack(av1n0xmn[t2+reg_low],i)-jack(av1n0ymn[t2+reg_low],i)+jack(av1n0zmn[t2+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
            elif nmom==24:
                x=x+((((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)+jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)+jack(av1n0zm[t1+reg_low],i)+jack(av1n0xn[t1+reg_low],i)-jack(av1n0yn[t1+reg_low],i)+jack(av1n0zn[t1+reg_low],i)+jack(av1n0xmn[t1+reg_low],i)-jack(av1n0ymn[t1+reg_low],i)+jack(av1n0zmn[t1+reg_low],i)+jack(av1n0xp[t1+reg_low],i)-jack(av1n0yp[t1+reg_low],i)+jack(av1n0zp[t1+reg_low],i)+jack(av1n0xmp[t1+reg_low],i)-jack(av1n0ymp[t1+reg_low],i)+jack(av1n0zmp[t1+reg_low],i)+jack(av1n0xnp[t1+reg_low],i)-jack(av1n0ynp[t1+reg_low],i)+jack(av1n0znp[t1+reg_low],i)+jack(av1n0xmnp[t1+reg_low],i)-jack(av1n0ymnp[t1+reg_low],i)+jack(av1n0zmnp[t1+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)+jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)+jack(av1n0zm[t2+reg_low],i)+jack(av1n0xn[t2+reg_low],i)-jack(av1n0yn[t2+reg_low],i)+jack(av1n0zn[t2+reg_low],i)+jack(av1n0xmn[t2+reg_low],i)-jack(av1n0ymn[t2+reg_low],i)+jack(av1n0zmn[t2+reg_low],i)+jack(av1n0xp[t2+reg_low],i)-jack(av1n0yp[t2+reg_low],i)+jack(av1n0zp[t2+reg_low],i)+jack(av1n0xmp[t2+reg_low],i)-jack(av1n0ymp[t2+reg_low],i)+jack(av1n0zmp[t2+reg_low],i)+jack(av1n0xnp[t2+reg_low],i)-jack(av1n0ynp[t2+reg_low],i)+jack(av1n0znp[t2+reg_low],i)+jack(av1n0xmnp[t2+reg_low],i)-jack(av1n0ymnp[t2+reg_low],i)+jack(av1n0zmnp[t2+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
    '''
        covmat[t1][t2]=x
        covmat[t2][t1]=x  
        
        
def chi(a):
    return (nconf-1-reg_low-cut)/(nconf-reg_low-cut)*np.dot(np.transpose([i-a for i in avn0[reg_low:reg_up]]),np.matmul(np.linalg.inv(covmat),[i-a for i in avn0[reg_low:reg_up]]))

mbar=minimize(chi,0.1,method='Nelder-Mead', tol=1e-6)


def jackmass(t1,i):
    #+jack(av1n0xm[t1],i)-jack(av1n0ym[t1],i)+jack(av1n0zm[t1],i)
    return (((sum_with_exceptions_jack(av1n0, t1, i)))/(1/3*np.sqrt((jack(avdx[t1],i)+jack(avdy[t1],i)+jack(avdz[t1],i))*jack(avb[t1],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1))*np.exp(-mb*(30-(t1)))))*pre

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

ax.set_xlabel('time Steps')
ax.set_ylabel('A0')
ax.errorbar(list(range(dt))[1:dt], avn0[1:dt], yerr=errn0[1:dt],fmt='x', label='nsq={}'.format(nsq))
plt.axhline(y = mbar.x[0], color = 'r', linestyle = 'dashed')
plt.fill_between(list(range(dt))[reg_low:reg_up], mbar.x[0]+sigma, mbar.x[0]-sigma, color='r',alpha=0.2)
#ax.set_yscale('log')
plt.savefig('./Fits/A0-Av-nsq{}-Fit.png'.format(nsq))

df3 = pd.DataFrame(columns=['EffectiveMass','Error','RegUp','RegLow'])
df3['EffectiveMass']=mbar.x
df3['Error']=sigma  
df3['RegUp']=reg_up
df3['RegLow']=reg_low    
df3.to_csv('./Fits/A0-Av-nsq{}-Fit.csv'.format(nsq), sep='\t')






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
'''



'''
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
elif len(mom[nsq])==16:
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
    

    nmom=16
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

elif len(mom[nsq])==48:
    dset1xn0q=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][0])]
    dset1xn0bq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][0])]
    dset1yn0q=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][1])]
    dset1yn0bq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][1])]
    dset1zn0q=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][2])]
    dset1zn0bq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][2])]   
    dset1xn0mq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][3])]
    dset1xn0bmq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][3])]
    dset1yn0mq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][4])]
    dset1yn0bmq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][4])]
    dset1zn0mq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][5])]
    dset1zn0bmq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][5])]   
    
    dset1xn0nq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][6])]
    dset1xn0bnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][6])]
    dset1yn0nq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][7])]
    dset1yn0bnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][7])]
    dset1zn0nq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][8])]
    dset1zn0bnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][8])]   
    dset1xn0mnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][9])]
    dset1xn0bmnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][9])]
    dset1yn0mnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][10])]
    dset1yn0bmnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][10])]
    dset1zn0mnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][11])]
    dset1zn0bmnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][11])]   
    
    dset1xn0pq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][12])]
    dset1xn0bpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][12])]
    dset1yn0pq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][13])]
    dset1yn0bpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][13])]
    dset1zn0pq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][14])]
    dset1zn0bpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][14])]   
    dset1xn0mpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][15])]
    dset1xn0bmpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][15])]
    dset1yn0mpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][16])]
    dset1yn0bmpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][16])]
    dset1zn0mpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][17])]
    dset1zn0bmpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][17])]   
    
    dset1xn0npq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][18])]
    dset1xn0bnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][18])]
    dset1yn0npq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][19])]
    dset1yn0bnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][19])]
    dset1zn0npq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][20])]
    dset1zn0bnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][20])]   
    dset1xn0mnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][21])]
    dset1xn0bmnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][21])]
    dset1yn0mnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][22])]
    dset1yn0bmnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][22])]
    dset1zn0mnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][23])]
    dset1zn0bmnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][23])]   
    
    dset1xn0q=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][0])]
    dset1xn0bq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][0])]
    dset1yn0q=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][1])]
    dset1yn0bq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][1])]
    dset1zn0q=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][2])]
    dset1zn0bq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][2])]   
    dset1xn0mq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][3])]
    dset1xn0bmq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][3])]
    dset1yn0mq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][4])]
    dset1yn0bmq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][4])]
    dset1zn0mq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][5])]
    dset1zn0bmq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][5])]   
    
    dset1xn0nq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][6])]
    dset1xn0bnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][6])]
    dset1yn0nq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][7])]
    dset1yn0bnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][7])]
    dset1zn0nq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][8])]
    dset1zn0bnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][8])]   
    dset1xn0mnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][9])]
    dset1xn0bmnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][9])]
    dset1yn0mnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][10])]
    dset1yn0bmnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][10])]
    dset1zn0mnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][11])]
    dset1zn0bmnq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][11])]   
    
    dset1xn0pq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][12])]
    dset1xn0bpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][12])]
    dset1yn0pq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][13])]
    dset1yn0bpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][13])]
    dset1zn0pq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][14])]
    dset1zn0bpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][14])]   
    dset1xn0mpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][15])]
    dset1xn0bmpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][15])]
    dset1yn0mpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][16])]
    dset1yn0bmpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][16])]
    dset1zn0mpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][17])]
    dset1zn0bmpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][17])]   
    
    dset1xn0npq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][18])]
    dset1xn0bnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][18])]
    dset1yn0npq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][19])]
    dset1yn0bnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][19])]
    dset1zn0npq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][20])]
    dset1zn0bnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][20])]   
    dset1xn0mnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][21])]
    dset1xn0bmnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][21])]
    dset1yn0mnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][22])]
    dset1yn0bmnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][22])]
    dset1zn0mnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][23])]
    dset1zn0bmnpq=f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][23])]   
    nmom=24       


tmp1x = np.mean((np.real(dsets[0][k, :, j]) + np.real(dsetsb[0][k, :, dt-j]))) / 2
tmp1y = np.mean((np.real(dsets[1][k, :, j]) + np.real(dsetsb[1][k, :, dt-j]))) / 2
tmp1z = np.mean((np.real(dsets[2][k, :, j]) + np.real(dsetsb[2][k, :, dt-j]))) / 2
if nmom==6:
    tmp1xm = np.mean((np.real(dsets[3][k, :, j]) + np.real(dsetsb[3][k, :, dt-j]))) / 2
    tmp1ym = np.mean((np.real(dsets[4][k, :, j]) + np.real(dsetsb[4][k, :, dt-j]))) / 2
    tmp1zm = np.mean((np.real(dsets[5][k, :, j]) + np.real(dsetsb[5][k, :, dt-j]))) / 2
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
    
elif nmom==16:
    tmp1xm = np.mean((np.real(dset1xn0m[k, :, j]) + np.real(dset1xn0bm[k, :, dt-j]))) / 2
    tmp1ym = np.mean((np.real(dset1yn0m[k, :, j]) + np.real(dset1yn0bm[k, :, dt-j]))) / 2
    tmp1zm = np.mean((np.real(dset1zn0m[k, :, j]) + np.real(dset1zn0bm[k, :, dt-j]))) / 2
    
    tmp1xn = np.mean((np.real(dset1xn0n[k, :, j]) + np.real(dset1xn0bn[k, :, dt-j]))) / 2
    tmp1yn = np.mean((np.real(dset1yn0n[k, :, j]) + np.real(dset1yn0bn[k, :, dt-j]))) / 2
    tmp1zn = np.mean((np.real(dset1zn0n[k, :, j]) + np.real(dset1zn0bn[k, :, dt-j]))) / 2
    
    tmp1xmn = np.mean((np.real(dset1xn0mn[k, :, j]) + np.real(dset1xn0bmn[k, :, dt-j]))) / 2
    tmp1ymn = np.mean((np.real(dset1yn0mn[k, :, j]) + np.real(dset1yn0bmn[k, :, dt-j]))) / 2
    tmp1zmn = np.mean((np.real(dset1zn0mn[k, :, j]) + np.real(dset1zn0bmn[k, :, dt-j]))) / 2
    
    tmp1xp = np.mean((np.real(dset1xn0p[k, :, j]) + np.real(dset1xn0bp[k, :, dt-j]))) / 2
    tmp1yp = np.mean((np.real(dset1yn0p[k, :, j]) + np.real(dset1yn0bp[k, :, dt-j]))) / 2
    tmp1zp = np.mean((np.real(dset1zn0p[k, :, j]) + np.real(dset1zn0bp[k, :, dt-j]))) / 2
    
    tmp1xmp = np.mean((np.real(dset1xn0mp[k, :, j]) + np.real(dset1xn0bmp[k, :, dt-j]))) / 2
    tmp1ymp=0
    tmp1zmp=0
    
    tmp1xnp=0
    tmp1ynp=0
    tmp1znp=0
    
    tmp1xmnp=0
    tmp1ymnp=0
    tmp1zmnp=0
    
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
    tmp1xm = np.mean((np.real(dset1xn0m[k, :, j]) + np.real(dset1xn0bm[k, :, dt-j]))) / 2
    tmp1ym = np.mean((np.real(dset1yn0m[k, :, j]) + np.real(dset1yn0bm[k, :, dt-j]))) / 2
    tmp1zm = np.mean((np.real(dset1zn0m[k, :, j]) + np.real(dset1zn0bm[k, :, dt-j]))) / 2
    
    tmp1xn = np.mean((np.real(dset1xn0n[k, :, j]) + np.real(dset1xn0bn[k, :, dt-j]))) / 2
    tmp1yn = np.mean((np.real(dset1yn0n[k, :, j]) + np.real(dset1yn0bn[k, :, dt-j]))) / 2
    tmp1zn = np.mean((np.real(dset1zn0n[k, :, j]) + np.real(dset1zn0bn[k, :, dt-j]))) / 2
    tmp1xmn = np.mean((np.real(dset1xn0mn[k, :, j]) + np.real(dset1xn0bmn[k, :, dt-j]))) / 2
    tmp1ymn = np.mean((np.real(dset1yn0mn[k, :, j]) + np.real(dset1yn0bmn[k, :, dt-j]))) / 2
    tmp1zmn = np.mean((np.real(dset1zn0mn[k, :, j]) + np.real(dset1zn0bmn[k, :, dt-j]))) / 2

    tmp1xp = np.mean((np.real(dset1xn0p[k, :, j]) + np.real(dset1xn0bp[k, :, dt-j]))) / 2
    tmp1yp = np.mean((np.real(dset1yn0p[k, :, j]) + np.real(dset1yn0bp[k, :, dt-j]))) / 2
    tmp1zp = np.mean((np.real(dset1zn0p[k, :, j]) + np.real(dset1zn0bp[k, :, dt-j]))) / 2
    
    tmp1xmp = np.mean((np.real(dset1xn0mp[k, :, j]) + np.real(dset1xn0bmp[k, :, dt-j]))) / 2
    tmp1ymp = np.mean((np.real(dset1yn0mp[k, :, j]) + np.real(dset1yn0bmp[k, :, dt-j]))) / 2
    tmp1zmp = np.mean((np.real(dset1zn0mp[k, :, j]) + np.real(dset1zn0bmp[k, :, dt-j]))) / 2
    
    tmp1xnp = np.mean((np.real(dset1xn0np[k, :, j]) + np.real(dset1xn0bnp[k, :, dt-j]))) / 2
    tmp1ynp = np.mean((np.real(dset1yn0np[k, :, j]) + np.real(dset1yn0bnp[k, :, dt-j]))) / 2
    tmp1znp = np.mean((np.real(dset1zn0np[k, :, j]) + np.real(dset1zn0bnp[k, :, dt-j]))) / 2
    tmp1xmnp = np.mean((np.real(dset1xn0mnp[k, :, j]) + np.real(dset1xn0bmnp[k, :, dt-j]))) / 2
    tmp1ymnp = np.mean((np.real(dset1yn0mnp[k, :, j]) + np.real(dset1yn0bmnp[k, :, dt-j]))) / 2
    tmp1zmnp = np.mean((np.real(dset1zn0mnp[k, :, j]) + np.real(dset1zn0bmnp[k, :, dt-j]))) / 2

 
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
