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
nsq=1
##########


f = h5py.File("BsDsStar_M2.h5", "r")
f2 = h5py.File("./BsDsStar_2ptData_M2.h5", "r")

mom=[['final_state_GX/operator_GammaXGamma5/n2_0/0_0_0','final_state_GY/operator_GammaYGamma5/n2_0/0_0_0','final_state_GZ/operator_GammaZGamma5/n2_0/0_0_0'],
     ['final_state_GX/operator_GammaXGamma5/n2_1/0_1_0','final_state_GY/operator_GammaYGamma5/n2_1/1_0_0', 'final_state_GZ/operator_GammaZGamma5/n2_1/1_0_0','final_state_GX/operator_GammaXGamma5/n2_1/0_0_1','final_state_GY/operator_GammaYGamma5/n2_1/0_0_1', 'final_state_GZ/operator_GammaZGamma5/n2_1/0_1_0',
      'final_state_GX/operator_GammaXGamma5/n2_1/0_-1_0','final_state_GY/operator_GammaYGamma5/n2_1/-1_0_0', 'final_state_GZ/operator_GammaZGamma5/n2_1/-1_0_0','final_state_GX/operator_GammaXGamma5/n2_1/0_0_-1','final_state_GY/operator_GammaYGamma5/n2_1/0_0_-1', 'final_state_GZ/operator_GammaZGamma5/n2_1/0_-1_0'],
     ['final_state_GX/operator_GammaXGamma5/n2_2/0_1_1','final_state_GY/operator_GammaYGamma5/n2_2/1_0_1','final_state_GZ/operator_GammaZGamma5/n2_2/1_1_0',
      'final_state_GX/operator_GammaXGamma5/n2_2/0_-1_1','final_state_GY/operator_GammaYGamma5/n2_2/-1_0_1','final_state_GZ/operator_GammaZGamma5/n2_2/-1_1_0',
      'final_state_GX/operator_GammaXGamma5/n2_2/0_1_-1','final_state_GY/operator_GammaYGamma5/n2_2/1_0_-1','final_state_GZ/operator_GammaZGamma5/n2_2/1_-1_0',
      'final_state_GX/operator_GammaXGamma5/n2_2/0_-1_-1','final_state_GY/operator_GammaYGamma5/n2_2/-1_0_-1','final_state_GZ/operator_GammaZGamma5/n2_2/-1_-1_0'],
     [''],
     ['final_state_GX/operator_GammaXGamma5/n2_4/0_2_0','final_state_GY/operator_GammaYGamma5/n2_4/2_0_0', 'final_state_GZ/operator_GammaZGamma5/n2_4/2_0_0','final_state_GX/operator_GammaXGamma5/n2_4/0_0_2','final_state_GY/operator_GammaYGamma5/n2_4/0_0_2', 'final_state_GZ/operator_GammaZGamma5/n2_4/0_2_0',
      'final_state_GX/operator_GammaXGamma5/n2_4/0_-2_0','final_state_GY/operator_GammaYGamma5/n2_4/-2_0_0', 'final_state_GZ/operator_GammaZGamma5/n2_4/-2_0_0','final_state_GX/operator_GammaXGamma5/n2_4/0_0_-2','final_state_GY/operator_GammaYGamma5/n2_4/0_0_-2', 'final_state_GZ/operator_GammaZGamma5/n2_4/0_-2_0'],
     ['final_state_GX/operator_GammaXGamma5/n2_5/0_2_1','final_state_GY/operator_GammaYGamma5/n2_5/2_0_1','final_state_GZ/operator_GammaZGamma5/n2_5/2_1_0','final_state_GX/operator_GammaXGamma5/n2_5/0_1_2','final_state_GY/operator_GammaYGamma5/n2_5/1_0_2','final_state_GZ/operator_GammaZGamma5/n2_5/1_2_0',
      'final_state_GX/operator_GammaXGamma5/n2_5/0_-2_1','final_state_GY/operator_GammaYGamma5/n2_5/-2_0_1','final_state_GZ/operator_GammaZGamma5/n2_5/-2_1_0','final_state_GX/operator_GammaXGamma5/n2_5/0_-1_2','final_state_GY/operator_GammaYGamma5/n2_5/-1_0_2','final_state_GZ/operator_GammaZGamma5/n2_5/-1_2_0',
      'final_state_GX/operator_GammaXGamma5/n2_5/0_2_-1','final_state_GY/operator_GammaYGamma5/n2_5/2_0_-1','final_state_GZ/operator_GammaZGamma5/n2_5/2_-1_0','final_state_GX/operator_GammaXGamma5/n2_5/0_1_-2','final_state_GY/operator_GammaYGamma5/n2_5/1_0_-2','final_state_GZ/operator_GammaZGamma5/n2_5/1_-2_0']]


ptmom=['0_0_0','1_0_0','1_1_0','1_1_1', '2_0_0','2_1_0']

dsets=[f["/CHARM_PT_SEQ_SM10.36_s0.025/c0.340/dT26/{}/forward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]
dsetsb=[f["/CHARM_PT_SEQ_SM10.36_s0.025/c0.340/dT26/{}/backward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]

nmom=len(mom[nsq])

dsxn0=f2["/cl_SM10.36_SM10.36_0.025/c0.340/operator_GammaX/n2_{}/data".format(nsq)]
dsyn0=f2["/cl_SM10.36_SM10.36_0.025/c0.340/operator_GammaX/n2_{}/data".format(nsq)]
dszn0=f2["/cl_SM10.36_SM10.36_0.025/c0.340/operator_GammaX/n2_{}/data".format(nsq)]
bsn0=f2["/hl_SM10.36_SM10.36_0.025_m3.49_csw3.07_zeta1.76/operator_Gamma5/n2_0/data"]



#####Eff Masses
mdlist=[0.9122552490234405,0.9363156127929717,0.9604193115234406,0.9817083740234408,1.0055334472656285,1.0299121093750034]
md=mdlist[nsq]

# Constants

mb = 2.2526141357421965
md0 = 0.9122552490234405
pre = -1/(mb+md)

nconf=889
dt=26
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
        tmpb = np.mean(np.real(bsn0[k, :, dt-1-j]) + np.real(bsn0[k, :, ts-dt+1+j])) / 2 if j != dt-1 else np.mean(np.real(bsn0[k, :, 0]))

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
    #avn0[j] = pre * (((temp1x - temp1y + temp1z+temp1xm - temp1ym + temp1zm+temp1xn - temp1yn + temp1zn+temp1xmn - temp1ymn + temp1zmn) / nmom) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    avn0[j] = pre * (sum_with_exceptions(temp) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    

errn0=np.zeros(shape=(dt))
#errnn0=np.zeros(shape=(96))

for j in range(dt):
    x=0
    for i in range(nconf):
        x=x+((sum_with_exceptions_jack(av1n0, j, i)/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(dt-j))))*pre-avn0[j])**2
    errn0[j]=np.sqrt((nconf-1)/nconf*x) 
    #errnn0[j]=np.sqrt((98-1)/98*x)*10**(85)

avn0[np.isnan(avn0)] = 0
errn0[np.isnan(errn0)] = 0

#plt.plot(list(range(96)), np.absolute(avn0)[0:96],'ro')
plt.xlabel('tsme')
plt.ylabel(r'$\widetilde{A}_1$')
plt.errorbar(list(range(dt)), np.absolute(avn0)[0:dt], yerr=errn0[0:dt],ls='none',fmt='x',label='nsq={}'.format(nsq))
#plt.yscale('log')
plt.legend()
plt.savefig('./A1/A1-nsq{}.png'.format(nsq))

np.savetxt('./A1/A1-nsq{}.txt'.format(nsq), np.c_[np.absolute(avn0), errn0])

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
            x=x+(((sum_with_exceptions_jack(av1n0, t1+reg_low, i))/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*(((sum_with_exceptions_jack(av1n0, t2+reg_low, i))/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
            '''
            if nmom==3:
                x=x+((((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i))/3)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(30-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i))/3)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(30-(t2+reg_low)))))*pre-avn0[t2])            
            elif nmom==6:
                x=x+((((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)+jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)+jack(av1n0zm[t1+reg_low],i))/6)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(30-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)+jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)+jack(av1n0zm[t2+reg_low],i))/6)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(30-(t2+reg_low)))))*pre-avn0[t2])            
            elif nmom==12:
                x=x+((((jack(av1n0x[t1+reg_low],i)-jack(av1n0y[t1+reg_low],i)+jack(av1n0z[t1+reg_low],i)+jack(av1n0xm[t1+reg_low],i)-jack(av1n0ym[t1+reg_low],i)+jack(av1n0zm[t1+reg_low],i)+jack(av1n0xn[t1+reg_low],i)-jack(av1n0yn[t1+reg_low],i)+jack(av1n0zn[t1+reg_low],i)+jack(av1n0xmn[t1+reg_low],i)-jack(av1n0ymn[t1+reg_low],i)+jack(av1n0zmn[t1+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(30-(t1+reg_low)))))*pre-avn0[t1])*((((jack(av1n0x[t2+reg_low],i)-jack(av1n0y[t2+reg_low],i)+jack(av1n0z[t2+reg_low],i)+jack(av1n0xm[t2+reg_low],i)-jack(av1n0ym[t2+reg_low],i)+jack(av1n0zm[t2+reg_low],i)+jack(av1n0xn[t2+reg_low],i)-jack(av1n0yn[t2+reg_low],i)+jack(av1n0zn[t2+reg_low],i)+jack(av1n0xmn[t2+reg_low],i)-jack(av1n0ymn[t2+reg_low],i)+jack(av1n0zmn[t2+reg_low],i))/nmom)/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
            '''
        covmat[t1][t2]=x
        covmat[t2][t1]=x  
        
        
def chi(a):
    return (nconf-1-reg_low-cut)/(nconf-reg_low-cut)*np.dot(np.transpose([i-a for i in avn0[reg_low:reg_up]]),np.matmul(np.linalg.inv(covmat),[i-a for i in avn0[reg_low:reg_up]]))

mbar=minimize(chi,0.1,method='Nelder-Mead', tol=1e-6)


def jackmass(t1,i):
    #+jack(av1n0xm[t1],i)-jack(av1n0ym[t1],i)+jack(av1n0zm[t1],i)
    #return (((jack(av1n0x[t1],i)-jack(av1n0y[t1],i)+jack(av1n0z[t1],i)+jack(av1n0xm[t1],i)-jack(av1n0ym[t1],i)+jack(av1n0zm[t1],i)+jack(av1n0xn[t1],i)-jack(av1n0yn[t1],i)+jack(av1n0zn[t1],i)+jack(av1n0xmn[t1],i)-jack(av1n0ymn[t1],i)+jack(av1n0zmn[t1],i))/nmom)/(1/3*np.sqrt((jack(avdx[t1],i)+jack(avdy[t1],i)+jack(avdz[t1],i))*jack(avb[t1],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1))*np.exp(-mb*(30-(t1)))))*pre
    return (((sum_with_exceptions_jack(av1n0, t1, i)))/(1/3*np.sqrt((jack(avdx[t1],i)+jack(avdy[t1],i)+jack(avdz[t1],i))*jack(avb[t1],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1))*np.exp(-mb*(dt-(t1)))))*pre

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

ax.set_xlabel('tsme Steps')
ax.set_ylabel('A1')
ax.errorbar(list(range(dt))[1:dt], avn0[1:dt], yerr=errn0[1:dt],fmt='x', label='nsq={}'.format(nsq))
plt.axhline(y = mbar.x[0], color = 'r', linestyle = 'dashed')
plt.fill_between(list(range(dt))[reg_low:reg_up], mbar.x[0]+sigma, mbar.x[0]-sigma, color='r',alpha=0.2)
#ax.set_yscale('log')
plt.savefig('./Fits/A1-Av-nsq{}-Fit.png'.format(nsq))

df3 = pd.DataFrame(columns=['EffectiveMass','Error','RegUp','RegLow'])
df3['EffectiveMass']=mbar.x
df3['Error']=sigma  
df3['RegUp']=reg_up
df3['RegLow']=reg_low    
df3.to_csv('./Fits/A1-Av-nsq{}-Fit.csv'.format(nsq), sep='\t')






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



        tmp1x = np.mean((np.real(dset1xn0[k, :, j]) + np.real(dset1xn0b[k, :, dt-j]))) / 2
        tmp1y = np.mean((np.real(dset1yn0[k, :, j]) + np.real(dset1yn0b[k, :, dt-j]))) / 2
        tmp1z = np.mean((np.real(dset1zn0[k, :, j]) + np.real(dset1zn0b[k, :, dt-j]))) / 2
        if nmom==6:
            tmp1xm = np.mean((np.real(dset1xn0m[k, :, j]) + np.real(dset1xn0bm[k, :, dt-j]))) / 2
            tmp1ym = np.mean((np.real(dset1yn0m[k, :, j]) + np.real(dset1yn0bm[k, :, dt-j]))) / 2
            tmp1zm = np.mean((np.real(dset1zn0m[k, :, j]) + np.real(dset1zn0bm[k, :, dt-j]))) / 2
            tmp1xn=0
            tmp1yn=0
            tmp1zn=0
            
            tmp1xmn=0
            tmp1ymn=0
            tmp1zmn=0
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
            
'''