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
nsq=5
##########


f = h5py.File("BsDsStar_M2.h5", "r")
f2 = h5py.File("./BsDsStar_2ptData_M2.h5", "r")

#f["/CHARM_PT_SEQ_SM10.36_s0.025/c0.340/dT26/final_state_GX/operator_GammaXGamma5/n2_1/1_0_0/forward/data"]

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


dsets=[f["/CHARM_PT_SEQ_SM10.36_s0.025/c0.340/dT26/{}/forward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]
dsetsb=[f["/CHARM_PT_SEQ_SM10.36_s0.025/c0.340/dT26/{}/backward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]

nmom=len(mom[nsq])

#/cl_SM10.36_SM10.36_0.025/c0.340/operator_GammaX/n2_4/data
dsxn0=f2["/cl_SM10.36_SM10.36_0.025/c0.340/operator_GammaX/n2_{}/data".format(nsq)]
dsyn0=f2["/cl_SM10.36_SM10.36_0.025/c0.340/operator_GammaY/n2_{}/data".format(nsq)]
dszn0=f2["/cl_SM10.36_SM10.36_0.025/c0.340/operator_GammaZ/n2_{}/data".format(nsq)]
bsn0=f2["/hl_SM10.36_SM10.36_0.025_m3.49_csw3.07_zeta1.76/operator_Gamma5/n2_0/data"]



#####Eff Masses
mdlist=[0.9122552490234405,0.9363156127929717,0.9604193115234406,0.9817083740234408,1.0055334472656285,1.0299121093750034]
md=mdlist[nsq]

# Constants

mb = 2.2526141357421965
md0 = 0.9122552490234405
pre = md0 / (2 * md * mb)

#nconf=100
#nconf=889
nconf=800
dt=26
ts=64

findx = np.zeros(dt)
findy = np.zeros(dt)
findz = np.zeros(dt)
finb = np.zeros(dt)
fin3pt=np.zeros((dt,nmom))

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
        tmpb = np.mean(np.real(bsn0[k, :, j+1]) + np.real(bsn0[k, :, ts-1-j])) / 2 if j != 0 else np.mean(np.real(bsn0[k, :, 0]))

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
    findx[j]=tempdx
    findy[j]=tempdy
    findz[j]=tempdz
    finb[j]=tempb
    fin3pt[j]=temp
    #avn0[j] = pre * (((temp1x - temp1y + temp1z+temp1xm - temp1ym + temp1zm+temp1xn - temp1yn + temp1zn+temp1xmn - temp1ymn + temp1zmn+temp1xp - temp1yp + temp1zp+temp1xmp - temp1ymp + temp1zmp+temp1xnp - temp1ynp + temp1znp+temp1xmnp - temp1ymnp + temp1zmnp) / nmom) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    #avn0[j] = pre * (((sum_with_exceptions(temp))) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    #avn0[j] = pre * (((temp[0]+temp[1]+temp[2]+temp[3]+temp[4]+temp[5]) / nmom) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))

for j in range(dt):
    avn0[j] = pre * (((sum_with_exceptions(fin3pt[j,:]))) / (1/3 * np.sqrt((findx[j] + findy[j] + findz[j]) * finb[dt-1-j]))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
errn0=np.zeros(shape=(dt))
#errnn0=np.zeros(shape=(96))

for j in range(dt):
    x=0
    for i in range(nconf):
#        x=x+((((jack(av1n0x[j],i)-jack(av1n0y[j],i)+jack(av1n0z[j],i)+jack(av1n0xm[j],i)-jack(av1n0ym[j],i)+jack(av1n0zm[j],i)+jack(av1n0xn[j],i)-jack(av1n0yn[j],i)+jack(av1n0zn[j],i)+jack(av1n0xmn[j],i)-jack(av1n0ymn[j],i)+jack(av1n0zmn[j],i)+jack(av1n0xp[j],i)-jack(av1n0yp[j],i)+jack(av1n0zp[j],i)+jack(av1n0xmp[j],i)-jack(av1n0ymp[j],i)+jack(av1n0zmp[j],i)+jack(av1n0xnp[j],i)-jack(av1n0ynp[j],i)+jack(av1n0znp[j],i)+jack(av1n0xmnp[j],i)-jack(av1n0ymnp[j],i)+jack(av1n0zmnp[j],i))/nmom)/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(30-j))))*pre-avn0[j])**2
        #x=x+((((sum_with_exceptions_jack(av1n0, j, i))/nmom)/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(30-j))))*pre-avn0[j])**2
        #x=x+((((jack(av1n0[0][j],i)+jack(av1n0[1][j],i)+jack(av1n0[2][j],i)+jack(av1n0[3][j],i)+jack(av1n0[4][j],i)+jack(av1n0[5][j],i))/nmom)/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(30-j))))*pre-avn0[j])**2
        x=x+((((sum_with_exceptions_jack(av1n0, j, i)))/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[dt-1-j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(dt-j))))*pre-avn0[j])**2

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
'''
reg_low=20
reg_up=22

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
            x=x+(((sum_with_exceptions_jack(av1n0, t1+reg_low, i))/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[dt-1-(t1+reg_low)],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*(((sum_with_exceptions_jack(av1n0, t2+reg_low, i))/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[dt-1-(t2+reg_low)],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            

        covmat[t1][t2]=x
        covmat[t2][t1]=x  
        
        
def chi(a):
    return (nconf-1-reg_low-cut)/(nconf-reg_low-cut)*np.dot(np.transpose([i-a for i in avn0[reg_low:reg_up]]),np.matmul(np.linalg.inv(covmat),[i-a for i in avn0[reg_low:reg_up]]))

mbar=minimize(chi,0.1,method='Nelder-Mead', tol=1e-6)


def jackmass(t1,i):
    #+jack(av1n0xm[t1],i)-jack(av1n0ym[t1],i)+jack(av1n0zm[t1],i)
    return (((sum_with_exceptions_jack(av1n0, t1, i)))/(1/3*np.sqrt((jack(avdx[t1],i)+jack(avdy[t1],i)+jack(avdz[t1],i))*jack(avb[dt-1-t1],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1))*np.exp(-mb*(dt-(t1)))))*pre

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

