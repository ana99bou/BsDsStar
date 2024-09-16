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

def sum_with_exceptions(lst,nsq):
    total = 0
    if nsq==5:
        
        for i in range(0, len(lst)//2, 6):
            total += lst[i]+lst[i+1]+lst[i+2]
            total -= lst[i+3]+lst[i+4]+lst[i+5]
        
        for i in range(len(lst)//2,len(lst), 6):
            total += 1/2*(lst[i]+lst[i+1]+lst[i+2])
            total -= 1/2*(lst[i+3]+lst[i+4]+lst[i+5])
        
        return total/len(lst)
    else:
        for i in range(0, len(lst), 6):
            total += lst[i]+lst[i+1]+lst[i+2]
            total -= lst[i+3]+lst[i+4]+lst[i+5]
        return total/len(lst)


def sum_with_exceptions_jack(lst,nsq,j,i):
    total = 0
    if nsq==5:
        for k in range(0, len(lst)//2, 6):
            total += jack(lst[k][j],i)+jack(lst[k+1][j],i)+jack(lst[k+2][j],i)
            total -= jack(lst[k+3][j],i)+jack(lst[k+4][j],i)+jack(lst[k+5][j],i)
        for k in range(len(lst)//2,len(lst), 6):
            total += 1/2*(jack(lst[k][j],i)+jack(lst[k+1][j],i)+jack(lst[k+2][j],i))
            total -= 1/2*(jack(lst[k+3][j],i)+jack(lst[k+4][j],i)+jack(lst[k+5][j],i))
    else:    
        for k in range(0, len(lst), 6):
            total += jack(lst[k][j],i)+jack(lst[k+1][j],i)+jack(lst[k+2][j],i)
            total -= jack(lst[k+3][j],i)+jack(lst[k+4][j],i)+jack(lst[k+5][j],i)
    '''
    for i, num in enumerate(lst):
        if i % 3 == 2:
            total -= jack(num[j],i)  # Subtract every third element starting from the second
        else:
            total += jack(num[j],i)  # Add other elements
    '''
    return total/len(lst)

'''
def sum_with_exceptions_jack(lst,nsq,j,i):
    total = 0
    if nsq==5:
        n = len(lst)
        half = n // 2
        
        # Process first half of the list (no factor)
        for k in range(0, half, 6):
        #for k, num in enumerate(lst):
            total += (jack(lst[k][j],i)+jack(lst[k+1][j],i)+jack(lst[k+2][j],i))/3
            total -= (jack(lst[k+3][j],i)+jack(lst[k+4][j],i)+jack(lst[k+4][j],i))/3
        
        # Process second half of the list (with factor of 1/2)
        for k in range(half, n, 6):
            total += 1/2*(jack(lst[k][j],i)+jack(lst[k+1][j],i)+jack(lst[k+2][j],i))/3
            total -= 1/2*(jack(lst[k+3][j],i)+jack(lst[k+4][j],i)+jack(lst[k+4][j],i))/3
        
        return total/len(lst)
    else:
        total = 0
        for k in range(0,len(lst),6):
        #for k, num in enumerate(lst):
            total += (jack(lst[k][j],i)+jack(lst[k+1][j],i)+jack(lst[k+2][j],i))/3
            total -= (jack(lst[k+3][j],i)+jack(lst[k+4][j],i)+jack(lst[k+4][j],i))/3
        
        return total
'''


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

dsets=[f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/forward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]
dsetsb=[f["/CHARM_PT_SEQ_SM12.14_s0.02144/c0.248/dT30/{}/backward/data".format(mom[nsq][i])] for i in range(len(mom[nsq]))]

nmom=len(mom[nsq])


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
        
        tmp=[np.mean((np.imag(dsets[i][k, :, j]) - np.imag(dsetsb[i][k, :, dt-j]))) / 2 for i in range(nmom)]
        
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
    #avn0[j] = pre * (((temp1x + temp1y + temp1z-temp1xm - temp1ym - temp1zm +temp1xn +temp1yn + temp1zn-temp1xmn - temp1ymn - temp1zmn+1/2*(temp1xp + temp1yp + temp1zp-temp1xmp - temp1ymp - temp1zmp +temp1xnp +temp1ynp + temp1znp-temp1xmnp - temp1ymnp - temp1zmnp)) / nmom) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    avn0[j] = pre * (sum_with_exceptions(temp, nsq) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    #avn0[j] = pre * ((temp[0]+temp[1]+temp[2]-temp[3]-temp[4]-temp[5])/6 / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    
    #+temp1xn + temp1yn + temp1zn+temp1xmn + temp1ymn + temp1zmn

errn0=np.zeros(shape=(dt))
#errnn0=np.zeros(shape=(96))

for j in range(dt):
    x=0
    for i in range(nconf):
        #x=x+((((jack(av1n0x[j],i)+jack(av1n0y[j],i)+jack(av1n0z[j],i)-jack(av1n0xm[j],i)-jack(av1n0ym[j],i)-jack(av1n0zm[j],i)+jack(av1n0xn[j],i)+jack(av1n0yn[j],i)+jack(av1n0zn[j],i)-jack(av1n0xmn[j],i)-jack(av1n0ymn[j],i)-jack(av1n0zmn[j],i)+1/2*(jack(av1n0xp[j],i)+jack(av1n0yp[j],i)+jack(av1n0zp[j],i)-jack(av1n0xmp[j],i)-jack(av1n0ymp[j],i)-jack(av1n0zmp[j],i)+jack(av1n0xnp[j],i)+jack(av1n0ynp[j],i)+jack(av1n0znp[j],i)-jack(av1n0xmnp[j],i)-jack(av1n0ymnp[j],i)-jack(av1n0zmnp[j],i)))/nmom)/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(30-j))))*pre-avn0[j])**2
        #x=x+(((jack(av1n0[0][j],i)+jack(av1n0[1][j],i)+jack(av1n0[2][j],i)-jack(av1n0[3][j],i)-jack(av1n0[4][j],i)-jack(av1n0[5][j],i))/6/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(dt-j))))*pre-avn0[j])**2
        x=x+((sum_with_exceptions_jack(av1n0,nsq,j,i)/(1/3*np.sqrt((jack(avdx[j],i)+jack(avdy[j],i)+jack(avdz[j],i))*jack(avb[j],i))))*np.sqrt((4*md*mb)/(np.exp(-md*j)*np.exp(-mb*(dt-j))))*pre-avn0[j])**2
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
#plt.savefig('./V/V-nsq{}.png'.format(nsq))

#np.savetxt('./V/V-nsq{}.txt'.format(nsq), np.c_[np.absolute(avn0), errn0])

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
            x=x+(((sum_with_exceptions_jack(av1n0,nsq, t1+reg_low, i))/(1/3*np.sqrt((jack(avdx[t1+reg_low],i)+jack(avdy[t1+reg_low],i)+jack(avdz[t1+reg_low],i))*jack(avb[t1+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1+reg_low))*np.exp(-mb*(dt-(t1+reg_low)))))*pre-avn0[t1])*(((sum_with_exceptions_jack(av1n0,nsq, t2+reg_low, i))/(1/3*np.sqrt((jack(avdx[t2+reg_low],i)+jack(avdy[t2+reg_low],i)+jack(avdz[t2+reg_low],i))*jack(avb[t2+reg_low],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t2+reg_low))*np.exp(-mb*(dt-(t2+reg_low)))))*pre-avn0[t2])            
            
            
        covmat[t1][t2]=x
        covmat[t2][t1]=x  
        
        
def chi(a):
    return (nconf-1-reg_low-cut)/(nconf-reg_low-cut)*np.dot(np.transpose([i-a for i in avn0[reg_low:reg_up]]),np.matmul(np.linalg.inv(covmat),[i-a for i in avn0[reg_low:reg_up]]))

mbar=minimize(chi,0.1,method='Nelder-Mead', tol=1e-6)


def jackmass(t1,i):
    #+jack(av1n0xm[t1],i)-jack(av1n0ym[t1],i)+jack(av1n0zm[t1],i)
    return (((sum_with_exceptions_jack(av1n0, nsq,t1, i)))/(1/3*np.sqrt((jack(avdx[t1],i)+jack(avdy[t1],i)+jack(avdz[t1],i))*jack(avb[t1],i))))*np.sqrt((4*md*mb)/(np.exp(-md*(t1))*np.exp(-mb*(30-(t1)))))*pre

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

'''
