import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy.optimize import minimize
import scipy
from sympy.utilities.iterables import multiset_permutations


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

def sum_with_exceptions(lst,nsq):
    total = 0
    if nsq==5:
        
        for i in range(0, len(lst)//2, 6):
            total += lst[i]+lst[i+1]+lst[i+2]
            total -= lst[i+3]+lst[i+4]+lst[i+5]
            print('Plus'+str(i)+str(i+1)+str(i+2))
            print('Minus'+str(i+3)+str(i+4)+str(i+5))
        
        for i in range(len(lst)//2,len(lst), 6):
            total += 1/2*(lst[i]+lst[i+1]+lst[i+2])
            total -= 1/2*(lst[i+3]+lst[i+4]+lst[i+5])
            print('Plus 1/2'+str(i)+str(i+1)+str(i+2))
            print('Minus 1/2'+str(i+3)+str(i+4)+str(i+5))
        return total/len(lst)
    else:
        for i in range(0, len(lst), 6):
            total += lst[i]+lst[i+1]+lst[i+2]
            total -= lst[i+3]+lst[i+4]+lst[i+5]
        print(total/len(lst))
        return total/len(lst)
    
    '''
    elif nsq==3:
        for i in range(len(lst)):
            total += lst[i]
        return total/len(lst)
    '''


def sum_with_exceptions_jack(lst,nsq,j,i):
    total = 0
    if nsq==5:
        for k in range(0, len(lst)//2, 6):
            total += jack(lst[k][j],i)+jack(lst[k+1][j],i)+jack(lst[k+2][j],i)
            total -= jack(lst[k+3][j],i)+jack(lst[k+4][j],i)+jack(lst[k+5][j],i)
            print('Plus'+k)
            print('Minus'+k+3)
        for k in range(len(lst)//2,len(lst), 6):
            total += 1/2*(jack(lst[k][j],i)+jack(lst[k+1][j],i)+jack(lst[k+2][j],i))
            total -= 1/2*(jack(lst[k+3][j],i)+jack(lst[k+4][j],i)+jack(lst[k+5][j],i))
            print('Plus 1/2'+k)
            print('Minus 1/2'+k)
    else:    
        for k in range(0, len(lst), 6):
            total += jack(lst[k][j],i)+jack(lst[k+1][j],i)+jack(lst[k+2][j],i)
            total -= jack(lst[k+3][j],i)+jack(lst[k+4][j],i)+jack(lst[k+5][j],i)
            print('Plus')
            print('Minus')
    return total/len(lst)

'''
elif nsq==3:
    for k in range(len(lst)):
        total += jack(lst[k][j],i)
'''

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

f = h5py.File("BsDsStar.h5", "r")


nSqs=[0,1,2,3,4,5]
#nSq=5
cycl=['X','Y','Z']
momss=[[[0,0,0]],[[1,0,0],[-1,0,0]],[[1,1,0],[1,-1,0],[-1,-1,0]],[[1,1,1],[-1,1,1],[-1,-1,1],[-1,-1,-1]],[[2,0,0],[-2,0,0]],[[2,1,0],[2,-1,0],[-2,1,0],[-2,-1,0]]]

mom=[]

for nSq in nSqs:
    iSq=nSqs.index(nSq)
    moms=momss[iSq]
    temp=[]
    #print('[')
    for gPol in cycl:
        iV=cycl.index(gPol)
        for gJ in cycl:
            iJ=cycl.index(gJ)
            iMom=3-iV-iJ #this component must not be 0
            if(gPol==gJ): continue
            string='final_state_G'+gPol+'/operator_Gamma'+gJ+'/n2_'+str(nSq)+'/'
            #print(string)
            for m in moms:
                for p in multiset_permutations(m):
                    if(p[iMom]==0): continue
                    strFinal=string+str(p[0])+'_'+str(p[1])+'_'+str(p[2])
                    #print(strFinal)
                    temp.append(strFinal)
    #print('],')
    print(temp)
    mom.append(temp)
print(mom)

#########decide here which nsq
nsq=5
##########

'''
mom=[[],
     ['final_state_GX/operator_GammaY/n2_1/0_0_1','final_state_GY/operator_GammaZ/n2_1/1_0_0', 'final_state_GZ/operator_GammaX/n2_1/0_1_0','final_state_GY/operator_GammaX/n2_1/0_0_1','final_state_GZ/operator_GammaY/n2_1/1_0_0', 'final_state_GX/operator_GammaZ/n2_1/0_1_0',
      'final_state_GY/operator_GammaX/n2_1/0_0_-1','final_state_GZ/operator_GammaY/n2_1/-1_0_0', 'final_state_GX/operator_GammaZ/n2_1/0_-1_0','final_state_GX/operator_GammaY/n2_1/0_0_-1','final_state_GY/operator_GammaZ/n2_1/-1_0_0', 'final_state_GZ/operator_GammaX/n2_1/0_-1_0'],
     ['final_state_GY/operator_GammaZ/n2_2/1_1_0','final_state_GX/operator_GammaY/n2_2/0_1_1','final_state_GZ/operator_GammaX/n2_2/1_1_0','final_state_GX/operator_GammaZ/n2_2/1_1_0','final_state_GZ/operator_GammaY/n2_2/1_1_0','final_state_GY/operator_GammaX/n2_2/1_0_1','final_state_GY/operator_GammaZ/n2_2/1_0_1','final_state_GX/operator_GammaY/n2_2/1_0_1','final_state_GZ/operator_GammaX/n2_2/0_1_1','final_state_GX/operator_GammaZ/n2_2/0_1_1','final_state_GZ/operator_GammaY/n2_2/1_0_1','final_state_GY/operator_GammaX/n2_2/0_1_1',
      'final_state_GY/operator_GammaZ/n2_2/1_-1_0','final_state_GX/operator_GammaY/n2_2/0_-1_1','final_state_GZ/operator_GammaX/n2_2/-1_1_0','final_state_GX/operator_GammaZ/n2_2/-1_1_0','final_state_GZ/operator_GammaY/n2_2/1_-1_0','final_state_GY/operator_GammaX/n2_2/-1_0_1','final_state_GY/operator_GammaZ/n2_2/1_0_-1','final_state_GX/operator_GammaY/n2_2/-1_0_1','final_state_GZ/operator_GammaX/n2_2/0_1_-1','final_state_GX/operator_GammaZ/n2_2/0_1_-1','final_state_GZ/operator_GammaY/n2_2/1_0_-1','final_state_GY/operator_GammaX/n2_2/0_-1_1',
      'final_state_GX/operator_GammaZ/n2_2/1_-1_0','final_state_GZ/operator_GammaY/n2_2/-1_1_0','final_state_GY/operator_GammaX/n2_2/1_0_-1','final_state_GY/operator_GammaZ/n2_2/-1_1_0','final_state_GX/operator_GammaY/n2_2/0_1_-1','final_state_GZ/operator_GammaX/n2_2/1_-1_0','final_state_GX/operator_GammaZ/n2_2/0_-1_1','final_state_GZ/operator_GammaY/n2_2/-1_0_1','final_state_GY/operator_GammaX/n2_2/0_1_-1','final_state_GY/operator_GammaZ/n2_2/-1_0_1','final_state_GX/operator_GammaY/n2_2/1_0_-1','final_state_GZ/operator_GammaX/n2_2/0_-1_1',
      'final_state_GX/operator_GammaZ/n2_2/-1_-1_0','final_state_GZ/operator_GammaY/n2_2/-1_-1_0','final_state_GY/operator_GammaX/n2_2/-1_0_-1','final_state_GY/operator_GammaZ/n2_2/-1_-1_0','final_state_GX/operator_GammaY/n2_2/0_-1_-1','final_state_GZ/operator_GammaX/n2_2/-1_-1_0',
      'final_state_GX/operator_GammaZ/n2_2/0_-1_-1','final_state_GZ/operator_GammaY/n2_2/-1_0_-1','final_state_GY/operator_GammaX/n2_2/0_-1_-1','final_state_GY/operator_GammaZ/n2_2/-1_0_-1','final_state_GX/operator_GammaY/n2_2/-1_0_-1','final_state_GZ/operator_GammaX/n2_2/0_-1_-1'],
     ['final_state_GX/operator_GammaY/n2_3/1_1_1','final_state_GY/operator_GammaZ/n2_3/1_1_1','final_state_GZ/operator_GammaX/n2_3/1_1_1','final_state_GY/operator_GammaX/n2_3/1_1_1','final_state_GZ/operator_GammaY/n2_3/1_1_1','final_state_GX/operator_GammaZ/n2_3/1_1_1',
      'final_state_GX/operator_GammaY/n2_3/-1_1_1','final_state_GY/operator_GammaZ/n2_3/1_1_-1','final_state_GZ/operator_GammaX/n2_3/-1_1_1', 'final_state_GY/operator_GammaX/n2_3/-1_1_1','final_state_GZ/operator_GammaY/n2_3/1_-1_1','final_state_GX/operator_GammaZ/n2_3/-1_1_1',
      'final_state_GX/operator_GammaY/n2_3/1_-1_1','final_state_GY/operator_GammaZ/n2_3/1_-1_1','final_state_GZ/operator_GammaX/n2_3/1_1_-1', 'final_state_GY/operator_GammaX/n2_3/1_-1_1','final_state_GZ/operator_GammaY/n2_3/1_1_-1','final_state_GX/operator_GammaZ/n2_3/1_1_-1',
      'final_state_GY/operator_GammaX/n2_3/1_1_-1','final_state_GZ/operator_GammaY/n2_3/-1_1_1','final_state_GX/operator_GammaZ/n2_3/1_-1_1','final_state_GX/operator_GammaY/n2_3/1_1_-1','final_state_GY/operator_GammaZ/n2_3/-1_1_1','final_state_GZ/operator_GammaX/n2_3/1_-1_1', 
      'final_state_GX/operator_GammaY/n2_3/-1_-1_1','final_state_GY/operator_GammaZ/n2_3/1_-1_-1','final_state_GZ/operator_GammaX/n2_3/-1_1_-1','final_state_GY/operator_GammaX/n2_3/-1_-1_1','final_state_GZ/operator_GammaY/n2_3/1_-1_-1','final_state_GX/operator_GammaZ/n2_3/-1_1_-1',
      'final_state_GY/operator_GammaX/n2_3/-1_1_-1','final_state_GZ/operator_GammaY/n2_3/-1_1_-1','final_state_GX/operator_GammaZ/n2_3/1_-1_-1','final_state_GX/operator_GammaY/n2_3/-1_1_-1','final_state_GY/operator_GammaZ/n2_3/-1_1_-1','final_state_GZ/operator_GammaX/n2_3/-1_-1_1',
      'final_state_GY/operator_GammaX/n2_3/1_-1_-1','final_state_GZ/operator_GammaY/n2_3/-1_-1_1','final_state_GX/operator_GammaZ/n2_3/-1_-1_1','final_state_GX/operator_GammaY/n2_3/1_-1_-1','final_state_GY/operator_GammaZ/n2_3/-1_-1_1','final_state_GZ/operator_GammaX/n2_3/1_-1_-1',
      'final_state_GY/operator_GammaX/n2_3/-1_-1_-1','final_state_GZ/operator_GammaY/n2_3/-1_-1_-1','final_state_GX/operator_GammaZ/n2_3/-1_-1_-1','final_state_GX/operator_GammaY/n2_3/-1_-1_-1','final_state_GY/operator_GammaZ/n2_3/-1_-1_-1','final_state_GZ/operator_GammaX/n2_3/-1_-1_-1'],
     ['final_state_GX/operator_GammaY/n2_4/0_0_2','final_state_GY/operator_GammaZ/n2_4/2_0_0', 'final_state_GZ/operator_GammaX/n2_4/0_2_0','final_state_GY/operator_GammaX/n2_4/0_0_2','final_state_GZ/operator_GammaY/n2_4/2_0_0', 'final_state_GX/operator_GammaZ/n2_4/0_2_0',
      'final_state_GY/operator_GammaX/n2_4/0_0_-2','final_state_GZ/operator_GammaY/n2_4/-2_0_0', 'final_state_GX/operator_GammaZ/n2_4/0_-2_0','final_state_GX/operator_GammaY/n2_4/0_0_-2','final_state_GY/operator_GammaZ/n2_4/-2_0_0', 'final_state_GZ/operator_GammaX/n2_4/0_-2_0'],
     ['final_state_GX/operator_GammaY/n2_5/0_2_1','final_state_GY/operator_GammaZ/n2_5/1_2_0','final_state_GZ/operator_GammaX/n2_5/2_1_0','final_state_GX/operator_GammaZ/n2_5/2_1_0','final_state_GY/operator_GammaX/n2_5/2_0_1','final_state_GZ/operator_GammaY/n2_5/1_2_0','final_state_GX/operator_GammaY/n2_5/2_0_1','final_state_GY/operator_GammaZ/n2_5/1_0_2','final_state_GZ/operator_GammaX/n2_5/0_1_2','final_state_GX/operator_GammaZ/n2_5/0_1_2','final_state_GY/operator_GammaX/n2_5/0_2_1','final_state_GZ/operator_GammaY/n2_5/1_0_2',
      'final_state_GX/operator_GammaY/n2_5/0_-2_1','final_state_GY/operator_GammaZ/n2_5/1_-2_0','final_state_GZ/operator_GammaX/n2_5/-2_1_0','final_state_GX/operator_GammaZ/n2_5/-2_1_0','final_state_GY/operator_GammaX/n2_5/-2_0_1','final_state_GZ/operator_GammaY/n2_5/1_-2_0','final_state_GX/operator_GammaY/n2_5/-2_0_1','final_state_GY/operator_GammaZ/n2_5/1_0_-2','final_state_GZ/operator_GammaX/n2_5/0_1_-2','final_state_GX/operator_GammaZ/n2_5/0_1_-2','final_state_GY/operator_GammaX/n2_5/0_-2_1','final_state_GZ/operator_GammaY/n2_5/1_0_-2',
      'final_state_GX/operator_GammaZ/n2_5/2_-1_0','final_state_GY/operator_GammaX/n2_5/2_0_-1','final_state_GZ/operator_GammaY/n2_5/-1_2_0','final_state_GX/operator_GammaY/n2_5/0_2_-1','final_state_GY/operator_GammaZ/n2_5/-1_2_0','final_state_GZ/operator_GammaX/n2_5/2_-1_0','final_state_GX/operator_GammaZ/n2_5/0_-1_2','final_state_GY/operator_GammaX/n2_5/0_2_-1','final_state_GZ/operator_GammaY/n2_5/-1_0_2','final_state_GX/operator_GammaY/n2_5/2_0_-1','final_state_GY/operator_GammaZ/n2_5/-1_0_2','final_state_GZ/operator_GammaX/n2_5/0_-1_2',
      'final_state_GX/operator_GammaZ/n2_5/-2_-1_0','final_state_GY/operator_GammaX/n2_5/-2_0_-1','final_state_GZ/operator_GammaY/n2_5/-1_-2_0','final_state_GX/operator_GammaY/n2_5/0_-2_-1','final_state_GY/operator_GammaZ/n2_5/-1_-2_0','final_state_GZ/operator_GammaX/n2_5/-2_-1_0','final_state_GX/operator_GammaZ/n2_5/0_-1_-2','final_state_GY/operator_GammaX/n2_5/0_-2_-1','final_state_GZ/operator_GammaY/n2_5/-1_0_-2','final_state_GX/operator_GammaY/n2_5/-2_0_-1','final_state_GY/operator_GammaZ/n2_5/-1_0_-2','final_state_GZ/operator_GammaX/n2_5/0_-1_-2',
     'final_state_GX/operator_GammaY/n2_5/0_1_2','final_state_GY/operator_GammaZ/n2_5/2_1_0','final_state_GZ/operator_GammaX/n2_5/1_2_0','final_state_GX/operator_GammaZ/n2_5/1_2_0','final_state_GY/operator_GammaX/n2_5/1_0_2','final_state_GZ/operator_GammaY/n2_5/2_1_0','final_state_GX/operator_GammaY/n2_5/1_0_2','final_state_GY/operator_GammaZ/n2_5/2_0_1','final_state_GZ/operator_GammaX/n2_5/0_2_1','final_state_GX/operator_GammaZ/n2_5/0_2_1','final_state_GY/operator_GammaX/n2_5/0_1_2','final_state_GZ/operator_GammaY/n2_5/2_0_1',
     'final_state_GX/operator_GammaY/n2_5/0_-1_2','final_state_GY/operator_GammaZ/n2_5/2_-1_0','final_state_GZ/operator_GammaX/n2_5/-1_2_0','final_state_GX/operator_GammaZ/n2_5/-1_2_0','final_state_GY/operator_GammaX/n2_5/-1_0_2','final_state_GZ/operator_GammaY/n2_5/2_-1_0','final_state_GX/operator_GammaY/n2_5/-1_0_2','final_state_GY/operator_GammaZ/n2_5/2_0_-1','final_state_GZ/operator_GammaX/n2_5/0_2_-1','final_state_GX/operator_GammaZ/n2_5/0_2_-1','final_state_GY/operator_GammaX/n2_5/0_-1_2','final_state_GZ/operator_GammaY/n2_5/2_0_-1',
     'final_state_GX/operator_GammaZ/n2_5/1_-2_0','final_state_GY/operator_GammaX/n2_5/1_0_-2','final_state_GZ/operator_GammaY/n2_5/2_1_0','final_state_GX/operator_GammaY/n2_5/0_1_-2','final_state_GY/operator_GammaZ/n2_5/-2_1_0','final_state_GZ/operator_GammaX/n2_5/1_-2_0','final_state_GX/operator_GammaZ/n2_5/0_-2_1','final_state_GY/operator_GammaX/n2_5/0_1_-2','final_state_GZ/operator_GammaY/n2_5/-2_0_1','final_state_GX/operator_GammaY/n2_5/1_0_-2','final_state_GY/operator_GammaZ/n2_5/-2_0_1','final_state_GZ/operator_GammaX/n2_5/0_-2_1',
     'final_state_GX/operator_GammaZ/n2_5/-1_-2_0','final_state_GY/operator_GammaX/n2_5/-1_0_-2','final_state_GZ/operator_GammaY/n2_5/-2_-1_0','final_state_GX/operator_GammaY/n2_5/0_-1_-2','final_state_GY/operator_GammaZ/n2_5/-2_-1_0','final_state_GZ/operator_GammaX/n2_5/-1_-2_0',
     'final_state_GX/operator_GammaZ/n2_5/0_-2_-1','final_state_GY/operator_GammaX/n2_5/0_-1_-2','final_state_GZ/operator_GammaY/n2_5/-2_0_-1','final_state_GX/operator_GammaY/n2_5/-1_0_-2','final_state_GY/operator_GammaZ/n2_5/-2_0_-1','final_state_GZ/operator_GammaX/n2_5/0_-2_-1']]
'''

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


iNs=0

for nsqStr in mom:
    for threePt in nsqStr:
        a,b,c=parse_string(threePt)
        res=eps_over_k(a,b,c)
        print('n^2=',iNs, 'iV, iJ, iMom=', a,b,c, 'eps/k=', res)
    iNs+=1     

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
#dt=1
ts=96



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
    
        tmp1=[np.mean((np.imag(dsets[i][k, :, j])))  for i in range(nmom)]
        tmp1back=[np.mean((np.imag(dsetsb[i][k, :, dt-j])))  for i in range(nmom)]
    

        for l in range(nmom):
            av1n01[l][j,k]=tmp1[l]
            av1n01back[l][j,k]=tmp1back[l]
            temp1[l] +=tmp1[l]
            temp1back[l] +=tmp1back[l]
        
    #avn0[j] = pre * (sum_with_exceptions(temp, nsq) / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    for l in range(nmom):
        avn01[l][j]=temp1[l]/nconf
        avn01back[l][j]=temp1back[l]/nconf
    #print(sum_with_exceptions(temp1, nsq))
    print(sum_comps(temp1))
    #avn0[j] = pre * ((temp[0]+temp[1]+temp[2]-temp[3]-temp[4]-temp[5])/6 / (1/3 * np.sqrt((tempdx + tempdy + tempdz) * tempb))) * np.sqrt((4 * mb * md) / (np.exp(-md * j) * np.exp(-mb * (dt - j))))
    
    #+temp1xn + temp1yn + temp1zn+temp1xmn + temp1ymn + temp1zmn


errn01=np.zeros(shape=(nmom,dt))
errn01back=np.zeros(shape=(nmom,dt))

for l in range(nmom):
    for j in range(dt):
        x=0
        xb=0
        for i in range(nconf):
            x=x+(jack(av1n01[l][j],i)-avn01[l][j])**2
            xb=xb+(jack(av1n01back[l][j],i)-avn01back[l][j])**2
        errn01[l,j]=np.sqrt((nconf-1)/nconf*x) 
        errn01back[l,j]=np.sqrt((nconf-1)/nconf*xb) 

for l in range(3):
    plt.errorbar(list(range(dt)),-avn01[:][l],yerr=errn01[l],fmt='x')
    plt.errorbar(list(range(dt)),avn01back[:][l],yerr=errn01back[l],fmt='x')


#antisym-> minus
for l in range(3):
    plt.errorbar(list(range(dt)),avn01[:][l+3],yerr=errn01[l+3],fmt='x')
    plt.errorbar(list(range(dt)),-avn01back[:][l+3],yerr=errn01back[l+3],fmt='x')

# antisym+minus in mmom comp- -> plus
for l in range(3):
    plt.errorbar(list(range(dt)),-avn01[:][l+6],yerr=errn01[l+6],fmt='x')
    plt.errorbar(list(range(dt)),avn01back[:][l+6],yerr=errn01back[l+6],fmt='x')

#anti sym-> minus
for l in range(3):
    plt.errorbar(list(range(dt)),avn01[:][l+9],yerr=errn01[l+9],fmt='x')
    plt.errorbar(list(range(dt)),-avn01back[:][l+9],yerr=errn01back[l+9],fmt='x')
#print(avn01[:][0]-avn01back[:][0])
#plt.yscale("log")
#plt.axis((13,17,10**(-44),10**(-42)))
#plt.savefig('v-3pt-Test-zoom.pdf')


'''
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
