import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy.optimize import minimize
import scipy
from sympy.utilities.iterables import multiset_permutations


f = h5py.File("BsDsStar.h5", "r")


nSqs=[0,1,2,3,4,5]
cycl=['X','Y','Z']
momss=[[[0,0,0]],[[1,0,0],[-1,0,0]],[[1,1,0],[1,-1,0],[-1,-1,0]],[[1,1,1],[-1,1,1],[-1,-1,1],[-1,-1,-1]],[[2,0,0],[-2,0,0]],[[2,1,0],[2,-1,0],[-2,1,0],[-2,-1,0]]]

mom=[]

for nSq in nSqs:
    moms=momss[nSq]
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
    mom.append(temp)
print(mom)


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


iNs=0

for nsqStr in mom:
    for threePt in nsqStr:
        a,b,c=parse_string(threePt)
        res=eps_over_k(a,b,c)
        print('n^2=',iNs, 'iV, iJ, iMom=', a,b,c, 'eps/k=', res)
    iNs+=1     

ptmom=['0_0_0','1_0_0','1_1_0','1_1_1', '2_0_0','2_1_0']
