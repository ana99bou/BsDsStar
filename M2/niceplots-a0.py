#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 11:15:26 2024

@author: anastasiaboushmelev
"""

import numpy as np
import pandas as pd
#from bokeh.plotting import figure, show, output_file
import matplotlib.pyplot as plt

mb=1.9257122802734448
md0=0.73482632
md1=0.7458868408203149
md2=0.7567413326078149
md3=0.7673809814453149
md4=0.77745605
md5=0.787656860351565
pre1=md0/(2*md1*mb)
pre2=md0/(2*md2*mb)
pre3=md0/(2*md3*mb)
pre4=md0/(2*md4*mb)
pre5=md0/(2*md5*mb)

nsq3=pd.read_csv('./A0/A0-nsq3.txt',sep=' ',header=None)
nsq1=pd.read_csv('./A0/A0-nsq1.txt',sep=' ',header=None)
nsq2=pd.read_csv('./A0/A0-nsq2.txt',sep=' ',header=None)
nsq4=pd.read_csv('./A0/A0-nsq4.txt',sep=' ',header=None)
nsq5=pd.read_csv('./A0/A0-nsq5.txt',sep=' ',header=None)

'''
nsq3plt=pd.read_csv('./Fits/A0-Av-nsq3-Fit.csv',sep='\s')
nsq1plt=pd.read_csv('./Fits/A0-Av-nsq1-Fit.csv',sep='\s')
nsq2plt=pd.read_csv('./Fits/A0-Av-nsq2-Fit.csv',sep='\s')
nsq4plt=pd.read_csv('./Fits/A0-Av-nsq4-Fit.csv',sep='\s')
nsq5plt=pd.read_csv('./Fits/A0-Av-nsq5-Fit.csv',sep='\s')

x3, y3 = [1, 32], [nsq3plt['EffectiveMass'], nsq3plt['EffectiveMass']]
x1, y1 = [1, 32], [nsq1plt['EffectiveMass'], nsq1plt['EffectiveMass']]
x2, y2 = [1, 32], [nsq2plt['EffectiveMass'], nsq2plt['EffectiveMass']]
x4, y4 = [1, 32], [nsq4plt['EffectiveMass'], nsq4plt['EffectiveMass']]
x5, y5 = [1, 32], [nsq5plt['EffectiveMass'], nsq5plt['EffectiveMass']]



reg_low1=nsq1plt['RegLow']
reg_up1=nsq1plt['RegUp']
sigma1=nsq1plt['Error']
reg_low2=nsq2plt['RegLow']
reg_up2=nsq2plt['RegUp']
sigma2=nsq2plt['Error']
reg_low3=nsq3plt['RegLow']
reg_up3=nsq3plt['RegUp']
sigma3=nsq3plt['Error']
reg_low4=nsq4plt['RegLow']
reg_up4=nsq4plt['RegUp']
sigma4=nsq4plt['Error']
reg_low5=nsq5plt['RegLow']
reg_up5=nsq5plt['RegUp']
sigma5=nsq5plt['Error']
'''

plt.xlabel('Time')
plt.ylabel(r'$\widetilde{A}_0$')
#plt.plot(range(96),nsq1[0])
#plt.plot(range(26),avn0y[0:26])
#plt.plot(range(26),avn0z[0:26])
#plt.plot(range(26),avn0[0:26])
#plt.plot(range(96),avn00)
#plt.errorbar(list(range(26)), pre0*nsq0[0][0:26], yerr=pre0*nsq0[1][0:26],ls='none',fmt='x',label='$n^2=0$',color='g')

plt.errorbar(list(range(26)), nsq1[0][0:26], yerr=nsq1[1][0:26],ls='none',fmt='x',label='$n^2=1$',color='b')
plt.errorbar(list(range(26)), nsq2[0][0:26], yerr=nsq2[1][0:26],ls='none',fmt='x',label='$n^2=2$',color='orange')

plt.errorbar(list(range(26)), nsq3[0][0:26], yerr=nsq3[1][0:26],ls='none',fmt='x',label='$n^2=3$',color='brown')
plt.errorbar(list(range(26)), nsq4[0][0:26], yerr=nsq4[1][0:26],ls='none',fmt='x',label='$n^2=4$',color='red')
plt.errorbar(list(range(26)), nsq5[0][0:26], yerr=nsq5[1][0:26],ls='none',fmt='x',label='$n^2=5$',color='magenta')

#plt.plot(x0,y0,color='g')
#plt.fill_between(list(range(47))[int(reg_low0):int(reg_up0+1)], nsq0plt['EffectiveMass']+sigma0, nsq0plt['EffectiveMass']-sigma0, color='g',alpha=0.2)
'''
plt.plot(x1,y1, color='b')
plt.fill_between(list(range(47))[int(reg_low1):int(reg_up1+1)], nsq1plt['EffectiveMass']+sigma1, nsq1plt['EffectiveMass']-sigma1, color='b',alpha=0.2)
plt.plot(x2,y2,color='orange')
plt.fill_between(list(range(47))[int(reg_low2):int(reg_up2+1)], nsq2plt['EffectiveMass']+sigma2, nsq2plt['EffectiveMass']-sigma2, color='orange',alpha=0.2)

plt.plot(x3,y3,color='brown')
plt.fill_between(list(range(47))[int(reg_low2):int(reg_up2+1)], nsq3plt['EffectiveMass']+sigma3, nsq3plt['EffectiveMass']-sigma3, color='brown',alpha=0.2)
plt.plot(x4,y4,color='red')
plt.fill_between(list(range(47))[int(reg_low4):int(reg_up4+1)], nsq4plt['EffectiveMass']+sigma4, nsq4plt['EffectiveMass']-sigma4, color='red',alpha=0.2)

plt.plot(x5,y5,color='magenta')
plt.fill_between(list(range(47))[int(reg_low5):int(reg_up5+1)], nsq5plt['EffectiveMass']+sigma5, nsq5plt['EffectiveMass']-sigma5, color='magenta',alpha=0.2)
'''

plt.annotate(r'$\bf{preliminary}$',xy=(0.7,0.03),xycoords='axes fraction',fontsize=15,color='grey',alpha=.7)
#plt.axis((0,26,0.15,0.36))



#plt.yscale('log')
plt.legend()
plt.savefig('Niceplot-A0.pdf',transparent=True)
