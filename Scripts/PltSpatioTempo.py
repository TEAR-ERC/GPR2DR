#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 16:13:19 2020

@author: root
"""
from pandas import read_csv as rd
import numpy as np
import matplotlib.pyplot as plt

# model='../TestLSW-onfault/TestLSW-2e/'
# model='../TestLSW-onfault/LSW-2m/'
# model='../TestLSW-onfault/LSW-2n
folder ='/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/'
model = 'test2'
#folder ='/import/freenas-m-05-seissol/dli/ExaHyPE/ShearBand/'
#model = 'band3'
#model = 'nodam'
figout = folder + model+'_st'+'.png'
pos = '100 m'

ndt = 120;
nx = 301;
t = np.linspace(0,12.0,ndt)  

up1=np.zeros((nx,ndt))

for i in range(0,ndt):
    fname = folder+model+'/spatiotemporal_data.'+str(i)+'.csv'
    data1 = rd(fname,sep=',')
    u_all=data1['Q:1']
    x_all=data1["Points:0"]
    
    fname =folder+ model+'/spatiotemporal_data2.'+str(i)+'.csv'
    data1 = rd(fname,sep=',')
    u2_all=data1['Q:1']
    
    up1[:,i] = u_all[:]-u2_all[:]


#%%

t1 = [2,2.5]
x1 = [-2000,-2000.0+0.5*3464.0] 
    

plt.figure(figsize=(7.5,6.5))
plt.title(model)

ax = plt.subplot(111)

plt.pcolor(x_all[::-1],t,up1.transpose(),cmap='viridis',vmin=0,vmax=16)
plt.plot(x1,t1,c='white',label='cs')


plt.ylabel('time (s)')
plt.xlabel('position (m)')
plt.colorbar()
plt.ylim(0,8)
plt.show()
plt.savefig(figout,dpi=150)
