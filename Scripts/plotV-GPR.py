#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 10:28:14 2018

@author: dli
"""

import numpy as np
#import matplotlib.pyplot as plt
    
station1='body030st-030dp000'
station2='body030st-030dp000'
rank1 = '70'
folder1 = '/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/test2/'
folder2= '/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/nodam/'
rank2 = '55'

#   TPV 28 normal-fault 3 km, strike +/-15 km, depth 0 km
fin1 = open(folder1+'seismo_'+station1+'-rank-'+rank1+'.probe','r')
fin2 = open(folder2+'seismo_'+station1+'-rank-'+rank2+'.probe','r')

#fin1 = open('./TPV28/output2/seismogram_body030st-150dp000-rank-35.probe','r')    
#fin2 = open('./TPV28/waveq-tpv28-body030st-150dp000.dat','r')

dd1 = np.loadtxt(fin1,comments='#',skiprows=0,delimiter=',')                 
dd2 = np.loadtxt(fin2,comments='#',skiprows=0,delimiter=',')

plt.figure()
plt.title(station1)

ax = plt.subplot(411)
plt.title(station1)
plt.plot(dd1[:,0],dd1[:,3])  # x-axis
plt.plot(dd2[:,0],dd2[:,3])
plt.ylabel('u (m/s)')            
plt.legend(['u-plastic','u-elastic'])

plt.subplot(412,sharex=ax)
plt.plot(dd1[:,0],dd1[:,4]) # y-axis
plt.plot(dd2[:,0],dd2[:,4])
plt.ylabel('v (m/s)')            
plt.legend(['v-plastic','v-elastic'])

plt.subplot(413,sharex=ax)
plt.plot(dd1[:,0],dd1[:,26]) # x-axis
plt.plot(dd2[:,0],dd2[:,26])
plt.ylabel('U (m)')            
plt.legend(['U-plastic','U-elastic'])

plt.subplot(414,sharex=ax)
plt.plot(dd1[:,0],dd1[:,27]) # y-axis
plt.plot(dd2[:,0],dd2[:,27])
plt.ylabel('V (m)')  
plt.legend(['V-plastic','V-elastic'])

plt.xlabel('Time (s)')
plt.xlim(0,12)
plt.show()
plt.savefig('obv-'+station1+'.png',dpi=150)      
    
    

        
