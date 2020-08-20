# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 18:20:07 2020

@author: ph
"""
import seuif97

def  dichotomy(maxValue,minValue,press,sumheat,dQ,initVolume):  
  error1 = -1
  error2 = 1
  k = minValue
  count = 0
  while (abs(error1 - error2) > 0.00001):
    if (error1 * error2 < 0):
      maxValue = maxValue
      minValue = k
      k = (maxValue + minValue) / 2
    
      h2 = seuif97.pt2h(press,k)
      rho2 = 1 /seuif97.pt2v(press,k)
      
      h2up = seuif97.pt2h(press,maxValue)
      rho2up = 1 / seuif97.pt2v(press,maxValue)
      
      error1 = (sumheat - dQ * 1) -(initVolume * rho2 * h2)
      error2 = (sumheat - dQ * 1) -(initVolume * rho2up * h2up) 
      count = count + 1
            
    
    else:
      maxValue = k
      minValue = minValue
      k = (maxValue + minValue) / 2
      
      h2 = seuif97.pt2h(press,k)
      rho2 = seuif97.pt2v(press,k)
      
      h2down = seuif97.pt2h(press,minValue)
      rho2down = 1 / seuif97.pt2v(press,minValue)
      
      error1 = (sumheat - dQ * 1) -(initVolume * h2 * rho2)
      error2 = (sumheat - dQ * 1) -(initVolume * h2down * rho2down)
      count = count + 1
  return k 
    
# k = dichotomy(300,299,1,1324,0.1773,0.112)

  

    
