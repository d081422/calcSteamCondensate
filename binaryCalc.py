# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 18:20:07 2020

@author: ph
"""
from thermalCalculate import CalcOverHeadPipeThermalFlux,calculateDropPressure
import seuif97
import math

def dichotomy(maxValue,minValue,press,sumheat,dQ,initVolume):
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

# 计算支线过热蒸汽的输出温度和压强
def calcOverHeatSteam(tempe, maxValue, minValue, press,seconds,innerDiameter,outerDiameter,oneInsulaDiameter,twoInsulaDiameter,pipeLength,insulaType,sumheat):
  # 计算这一段管段产生的饱和水量
  # 就当是这条管道是稳态的已经饱和的蒸汽
  #dQ = calcheatLossPermeter(tempe,seconds)
  heatLossPerSquare = CalcOverHeadPipeThermalFlux(tempe, innerDiameter, outerDiameter, oneInsulaDiameter,twoInsulaDiameter, insulaType)
  dQ = heatLossPerSquare * math.pi * twoInsulaDiameter * 1e-3 * 3.6 * seconds / 3600  # 计算每米的热损，kJ
  initVolume = math.pi * pow((twoInsulaDiameter / 2), 2) * 1e-6  # 单位：m^3 ，单位体积
  h1 = seuif97.pt2h(press, tempe)  # 输入的焓值，kJ/kg
  rho = 1 / seuif97.pt2v(press, tempe)  # kg/m^3
  #sumheat = h1 * rho * initVolume  # 单位：kJ，管道中预留的总热量
  mass = rho * initVolume * 1e-3  # 单位：t，管道中总的蒸汽质量
  # 计算压降
  #dropPress = calculateDropPressure(innerDiameter, pipeLength, mass, press, tempe)  # Mpa
  # 二分法计算
  outTempe = dichotomy(maxValue, minValue, press, sumheat, dQ, initVolume)  # 出口的温度
  sumheat = sumheat - dQ
  #press = press - dropPress  # 出口的压强
  return outTempe, press, heatLossPerSquare,dQ,rho,h1,sumheat


#计算支线蒸汽饱和之后的温度和压强
def calcSatSteam(tempe,innerDiameter,outerDiameter,oneInsulaDiameter,twoInsulaDiameter,pipeLength,insulaType):
  heatLossPerSquare = CalcOverHeadPipeThermalFlux(tempe, innerDiameter, outerDiameter, oneInsulaDiameter,
                                               twoInsulaDiameter, insulaType)
  dQ = heatLossPerSquare * math.pi * twoInsulaDiameter * 1e-3 * 3.6 / 3600  # 计算每米的热损，kJ

  return dQ
