# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 19:38:43 2020

@author: ph
"""

from thermalCalculate import CalcOverHeadPipeThermalFlux,calculateDropPressure
import seuif97
import math
from iapws import iapws97
from binaryCalc import dichotomy
import matplotlib.pyplot as plt


def calccondensate(tempe,maxValue,minValue,seconds,press):
  # 计算这一段管段产生的饱和水量
  # 就当是这条管道是稳态的已经饱和的蒸汽

  innerDiameter = 150
  outerDiameter = 159
  oneInsulaThick = 60
  twoInsulaThick = 50
  oneInsulaDiameter =  outerDiameter + 2 * oneInsulaThick
  twoInsulaDiameter = oneInsulaDiameter + 2 * twoInsulaThick
  pipeLength = 87
  insulaType = ['GSL','GLASS']

  #得到单位面积散热损失，w/m^2
  heatLossPerSquare = CalcOverHeadPipeThermalFlux(tempe, innerDiameter, outerDiameter, oneInsulaDiameter, twoInsulaDiameter, insulaType)

  dQ  = heatLossPerSquare * math.pi * twoInsulaDiameter * 1e-3 * 3.6 * seconds / 3600 #计算每米的热损，kJ/s

  initVolume  = math.pi * pow((twoInsulaDiameter/2),2) * 1e-6   #单位：m^3 ，单位体积

  h1 = seuif97.pt2h(press, tempe)  # 输入的焓值，kJ/kg 
  rho = 1 / seuif97.pt2v(press, tempe)   # kg/m^3
  sumheat = h1 * rho * initVolume   #单位：kJ，管道中预留的总热量
  mass = rho * initVolume * 1e-3  #单位：t，管道中总的蒸汽质量
  
  #计算压降
  dropPress = calculateDropPressure(innerDiameter,pipeLength,mass,press,tempe)  # Mpa

  #二分法计算
  outTempe = dichotomy(maxValue,minValue,press,sumheat,dQ,initVolume)  #出口的温度
  press = press - dropPress  #出口的压强

  return outTempe, press


satT = iapws97._TSat_P(1) - 273.15   #单位：℃
tempe = 300
press = 1
count = 0
tempelist = []

while abs(tempe - satT) > 0.5:
  maxValue = tempe
  minValue = tempe - 1
  seconds = 1
  #计算出口的温度和压强
  tempe,press = calccondensate(tempe,maxValue,minValue,seconds,press)
  count = count + 1
  tempelist.append(round(tempe,5))
print(tempelist)

# 画图
plt.plot(tempelist)
plt.show()
  




