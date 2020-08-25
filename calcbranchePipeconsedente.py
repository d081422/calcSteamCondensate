# -*- coding: utf-8 -*-
from iapws import iapws97,IAPWS97
import matplotlib.pyplot as plt
from binaryCalc import calcOverHeatSteam,calcSatSteam
import math
import seuif97


satT = iapws97._TSat_P(1) - 273.15   #单位：℃
tempe = 300
press = 1
count = 0
tempelist,pressList = [],[]
sumTime = 0

innerDiameter = 150
outerDiameter = 159
oneInsulaThick = 60
twoInsulaThick = 50
pipeLength = 87
oneInsulaDiameter = outerDiameter + 2 * oneInsulaThick
twoInsulaDiameter = oneInsulaDiameter + 2 * twoInsulaThick
insulaType = ['GSL','GLASS']

#原管道中的蒸汽质量
initVolume = math.pi * pow((twoInsulaDiameter / 2), 2) * 1e-6  # 单位：m^3 ，单位体积
rho1 = 1 / seuif97.pt2v(press, tempe)  # kg/m^3
massSteam = rho1 * initVolume   # 单位：kg，管道中总的蒸汽质量
sumtime,heatLossPerSquareList,heatLossPerMeterList, rhoList, hList, sumHeatList = [],[],[],[],[],[]


# 初始的温度和压强下每米管道的蒸汽总热量
initVolume = math.pi * pow((twoInsulaDiameter / 2), 2) * 1e-6  # 单位：m^3 ，单位体积
h1 = seuif97.pt2h(press, tempe)  # 输入的焓值，kJ/kg
rho = 1 / seuif97.pt2v(press, tempe)  # kg/m^3
sumheat = h1 * rho * initVolume  # 单位：kJ，管道中预留的总热量


# 过热蒸汽到饱和蒸汽的过渡计算
while abs(tempe - satT) > 0.5:
  maxValue = tempe
  minValue = tempe - 1
  seconds = 1
  #计算出口的温度和压强
  tempe, press,heatLossPerSquare,dQ,rho,h1,sumheat = calcOverHeatSteam(tempe, maxValue, minValue, press,seconds,innerDiameter,
                                                                       outerDiameter,oneInsulaDiameter,twoInsulaDiameter,pipeLength,insulaType,sumheat)
  count = count + 1
  rho = 1 / seuif97.pt2v(press, tempe)
  h1 = seuif97.pt2h(press, tempe)  # 输入的焓值，kJ/kg
  tempelist.append(round(tempe,11))
  pressList.append(round(press,15))
  heatLossPerSquareList.append(round(heatLossPerSquare,15))
  heatLossPerMeterList.append(round(dQ,12))
  rhoList.append(round(rho,11))
  hList.append(round(h1,11))
  sumHeatList.append(round(sumheat,11))
  sumTime = sumTime + seconds
  sumtime.append(sumTime)

#绘图
from pylab import mpl
mpl.rcParams['font.sans-serif'] = ['SimHei']
# plt.plot(sumtime,heatLossPerSquareList)
# plt.xlabel("时间(S)")
# plt.ylabel("最外层保温单位面积散热损失(w/m^2)")

# plt.plot(sumtime,heatLossPerMeterList)
# plt.xlabel("时间(S)")
# plt.ylabel("最外层保温单位长度散热损失(w/m)")

# plt.plot(sumtime,rhoList)
# plt.xlabel("时间(S)")
# plt.ylabel("密度(kg/m^3)")


# plt.plot(sumtime,hList)
# plt.xlabel("时间(S)")
# plt.ylabel("焓值(kJ/kg)")

plt.plot(sumtime,tempelist)
plt.xlabel("时间(S)")
plt.ylabel("温度(℃)")

# plt.plot(sumtime,sumHeatList)
# plt.xlabel("时间(S)")
# plt.ylabel("总热量(℃)")

plt.show()



#蒸汽变为饱和蒸汽以后的计算
#计算凝结水量
satWaterH = IAPWS97(P = press, x=0).h  # 饱和水焓，kJ/kg
satWaterH = IAPWS97(P = press,x = 0).rho # 饱和水密度，kg/m^3
satVaporH = IAPWS97(P = press, x=1).h  #饱和蒸汽的焓，kJ/kg

# 每米的散热损失
#print("饱和蒸汽温度：",tempe)
seconds = 1
dQ = calcSatSteam(satT,innerDiameter,outerDiameter,oneInsulaDiameter,twoInsulaDiameter,pipeLength,insulaType)  # kJ/s
# print(dQ)
#计算每秒产生的凝结水量
waterMass = dQ / satWaterH # kg/s
print(waterMass)

seconds = (massSteam / waterMass) / 60
print(seconds)
#计算原管道中的蒸汽质量
massSteam = massSteam - waterMass #饱和蒸汽质量，kg

print(massSteam)


  




