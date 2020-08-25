# -*- coding: utf-8 -*-
import math
import seuif97 


#根据材料类型计算保温材料的导热系数

def CalcTypeInsula(insulaType,tm):
  #增加判断进行计算
  if insulaType == 'NAMI':
    #计算纳米气凝胶的导热系数
    conduct = 0.0175 * (1 + 0.0004483516 * tm + 0.000003164835 * pow(tm,2));  
  elif insulaType =='GSL':
    #计算硅酸铝保温材料的导热系数
    if  (tm <= 400):
      conduct = 0.056 + 0.0002 * ( tm - 70);
    else:
      T_medial_mid = 400;
      lambda_L = 0.056 + 0.0002 * ( T_medial_mid - 70);
      conduct = lambda_L + 0.00036 * (tm - 400);
  elif insulaType == 'GLASS':
    #计算高温玻璃棉的导热系数
    conduct = 0.042 + 0.00017 * (tm - 70);
  else:
    print("error")
  return conduct


# 计算架空管道的单位面积热损
def CalcOverHeadPipeThermalFlux(tSteam, innerDiameter, outerDiameter, insulaOneDiameter, insulaTwoDiameter, insulaType):
  tBetween = 50
  tSurf = 30
  tAmbient = 25
  windSpeed= 2.35
  blackness = 0.3
  pipeThermalConduct = 38.65
  num2 = pow((273 + tAmbient)/100,4)
  num3 = pow(windSpeed,0.6)
  num4 = pow(insulaTwoDiameter,0.4)
  tNum1 = math.log(insulaOneDiameter / outerDiameter)
  tNum2 = math.log(insulaTwoDiameter / insulaOneDiameter)

  #初设tSurf和tBetween，迭代计算新的tSuf和tBetween,但是迭代次数只有5次，需要重新考虑
  #架空蒸汽管道
  for i in range(5):
    
    #step1:计算保温层间的平均温度
    tOneMedial = (tSteam + tBetween) / 2
    tTwoMedial = (tBetween + tSurf)  / 2

	#step2:计算第一层保温材料的热导率
    # step2: 计算内层保温材料的导热系数——硅酸铝管壳
    firstThermConduct = CalcTypeInsula(insulaType[0],tOneMedial)
      
    # step3: 计算外层保温的导热系数——高温玻璃棉板
    secondThermConduct = CalcTypeInsula(insulaType[1],tOneMedial)
  
	#step4:计算保温层外壁与环境间的对流换热系数
    num1 = pow((273 + tSurf)/100,4)
    alpha = 5.67 * blackness / (tSurf - tAmbient) * (num1 - num2) + 72.81 * num3 / num4

	#step5: 计算保温层外表面温度
    tNum3 = alpha * insulaTwoDiameter
    firstTsurfNum = 1 / firstThermConduct * tNum1
    secondTsurfNum = 1 / secondThermConduct * tNum2
    thirdTsurfNum = 2000 / tNum3
    tSurfNum = firstTsurfNum * tAmbient + secondTsurfNum * tAmbient + thirdTsurfNum * tSteam
    tSurfDem = firstTsurfNum + secondTsurfNum + thirdTsurfNum
    tSurf = tSurfNum / tSurfDem

	#step6: 计算保温层间的温度
    tBetweenNum = firstTsurfNum * tAmbient + secondTsurfNum * tSteam + thirdTsurfNum * tSteam
    tBetween = tBetweenNum / tSurfDem
	
  #step7：计算管段单位面积散热损失以及焓降
  x_1 = outerDiameter * 1e-3 /(2 * pipeThermalConduct) * math.log(outerDiameter / innerDiameter)
  x_2 = insulaOneDiameter * 1e-3 / (2 * firstThermConduct) * math.log(insulaOneDiameter / outerDiameter)
  x_3 = insulaTwoDiameter * 1e-3 / (2 * secondThermConduct) * math.log(insulaTwoDiameter / insulaOneDiameter)
  x_4 = 1 / alpha
  #x_5= 1/alpha2;
  heatLossOneSquare = (tSteam - tAmbient) / (x_1 + x_2 + x_3 + x_4)
  return heatLossOneSquare


def calculateDropPressure(innerDiameter,pipeLength,flowMass,press,tempe):
  #flowRate:流速(m/s);
  #innerDiameter: 管道内径（mm）
  #mu :介质的动力粘度（pa*s）
  #V_steam:介质的比容（m^3/kg）
  press = press + 0.1
  absRough= 0.0457
  KinematicViscosity = seuif97.pt(press, tempe,25) # 运动粘度
  volume = seuif97.pt2v(1, 300)
  dynamicViscosity  = KinematicViscosity / volume  #Pas
  rhoSteam  = 1 / volume
  
  #计算雷诺数，确定管道的摩擦阻力系数
  reynold = flowMass * innerDiameter * 1e-3 / (dynamicViscosity * volume)  #雷诺数
  a = -2 * math.log10(absRough /(innerDiameter  * 3.7) + 12 / reynold);
  b = -2 * math.log10(absRough /(innerDiameter  * 3.7) + 2.51 * a / reynold);
  mediumValueNum1 = pow((a - 4.781),2)
  mediumValueNum2 = b- 2 * a + 4.781
  lambda1 = pow(4.781-mediumValueNum1 / mediumValueNum2,-2)
  #k = 1.15; % 修正系数
  dropPressure =  lambda1 * pipeLength * rhoSteam * pow(flowMass,2) * 1e-6 /(2 * innerDiameter * 1e-3 ) #MPa
  return dropPressure


def calculatDropTempe(QLoss,twoInsulaDiameter,pipeLength,flowMass,tempe,press):
  press = press + 0.1
  #inputEnthalpy = seuif97.pt2h(press, tempe)
  hDropMainPipe = math.pi * twoInsulaDiameter * 1e-3 * QLoss * pipeLength * 3.6/ (flowMass * 1000)
  heatLoss = hDropMainPipe * flowMass * 1e-3  # 单位GJ
  #inputHeat = inputEnthalpy * flowMass * 1e-3
  return heatLoss
  

