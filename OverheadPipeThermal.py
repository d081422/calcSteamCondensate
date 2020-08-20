# -*- coding: utf-8 -*-
import math 

def CalcThermalFlux(tSteam, innerDiameter, outerDiameter, insulaOneDiameter, insulaTwoDiameter,insulaFirstType, insulaSecondType):
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
    tOneMedial = (tSteam + tBetween) / 2;
    tTwoMedial = (tBetween + tSurf)  / 2;

	#step2:计算第一层保温材料的热导率
    if  (tOneMedial <= 400):
      firstThermConduct = 0.056 + 0.0002 * ( tOneMedial - 70);
    else:
      T_medial_mid = 400;
      lambda_L = 0.056 + 0.0002 * ( T_medial_mid - 70);
      firstThermConduct = lambda_L + 0.00036 * (tOneMedial - 400);

	#step3: 计算二层保温材料导热系数
    secondThermConduct = 0.042 + 0.00017 * (tTwoMedial - 70);
    
	#step4:计算保温层外壁与环境间的对流换热系数
    num1 = pow((273 + tSurf)/100,4);
    alpha = 5.67 * blackness / (tSurf - tAmbient) * (num1 - num2) + 72.81 * num3 / num4;

	#step5: 计算保温层外表面温度
    tNum3 = alpha * insulaTwoDiameter;
    firstTsurfNum = 1 / firstThermConduct * tNum1;
    secondTsurfNum = 1 / secondThermConduct * tNum2;
    thirdTsurfNum = 2000 / tNum3;
    tSurfNum = firstTsurfNum * tAmbient + secondTsurfNum * tAmbient + thirdTsurfNum * tSteam;
    tSurfDem = firstTsurfNum +  secondTsurfNum + thirdTsurfNum;
    tSurf = tSurfNum / tSurfDem;

	#step6: 计算保温层间的温度
    tBetweenNum = firstTsurfNum * tAmbient + secondTsurfNum * tSteam + thirdTsurfNum * tSteam;
    tBetween = tBetweenNum / tSurfDem;
	
  #计算管段单位面积散热损失以及焓降
  x_1 = outerDiameter * 1e-3 /(2 * pipeThermalConduct) * math.log(outerDiameter / innerDiameter);
  x_2 = insulaOneDiameter * 1e-3 / (2 * firstThermConduct) * math.log(insulaOneDiameter / outerDiameter);
  x_3 = insulaTwoDiameter * 1e-3 / (2 * secondThermConduct) * math.log(insulaTwoDiameter / insulaOneDiameter);
  x_4 = 1 / alpha;
  #x_5= 1/alpha2;
  heatLoss = (tSteam - tAmbient) / (x_1 + x_2 + x_3 + x_4);
  return heatLoss;
