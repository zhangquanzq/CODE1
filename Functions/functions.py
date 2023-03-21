# coding=utf-8
import os
import sys
import re

import arcpy
import numpy as np
import calendar
from os import path
from numpy.polynomial import polynomial as poly
from sklearn import ensemble


def ymd2yday(ymd):
    # Convert YearMonthDay to YearDayOfYear. Both input and output are strings.
    # Format of input: 'yyyymmdd'. For example: '20100503'.
    # Format of output: 'yyyyday'. For example: '2010123'
    year = ymd[0:4]
    month = int(ymd[4:6])
    day = int(ymd[6:8])

    month_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if calendar.isleap(int(year)):
        month_day[1] = 29
    monthSum = sum(month_day[0:month - 1])
    dday = monthSum + day
    yday = '{0}{1:03d}'.format(year, dday)

    return yday


def yday2ymd(year_day):
    year = int(year_day[0:4])
    day = int(year_day[4:7])
    month_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if calendar.isleap(year):
        month_day[1] = 29

    i = 0
    while day > 0:
        day = day - month_day[i]
        i = i + 1

    return '{0}{1:>02}{2:>02}'.format(year, i, day + month_day[i - 1])


def resamplePixelSeries(rasterBufferPath, regionExtentPath, resampleDir, cellsizeList, region,
                        downscale='initial'):
    # 限定downscale参数的取值范围.
    if downscale not in ['initial', 'previous']:
        raise ValueError(u"downscale参数的取值必须为字符串 'initial' 或 'previous' 之一.")

    # 按分辨率序列重采样.
    rasterBufferName = path.basename(rasterBufferPath)
    rasterRegionName = rasterBufferName.replace('Buffer', region)
    rasterBufferDir = path.dirname(rasterBufferPath)
    rasterTempPathList = []
    for i in range(len(cellsizeList) - 1):
        cellsizeNow = cellsizeList[i]
        cellsizeNext = cellsizeList[i + 1]

        rasterBufferNext = rasterBufferName.replace('0.01', cellsizeNext)
        rasterBufferResamplePath = path.join(rasterBufferDir, rasterBufferNext)

        rasterNext = rasterRegionName.replace('0.01', cellsizeNext)
        rasterResamplePath = path.join(resampleDir, rasterNext)

        # 两种重采样方式. 1, 上一步重采样结果作为下一步输入数据, 2, 原始MODIS数据生成所有分辨率序列中的数据.
        if downscale == 'initial':
            rasterPath = rasterBufferPath
        elif downscale == 'previous':
            rasterPath = rasterBufferPath if cellsizeNow == cellsizeList[0] else \
                rasterBufferPath.replace('0.01', cellsizeNow)

        arcpy.Resample_management(rasterPath, rasterBufferResamplePath, cellsizeNext, 'BILINEAR')
        arcpy.Clip_management(rasterBufferResamplePath, '#', rasterResamplePath, regionExtentPath)
        rasterTempPathList.append(rasterBufferResamplePath)

    # 删除Buffer过程数据.
    for rasterTempPath in rasterTempPathList:
        arcpy.Delete_management(rasterTempPath)


def ScalingRegress(largePixelImagePath, smallPixelImagePath, fittedSmallPixelImagePath, neighborN):
    # Get the arrays storing the pixel values and pixel coordinates from the image with large pixel.
    largePixelRaster = arcpy.Raster(largePixelImagePath)
    nodata = largePixelRaster.noDataValue
    extentLargePixel = largePixelRaster.extent
    cellsizeXLarge = largePixelRaster.meanCellWidth
    cellsizeYLarge = largePixelRaster.meanCellHeight
    largePixelXVector = np.arange(extentLargePixel.XMin + cellsizeXLarge / 2.0,
                                  extentLargePixel.XMax, cellsizeXLarge)
    largePixelYVector = np.arange(extentLargePixel.YMax - cellsizeYLarge / 2.0,
                                  extentLargePixel.YMin, - cellsizeYLarge)
    largePixelXArray, largePixelYArray = np.meshgrid(largePixelXVector, largePixelYVector)

    largePixelArray = arcpy.RasterToNumPyArray(largePixelImagePath)
    largePixelNanIndex = (largePixelArray == nodata)
    largePixelArray = largePixelArray.astype(float)
    largePixelArray[largePixelNanIndex] = np.nan
    largePixelXArray[largePixelNanIndex] = np.nan
    largePixelYArray[largePixelNanIndex] = np.nan

    # Get the arrays storing the pixel values and pixel coordinates from the image with small pixel.
    smallPixelRaster = arcpy.Raster(smallPixelImagePath)
    nodata = smallPixelRaster.noDataValue
    extentSmallPixel = smallPixelRaster.extent
    cellsizeXSmall = smallPixelRaster.meanCellWidth
    cellsizeYSmall = smallPixelRaster.meanCellHeight
    smallPixelXVector = np.arange(extentSmallPixel.XMin + cellsizeXSmall / 2.0,
                                  extentSmallPixel.XMax, cellsizeXSmall)
    smallPixelYVector = np.arange(extentSmallPixel.YMax - cellsizeYSmall / 2.0,
                                  extentSmallPixel.YMin, - cellsizeYSmall)
    smallPixelXArray, smallPixelYArray = np.meshgrid(smallPixelXVector, smallPixelYVector)

    smallPixelArray = arcpy.RasterToNumPyArray(smallPixelImagePath)
    smallPixelNanIndex = (smallPixelArray == nodata)
    smallPixelArray = smallPixelArray.astype(float)
    smallPixelArray[smallPixelNanIndex] = np.nan
    smallPixelXArray[smallPixelNanIndex] = np.nan
    smallPixelYArray[smallPixelNanIndex] = np.nan

    # Find the top N nearest pixels around each pixel in the large-pixel image, and construct the
    #   regression model between these nearest pixels and their nearest pixels in the small-pixel
    #   image. N equals neighborN, in which the each pixel itself is included.
    fittedSmallPixelArray = np.zeros_like(smallPixelArray) * np.nan
    for ii in range(len(largePixelYVector)):
        for jj in range(len(largePixelXVector)):
            # Skip the pixels of nan values.
            largePixelValue = largePixelArray[ii, jj]
            if np.isnan(largePixelValue):
                continue

            # Get the image values and coordinates of top N nearest pixels around each pixel in the
            #   large-pixel image.
            largePixelX = largePixelXArray[ii, jj]
            largePixelY = largePixelYArray[ii, jj]
            distanceXArray = largePixelX - largePixelXArray
            distanceYArray = largePixelY - largePixelYArray
            distanceArray = np.sqrt(distanceXArray ** 2 + distanceYArray ** 2)
            distanceMax = np.sort(distanceArray, axis=None)[neighborN - 1]
            neighborPixelIndex = (distanceArray <= distanceMax)
            neighborPixelXVector = largePixelXArray[neighborPixelIndex]
            neighborPixelYVector = largePixelYArray[neighborPixelIndex]
            neighborPixelVector = largePixelArray[neighborPixelIndex]

            # Get the nearest pixels of these top N nearest pixels from the small-pixel image.
            nearestPixelVector = np.zeros_like(neighborPixelVector)
            for j in range(len(neighborPixelVector)):
                distanceXArray = neighborPixelXVector[j] - smallPixelXArray
                distanceYArray = neighborPixelYVector[j] - smallPixelYArray
                distanceArray = np.sqrt(distanceXArray ** 2 + distanceYArray ** 2)
                nearestPixelIndex = (distanceArray == np.nanmin(distanceArray))
                nearestPixelVector[j] = smallPixelArray[nearestPixelIndex][0]

            # Regress the image values in two cellsize levels.
            p = np.polyfit(nearestPixelVector, neighborPixelVector, 1)
            distanceXIndexArray = np.abs(smallPixelXArray - largePixelX) <= cellsizeXLarge / 2.0
            distanceYIndexArray = np.abs(smallPixelYArray - largePixelY) <= cellsizeYLarge / 2.0
            distanceIndexArray = (distanceXIndexArray & distanceYIndexArray)
            fittedSmallPixelArray[distanceIndexArray] = \
                p[0] * smallPixelArray[distanceIndexArray] + p[1]

    # Export the adjusted image.
    lowerLeft = arcpy.Point(extentSmallPixel.XMin, extentSmallPixel.YMin)
    reference = smallPixelRaster.spatialReference
    fittedSmallPixelRaster = arcpy.NumPyArrayToRaster(fittedSmallPixelArray, lowerLeft,
                                                      cellsizeXSmall, cellsizeYSmall, nodata)
    arcpy.DefineProjection_management(fittedSmallPixelRaster, reference)
    fittedSmallPixelRaster.save(fittedSmallPixelImagePath)


# TsHARP降尺度模型.
def tsharp(lstNowPath, ndviNowPath, ndviNextPath, predictLstNextPath, cellsizeNext):
    # 当前分辨率下, 用于建立回归关系的LST和NDVI(FVC)数据层.
    lstNowLayer = arcpy.RasterToNumPyArray(lstNowPath, nodata_to_value=np.nan)

    ndviNowLayer = arcpy.RasterToNumPyArray(ndviNowPath) / 10000.0
    fvcNowLayer = 1 - np.power(1 - ndviNowLayer, 0.625)

    # 下一分辨率级别的NDVI(FVC)数据层.
    ndviNextLayer = arcpy.RasterToNumPyArray(ndviNextPath) / 10000.0
    fvcNextLayer = 1 - np.power(1 - ndviNextLayer, 0.625)

    # 排除无效值像元, 非植被区(水体, 等), 以及定义为NDVI <= 0的像元.
    validIndexLayer = np.logical_not(np.isnan(ndviNowLayer) | np.isnan(lstNowLayer) |
                                     np.less_equal(ndviNowLayer, 0))
    fvcNowVector, lstNowVector = fvcNowLayer[validIndexLayer], lstNowLayer[validIndexLayer]

    # 建立回归关系, 获取回归参数.
    p1 = poly.polyfit(fvcNowVector, lstNowVector, 1)

    # 计算LST残差.
    residualLayer = np.zeros_like(lstNowLayer) * np.nan
    residualLayer[validIndexLayer] = lstNowVector - poly.polyval(fvcNowVector, p1)

    # 将残差重采样到下一分辨率级别.
    ndviNowRas = arcpy.sa.Raster(ndviNowPath)
    residualName = path.basename(lstNowPath).replace('.tif', '_residual.tif')
    residualPath = path.join(path.dirname(predictLstNextPath), residualName)
    residualRas = arcpy.NumPyArrayToRaster(residualLayer, ndviNowRas.extent.lowerLeft,
                                           ndviNowRas.meanCellWidth, ndviNowRas.meanCellHeight)
    arcpy.env.snapRaster = ndviNextPath
    arcpy.Resample_management(residualRas, residualPath, cellsizeNext, 'NEAREST')
    arcpy.DefineProjection_management(residualPath, ndviNowRas.spatialReference)

    # 计算下一分辨率级别的LST.
    lstPredictLayer = poly.polyval(fvcNextLayer, p1)
    ndviNextRas = arcpy.sa.Raster(ndviNextPath)
    lstPredictRas = arcpy.NumPyArrayToRaster(lstPredictLayer, ndviNextRas.extent.lowerLeft,
                                             ndviNextRas.meanCellWidth, ndviNextRas.meanCellHeight)
    lstPredictRas = lstPredictRas + arcpy.sa.Raster(residualPath)
    lstPredictRas.save(predictLstNextPath)
    arcpy.DefineProjection_management(predictLstNextPath, ndviNextRas.spatialReference)


# RF降尺度模型.
def rfDownscaling(lstNowPath, factorNowPathList, factorNextPathList, predictLstNextPath,
                  cellsizeNext):
    # 判断两个分辨率级别的环境变量数量是否一致.
    factorNowN, factorNextN = len(factorNowPathList), len(factorNowPathList)
    if factorNowN != factorNextN:
        print('降尺度前后的环境变量个数不一致!')
        os.system("pause")
        sys.exit()

    # 当前分辨率下, 用于建立回归关系的LST数据层.
    lstNowLayer = arcpy.RasterToNumPyArray(lstNowPath, nodata_to_value=np.nan)

    # 当前分辨率下, 用于建立回归关系的环境变量数据层.
    # factorNowPathList = [ndviNowPath, elevNowPath, slpNowPath]
    rowNowN = arcpy.sa.Raster(factorNowPathList[0]).height
    colNowN = arcpy.sa.Raster(factorNowPathList[0]).width
    factorNowArray = np.zeros((rowNowN, colNowN, factorNowN)) * np.nan
    ndviNowIndex = 0
    for i in range(factorNowN):
        factorNowPath = factorNowPathList[i]
        if 'NDVI' in factorNowPath:
            factorLayer = arcpy.RasterToNumPyArray(factorNowPath) / 10000.0
            factorNowArray[:, :, i] = 1 - np.power(1 - factorLayer, 0.625)
            ndviNowIndex = i
        else:
            factorNowArray[:, :, i] = arcpy.RasterToNumPyArray(factorNowPath)

    # 下一分辨率级别的环境变量数据层.
    # factorNextPathList = [ndviNextPath, elevNextPath, slpNextPath]
    rowNextN = arcpy.sa.Raster(factorNextPathList[0]).height
    colNextN = arcpy.sa.Raster(factorNextPathList[0]).width
    factorNextArray = np.zeros((rowNextN, colNextN, factorNextN)) * np.nan
    ndviNextIndex = 0
    for i in range(factorNextN):
        factorNextPath = factorNextPathList[i]
        if 'NDVI' in factorNextPath:
            factorLayer = arcpy.RasterToNumPyArray(factorNextPath) / 10000.0
            factorNextArray[:, :, i] = 1 - np.power(1 - factorLayer, 0.625)
            ndviNextIndex = i
        else:
            factorNextArray[:, :, i] = arcpy.RasterToNumPyArray(factorNextPath)

    # 排除当前分辨率级别的无效值像元, 以及非植被区(水体, 等).
    validIndexNowLayer = np.logical_not(np.isnan(factorNowArray[:, :, ndviNowIndex]) |
                                        np.isnan(lstNowLayer))
    lstNowVector = lstNowLayer[validIndexNowLayer]
    factorNowMatrix = np.zeros((len(lstNowVector), factorNowN)) * np.nan
    for i in range(factorNowN):
        factorNowMatrix[:, i] = factorNowArray[:, :, i][validIndexNowLayer]

    # 排除下一分辨率级别的无效值像元.
    validIndexNextLayer = np.logical_not(np.isnan(factorNextArray[:, :, ndviNextIndex]))
    factorNextMatrix = np.zeros((np.sum(validIndexNextLayer), factorNowN)) * np.nan
    for i in range(factorNextN):
        factorNextMatrix[:, i] = factorNextArray[:, :, i][validIndexNextLayer]

    # 建立当前分辨率级别的RF回归关系.
    rfRegModel = ensemble.RandomForestRegressor(n_estimators=100, oob_score=True)
    rfRegModel.fit(factorNowMatrix, lstNowVector)

    # 计算LST残差.
    residualLayer = np.zeros_like(lstNowLayer) * np.nan
    residualLayer[validIndexNowLayer] = lstNowVector - rfRegModel.predict(factorNowMatrix)

    # 将残差重采样到下一分辨率级别.
    ndviNowRas = arcpy.sa.Raster(factorNowPathList[ndviNowIndex])
    residualName = path.basename(lstNowPath).replace('.tif', '_residual.tif')
    residualPath = path.join(path.dirname(predictLstNextPath), residualName)
    residualRas = arcpy.NumPyArrayToRaster(residualLayer, ndviNowRas.extent.lowerLeft,
                                           ndviNowRas.meanCellWidth, ndviNowRas.meanCellHeight)
    arcpy.env.snapRaster = factorNextPathList[ndviNextIndex]
    arcpy.Resample_management(residualRas, residualPath, cellsizeNext, 'BILINEAR')
    arcpy.DefineProjection_management(residualPath, ndviNowRas.spatialReference)

    # 计算下一分辨率级别的LST.
    lstPredictLayer = np.zeros((rowNextN, colNextN)) * np.nan
    lstPredictLayer[validIndexNextLayer] = rfRegModel.predict(factorNextMatrix)
    ndviNextRas = arcpy.sa.Raster(factorNextPathList[ndviNextIndex])
    lstPredictRas = arcpy.NumPyArrayToRaster(lstPredictLayer, ndviNextRas.extent.lowerLeft,
                                             ndviNextRas.meanCellWidth,
                                             ndviNextRas.meanCellHeight)
    lstPredictRas = lstPredictRas + arcpy.sa.Raster(residualPath)
    lstPredictRas.save(predictLstNextPath)
    arcpy.DefineProjection_management(predictLstNextPath, ndviNextRas.spatialReference)


# 使用TsHARP模型对LST进行逐级降尺度(直接降尺度时逐级降尺度的特例).
def tsharpStepwise(initialLstPath, ndviDir, predictLstPath, cellsizeList):
    # 根据分辨率序列判断降尺度策略.
    strategy = 'Stepwise' if len(cellsizeList) == 2 else 'Straight'

    # 从文件名中获取数据属性信息.
    nameParts = re.split('[._]', path.basename(initialLstPath))
    sourceLst = nameParts[0]
    daynight = nameParts[3] if 'AMSR' in sourceLst else nameParts[3]
    region = nameParts[4] if 'AMSR' in sourceLst else nameParts[4]
    dateNum = nameParts[2] if 'AMSR' in sourceLst else nameParts[1][2:]

    predictLstStr = '_'.join([sourceLst, 'LST', daynight, region, str(dateNum)])

    # 创建存放指定降尺度方法数据的文件夹.
    regionDateDir = path.dirname(initialLstPath)
    regionMethodDir = path.join(regionDateDir, 'TsHARP_{0}'.format(strategy))
    if not path.exists(regionMethodDir):
        os.mkdir(regionMethodDir)

    for i in range(len(cellsizeList) - 1):
        # 当前和下一分辨率级别.
        cellsizeNow, cellsizeNext = cellsizeList[i], cellsizeList[i + 1]

        # 判断当前分辨率级别的降尺度结果是否存在.
        predictLstNextName = '{0}_{1}.tif'.format(predictLstStr, cellsizeNext)
        predictLstNextPath = path.join(regionMethodDir, predictLstNextName)
        if not path.exists(predictLstNextPath):
            # 读取当前分辨率级别需要降尺度的LST数据.
            predictLstNowName = '{0}_{1}.tif'.format(predictLstStr, cellsizeNow)
            predictLstNowPath = path.join(regionMethodDir, predictLstNowName)
            lstNowPath = initialLstPath if cellsizeNow == cellsizeList[
                0] else predictLstNowPath

            # 当前和下一分辨率级别的NDVI(FVC)数据路径.
            ndviNowName = '*NDVI*{0}.tif'.format(cellsizeNow)
            ndviNowPath = glob(path.join(regionDateDir, ndviNowName))[0]
            ndviNextName = path.basename(ndviNowPath).replace(cellsizeNow, cellsizeNext)
            ndviNextPath = path.join(path.dirname(ndviNowPath), ndviNextName)

            tsharp(lstNowPath, ndviNowPath, ndviNextPath, predictLstNextPath, cellsizeNext)

