


cladeNameList=["4.2"]



inputRegionFile=r"...\script_and_testData\testData\strain_geographical_labels_to_large_region_mapping.txt"
outputFolder=r"...\script_and_testData\testData\output\calculate_geographical_origin_FEAST\\"


import os
from tqdm import tqdm
import pandas as pd
import copy
import numpy as np



if os.path.exists(outputFolder) !=True:
    os.mkdir(outputFolder)


inputRegionDict={}
with open(inputRegionFile,"r") as input:
    for ll in input:
        llx=ll.strip().split("\t")
        if llx[0] != "Run_accession":
            inputRegionDict[llx[0]] =[llx[2],llx[3]]


for ii in cladeNameList:
    print(f"{ii} start #######")
    inputlistfile=r"...\script_and_testData\testData\sampleList\\"+ii+".txt"

    strainlist = []
    with open(inputlistfile, "r") as input:
        for l in input:
            strainlist.append(l.strip())

    smallRegionDict = {}
    for qq in strainlist:
        if inputRegionDict[qq][0] not in smallRegionDict.keys():
            smallRegionDict[inputRegionDict[qq][0]] = [qq]
        else:
            smallRegionDict[inputRegionDict[qq][0]].append(qq)

    # 读取step1_1的矩阵
    file_path = r'...\script_and_testData\testData\output\calculate_geographical_origin_FEAST\\'+ii+"_SNP_allCDS_matrix.txt"
    df = pd.read_csv(file_path, delim_whitespace=True)
    # 提取index作为基因名称列表
    genes_list = df.index.tolist()
    # 将列转换为字典
    data_dict = df.to_dict(orient='list')

    smallRegion_data_dict={}
    smallRegion_MeanData_dict={}



    for ss in smallRegionDict.keys():
        tmpDataList = []
        tmpDataList = [0] * len(genes_list)  # 添加genes对应数量的0占位
        for tt in data_dict.keys():
            if tt in smallRegionDict[ss]:
                oldDataList = copy.deepcopy(tmpDataList)
                del tmpDataList
                addDataList = data_dict[tt]
                newDataList = np.add(oldDataList,addDataList).tolist() # 两个list对应位置的值依次相加
                tmpDataList = copy.deepcopy(newDataList)
                del oldDataList
                del addDataList
                del newDataList
            else:
                continue

        smallRegion_data_dict[ss]=tmpDataList
        sampleNum = len(smallRegionDict[ss])
        MeanDataList= [ round(x / sampleNum, 0) for x in tmpDataList ] # list中的每个值求均值，然后四舍五入为整数
        smallRegion_MeanData_dict[ss] = MeanDataList
        del tmpDataList
        del MeanDataList
        del sampleNum

    resultDF=pd.DataFrame(smallRegion_MeanData_dict,index=genes_list)
    resultDF.to_csv(outputFolder + ii + '_SNP_allCDS_smallRegion_matrix.txt', sep='\t', index=True)


print("finished!!!!!")














