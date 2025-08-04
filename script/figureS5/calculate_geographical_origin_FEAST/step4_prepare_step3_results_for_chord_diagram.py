
cladeNameList=["4.2"]

inputFolder=r"...\script_and_testData\testData\output\calculate_geographical_origin_FEAST\\"
outputFolder=r"...\script_and_testData\testData\output\calculate_geographical_origin_FEAST\\"


import os
from tqdm import tqdm
import pandas as pd


if os.path.exists(outputFolder) !=True:
    os.mkdir(outputFolder)



for ii in cladeNameList:
    with open(inputFolder+ii+"_source_contributions_matrix.txt","r") as input:
        lines=input.readlines()
        sourceName=lines[0].strip()
        sinkDict={}
        continentList=[]
        for ll in lines[1:]:
            llx=ll.strip().split("\t")

            if "Australia" in ll:
                largeRegionName = "Oceania"
                regionName = " ".join(llx[0].split("_")[0:-1])
                continentName = "Oceania"
            elif "Africa" in ll:
                largeRegionName = "Africa"
                regionName = " ".join(llx[0].split("_")[0:-1])
                continentName = "Africa"
            else:  # 其他大地区的标签都是两个单词的
                largeRegionName="_".join(llx[0].split("_")[-2:])
                regionName = " ".join(llx[0].split("_")[0:-2])
                continentName = largeRegionName.split("_")[1]
            regionName=regionName.replace('"','')
            largeRegionName = largeRegionName.replace('"','')
            continentName = continentName.replace('"','')
            continentList.append(continentName)


            if largeRegionName not in sinkDict.keys():
                sinkDict[largeRegionName]= {regionName: [round(float(x),4) for x in llx[1:] ] }  # 取四位小数
            else:
                sinkDict[largeRegionName][regionName] = [round(float(x),4) for x in llx[1:] ]

    # print(sorted(sinkDict.keys()))
    resultList=[]
    sourceNameList=sourceName.split()
    sourceNameList2=[]
    for zz in sourceNameList:
        if "China" in zz:
            sourceNameList2.append( " ".join(zz.split("_")[:-2]).replace('"','')+"[China]" )
        else:
            sourceNameList2.append(" ".join(zz.split("_")[:-2]).replace('"', '') )
    sourceNameList2[-1]="Unknown"
    sourceName2="\t".join(sourceNameList2)
    resultList.append("\t"+sourceName2)
    for qq in sorted(continentList):
        for kk in sorted(sinkDict.keys()):
            if qq in kk:
                for cc in sinkDict[kk].keys():
                    if "China" in kk:
                        name = cc+"[China]"
                    else:
                        name=cc
                    valueList=[ str(v) for v in sinkDict[kk][cc] ]
                    addLine = name+"\t"+"\t".join(valueList)
                    if addLine not in resultList:
                        resultList.append(addLine)
            else:
                continue

    with open(outputFolder+ii+"tmp.txt","w") as output:
        for r in resultList:
            output.write(r+"\n")

    df=pd.read_csv(outputFolder+ii+"tmp.txt",sep="\t",header=None)
    df_T=df.transpose()
    df_T.to_csv(outputFolder+ii+"_source_sink_matrix.txt",header=None,index=False,sep="\t")
    os.remove(outputFolder+ii+"tmp.txt")

print("finished!!!!!")



