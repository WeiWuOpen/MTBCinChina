

sourceList=["West_Europe"]
cladeNameList=["4.2"]


inputFolder=r"...\script_and_testData\testData\output\calculate_geographical_origin_FEAST\\"
inputRegionFile=r"...\script_and_testData\testData\strain_geographical_labels_to_large_region_mapping.txt"
outputFolder=r"...\script_and_testData\testData\output\calculate_geographical_origin_FEAST\\"


import os
from tqdm import tqdm
import pandas as pd



if os.path.exists(outputFolder) !=True:
    os.mkdir(outputFolder)


inputRegionDict={}
with open(inputRegionFile,"r") as input:
    for ll in input:
        llx=ll.strip().split("\t")
        if llx[0] != "Run_accession":
            if llx[2] not in inputRegionDict.keys():
                inputRegionDict[llx[2]] = llx[3] # {smallRegion:LargeRegion}


# print(inputRegionDict)

for ii in cladeNameList:
    print(f"{ii} start #######")
    inputStep1File=inputFolder+ii+"_SNP_allCDS_smallRegion_matrix.txt"
    with open(inputStep1File,"r") as input:
        lines=input.readlines()
        regionList=lines[0].strip().split("\t")


    outPutList=[]
    sinkNum = 0  # 注意；sink对应样本的编号需要各不相同，source对应样本的编号需要每次从头开始蝙
    for ss in tqdm(regionList):
        largeRegion = inputRegionDict[ss]
        if largeRegion in sourceList:
            outputLine = ss + "\t" + largeRegion + "\t" + "Source" + "\t" + "NA"
        else:
            sinkNum +=1
            outputLine = ss + "\t" + largeRegion + "\t" + "Sink" + "\t" + str(sinkNum)
        outPutList.append(outputLine)
        del outputLine

    with open(outputFolder+ii+"_sourceSink_metadata_"+"singleSource_mergeSmallregion_"+"_".join(sourceList)+".txt","w") as output:
        output.write("SampleID	Env	SourceSink	id"+"\n")
        for oo in outPutList:
            output.write(oo+"\n")

print("finished!!!!")






