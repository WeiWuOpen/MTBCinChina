


inputvcffolder=r"...\script_and_testData\testData\vcf\\"
inputAnnotationFile=r"...\script_and_testData\testData\H37Rv.annotation_all.txt"

cladeNameList=["4.2"]

outputFolder=r"...\script_and_testData\testData\output\\"


import os
from tqdm import tqdm
import pandas as pd


if os.path.exists(outputFolder) !=True:
    os.mkdir(outputFolder)


allCDSDict = {}
with open(inputAnnotationFile, "r") as input:
    for l in input:
        lx = l.strip().split()
        if lx[2] =="CDS":
            allCDSDict[int(lx[3])] = [lx[0],int(lx[4])]


for ii in cladeNameList:
    print(f"{ii} start #######")
    inputlistfile=r"G:\aaaworksite\mtbc基因组数据总结\\"+ii+".txt"

    strainlist=[]
    with open(inputlistfile,"r") as input:
        for l in input:
            strainlist.append(l.strip())

    filterStrainlist=[]
    resultDict={}
    for ss in tqdm(strainlist):
        try:
            poslistSNP=[]
            cdsSNPnumDict={}
            for p in allCDSDict.keys():
                cdsSNPnumDict[allCDSDict[p][0]]=0
            with open(inputvcffolder + ss + ".recode.ancestralREFseq.vcf", "r") as input:
                for ll in input:
                    llx = ll.strip().split()
                    if llx[0] != "#":
                        if llx[0] == "aMTBC":
                            poslistSNP.append(int(llx[1]))
            filterStrainlist.append(ss)
            sorted_poslistSNP=sorted(poslistSNP)
            lastIndex=0
            cdsStartPosList=sorted(allCDSDict.keys())
            for pp in sorted_poslistSNP:
                for si in range(lastIndex,len(cdsStartPosList)):
                    lpos=cdsStartPosList[si]
                    rpos=allCDSDict[lpos][1]
                    if si != len(cdsStartPosList)-1:
                        nextlpos=cdsStartPosList[si+1]
                        nextrpos=allCDSDict[nextlpos][1]
                        if pp < lpos:
                            print("wrong in SNP pos !!!!!")
                        elif lpos <= pp <= rpos:
                            lastIndex=si
                            cdsSNPnumDict[allCDSDict[lpos][0]] +=1
                            break
                        elif rpos < pp < nextlpos:
                            break  # SNP在基因间区
                        else:
                            continue
                    else: # 最后一个CDS
                        if pp < lpos:
                            print("wrong in SNP pos !!!!!")
                        elif lpos <= pp <= rpos:
                            lastIndex=si
                            cdsSNPnumDict[allCDSDict[lpos][0]] +=1
                            break
                        else:
                            continue



            resultDict[ss]=cdsSNPnumDict

        except FileNotFoundError:
            print(ss + ".recode.ancestralREFseq.vcf is NOT found !!")

    strainNum=len(filterStrainlist)
    print(f"{strainNum} / {len(strainlist)}  sample exist !!!!!")
    resultDF=pd.DataFrame(resultDict)
    resultDF.to_csv(outputFolder+ii+'_SNP_allCDS_matrix.txt', sep='\t', index=True)









