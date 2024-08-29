

inputFolder=r"...\script_and_testData\testData\output\information_entropy\validation_r\\"
outputFolder=r"...\script_and_testData\testData\output\information_entropy\validation_r\\"

geneGroupList=["all"]




import numpy as np
import os

for type in geneGroupList:
    print(type+" start ######")
    fileList=os.listdir(inputFolder+type+"\\")
    resultDict={}
    for ff in fileList:
        inputklList=[]
        with open(inputFolder+type+"\\"+ff,"r") as input:
            for ll in input:
                llx=ll.strip().split()
                if llx[0] !="Cycle":
                    inputklList.append(float(llx[1]))
        mean=np.mean(inputklList)
        std=np.std(inputklList)
        resultDict[ff.split("_")[0]]=[mean,std]
        del mean
        del std
        del inputklList
    with open(outputFolder+type+"rKL_MEANandSTD.txt","w") as output:
        output.write("Gene"+"\t"+"rKLmean"+"\t"+"rKLstd"+"\n")
        for rr in resultDict.keys():
            output.write(rr+"\t"+str(resultDict[rr][0])+"\t"+str(resultDict[rr][1])+"\n")\

print("finished!!!!!")




