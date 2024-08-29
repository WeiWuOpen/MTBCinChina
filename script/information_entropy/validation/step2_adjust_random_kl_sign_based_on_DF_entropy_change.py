

inputFolder=r"...\script_and_testData\testData\output\information_entropy\validation\\"
outputFolder=r"...\script_and_testData\testData\output\information_entropy\validation_r\\"

geneGroupList=["all"]




import numpy as np
import os

if not os.path.exists(outputFolder):
    os.mkdir(outputFolder)



for type in geneGroupList:
    outputFolder2=outputFolder+type+"\\"
    if not os.path.exists(outputFolder2):
        os.mkdir(outputFolder2)

    print(type + " start ######")
    fileList = os.listdir(inputFolder + type + "\\")
    resultDict = {}
    for ff in fileList:
        inputklDict={}
        with open(inputFolder + type + "\\" + ff, "r") as input:
            for ll in input:
                llx = ll.strip().split()
                if llx[0] != "Cycle":
                    if float(llx[-1]) >=0:
                        inputklDict[llx[0]]= float(llx[1])
                    else:
                        inputklDict[llx[0]]= -float(llx[1])
        with open(outputFolder2+ff.replace(".txt","")+"_rkl.txt","w") as output:
            output.write("Cycle"+"\t"+"rkl"+"\n")
            for rr in inputklDict.keys():
                # if inputklDict[rr] <0 :
                    # print(str(inputklDict[rr]))
                output.write(rr+"\t"+str(inputklDict[rr])+"\n")


print("finished!!!!!")
