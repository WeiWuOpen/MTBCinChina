



cladeGroupList=[   ["4.2bNOchina","4.2cnONLYchina"] ]


inputFolder=r"...\script_and_testData\testData\output\information_entropy\\"
outputFolder=r"...\script_and_testData\testData\output\information_entropy\\"

import os

import math




if not os.path.exists(outputFolder):
    os.mkdir(outputFolder)


for cp in cladeGroupList:
    klvalueDict1={}
    with open(inputFolder+cp[0]+"_aMTBC_DFandKLvalues.txt","r") as input:
        for ll in input:
            llx=ll.strip().split("\t")
            if llx[1] =="CDS": # 只考虑CDS
                DFvalue1=float(llx[-1]) # 熵变值
                if DFvalue1 >= 0:
                    KLvalue1 = float(llx[4])
                else:
                    KLvalue1 = -float(llx[4]) # 对于熵减的，KL值取为负数
                klvalueDict1[llx[0]] = KLvalue1
            else:
                continue

    klvalueDict2={}
    with open(inputFolder+cp[1]+"_aMTBC_DFandKLvalues.txt","r") as input:
        for ll in input:
            llx=ll.strip().split("\t")
            if llx[1] =="CDS": # 只考虑CDS
                DFvalue2=float(llx[-1]) # 熵变值
                if DFvalue2 >= 0:
                    KLvalue2 = float(llx[4])
                else:
                    KLvalue2 = -float(llx[4]) # 对于熵减的，KL值取为负数
                klvalueDict2[llx[0]] = KLvalue2
            else:
                continue

    del KLvalue1
    del KLvalue2

    ################################################################
    # 寻找两个对比clade中全基因组的CDS的熵减最大CDS的kl值作为原点
    tmpList1 = []
    with open(inputFolder + cp[0] + "_aMTBC_DFandKLvalues.txt", "r") as input:
        for lll in input:
            lllx = lll.strip().split()
            if lllx[0] != "Gene":  # 去除标题
                    if lllx[1] == "CDS":
                        if float(lllx[-1]) < 0:  # 熵减的
                            tmpList1.append(float(lllx[4]))

    clade1zerovalue = max(tmpList1)
    del lll

    tmpList2 = []
    with open(inputFolder +cp[1] + "_aMTBC_DFandKLvalues.txt", "r") as input:
        for lll in input:
            lllx = lll.strip().split()
            if lllx[0] != "Gene":  # 去除标题
                    if lllx[1] == "CDS":
                        if float(lllx[-1]) < 0:  # 熵减的
                            tmpList2.append(float(lllx[4]))

    clade2zerovalue = -max(tmpList2)  # 因为是熵减，取负值

    zerovalue=min(clade1zerovalue,clade2zerovalue)
    print(zerovalue)

    resultDict={}
    for kk in klvalueDict1.keys():
        KLvalue1=klvalueDict1[kk]
        KLvalue2=klvalueDict2[kk]


        difValue= (KLvalue2-zerovalue) / (KLvalue1-zerovalue)


        resultDict[kk] = kk + "\t" + str(KLvalue2) + "\t" + str(difValue)

        del KLvalue1
        del KLvalue2

    with open(outputFolder+cp[1]+"_KLandDifKL.txt","w") as output:
        output.write("Gene"+"\t"+"KL_"+ cp[1]  +  "\t"+"MKL_"+ cp[1]+"\n")
        for rr in resultDict.keys():
            output.write(resultDict[rr]+"\n")


print("finished!!!!")

