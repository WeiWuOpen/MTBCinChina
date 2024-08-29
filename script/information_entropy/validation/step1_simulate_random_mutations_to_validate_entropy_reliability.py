

geneGroupList = ["all"]
cladeNameList=["allONLYchina"]


randomCycleNum=100


import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import math


for cladeName in cladeNameList:
    print(cladeName +  " is working #######")
    inputlistFile = r"...\script_and_testData\testData\sampleList\\" + cladeName + ".txt"
    inputvcfFolder = r"...\script_and_testData\testData\vcf\\"
    inputuncoverFolder = r"...\script_and_testData\testData\uncoverBed\\"
    inputrefseqFile = r"...\script_and_testData\testData\H37Rv_r.fasta"
    inputrefseqSNPFile_new = r"...\script_and_testData\testData\aMTBC_noCanetti_withprd.fasta"
    inputrefseqPOSFile_new = r"...\script_and_testData\testData\all_exceptAnimal_withCanetti_refH37Rv_withprd.pos.txt"

    outputFolder = r"...\script_and_testData\testData\output\information_entropy\validation\\"

    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)

    refseqdict = {}
    with open(inputrefseqFile, "r") as input:
        for l in input:
            if l.strip()[0] != ">":
                refseq = list(l.strip())
    for i in range(1, len(refseq) + 1):
        refseqdict[int(i)] = refseq[i - 1]

    newRefseqdict = {}
    newRefPoslist = []
    with open(inputrefseqSNPFile_new, "r") as input, open(inputrefseqPOSFile_new, "r") as inputpos:
        for l in input:
            if l.strip()[0] != ">":
                newrefseq = list(l.strip())
        for ll in inputpos:
            newRefPoslist.append(ll.strip())
    for i in range(0, len(newrefseq)):
        newRefseqdict[int(newRefPoslist[i])] = newrefseq[i]

    allLongNewRefseqdict = {}
    for k in refseqdict.keys():
        if k in newRefseqdict.keys() and newRefseqdict[k] != "-":
            allLongNewRefseqdict[k] = newRefseqdict[k]
        else:
            allLongNewRefseqdict[k] = refseqdict[k]

    anti_allLongNewRefseqdict = {}
    baseDict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    for baseindex in allLongNewRefseqdict.keys():
        anti_allLongNewRefseqdict[baseindex] = baseDict[allLongNewRefseqdict[baseindex]]


    ################################################################################
    # 相对熵（KL散度）
    def KL(pA,pT,pC,pG,qA,qT,qC,qG): # 前者为参考
        # import math
        e=math.e

        result=0
        pp=0
        pValueList=[pA,pT,pC,pG]
        for iii in pValueList:
            pp+=1
            if iii !=0:
                if pp==1:
                    result += -(pA*math.log(qA,e))  -  (-(pA * math.log(pA, e)))
                elif pp==2:
                    result += -(pT*math.log(qT,e))  -  (-(pT * math.log(pT, e)))
                elif pp==3:
                    result += -(pC*math.log(qC,e))  -  (-(pC * math.log(pC, e)))
                elif pp==4:
                    result += -(pG*math.log(qG,e))  -  (-(pG * math.log(pG, e)))
            else:
                continue
        del pp
        return result

    # 交叉熵
    def H(pA, pT, pC, pG, qA, qT, qC, qG):  # 前者为参考
        # import math
        e = math.e

        result=0
        pp=0
        pValueList=[pA,pT,pC,pG]
        for iii in pValueList:
            pp+=1
            if iii !=0:
                if pp==1:
                    result += -(pA*math.log(qA,e))
                elif pp==2:
                    result += -(pT*math.log(qT,e))
                elif pp==3:
                    result += -(pC*math.log(qC,e))
                elif pp==4:
                    result += -(pG*math.log(qG,e))
            else:
                continue
        del pp
        return result

    # 改变前后熵的差值（值的正负用于判断变化方向）
    def DF(pA,pT,pC,pG,qA,qT,qC,qG): # 前者为参考
        # import math
        e=math.e

        result=0
        pp=0
        pValueList=[pA,pT,pC,pG]
        for iii in pValueList:
            pp+=1
            if iii !=0:
                if pp==1:
                    result +=  (-(qA * math.log(qA, e))) - (-(pA*math.log(pA,e)))
                elif pp==2:
                    result +=  (-(qT * math.log(qT, e))) - (-(pT*math.log(pT,e)))
                elif pp==3:
                    result +=  (-(qC * math.log(qC, e))) - (-(pC*math.log(pC,e)))
                elif pp==4:
                    result +=  (-(qG * math.log(qG, e))) - (-(pG*math.log(pG,e)))
            else:
                continue
        del pp

        return result

    ########################################################################################


    for geneGroupName in geneGroupList:
        print(geneGroupName + " is starting #####")

        inputGenefile = r"...\script_and_testData\testData\geneGroup\\" + geneGroupName + ".txt"


        outputFolder2=outputFolder+geneGroupName +"\\"
        if not os.path.exists(outputFolder2):
            os.mkdir(outputFolder2)

        strainlist = []
        with open(inputlistFile, "r") as input:
            for l in input:
                strainlist.append(l.strip())
        sampleNum = len(strainlist)

        geneDict = {}
        with open(inputGenefile, "r") as input:
            for ll in input:
                llx = ll.strip().split()
                geneDict[llx[0]] = [int(llx[2]), int(llx[3])]
        del llx
        del ll

        vcfdict = {}
        beddict = {}
        for i in strainlist:
                vcfdict[i] = []
                beddict[i] = []
                with open(inputvcfFolder + i + ".recode.ancestralREFseq.vcf", "r") as inputvcf, \
                    open(inputuncoverFolder + i + ".uncover.bed", "r") as inputbed:
                    for line in inputbed:
                        linex = line.strip().split()
                        beddict[i].append(linex[1] + "_" + linex[2])
                    for l in inputvcf:
                        if l.strip()[0] != "#":
                            vcfdict[i].append(l.strip().split()[1] + "_" + l.strip().split()[3] + "_" + l.strip().split()[4])

        randomSNPDict = {}
        for ss in tqdm(strainlist):
            strainGeneNumDict = {}
            snpPOSDict = {}
            for p in vcfdict[ss]:
                uncoverlist = []
                for u in beddict[ss]:
                    if int(u.split("_")[0]) <= int(p.split("_")[0]) <= int(u.split("_")[1]):
                        uncoverlist.append("uncover")
                        break
                    else:
                        continue
                if len(uncoverlist) != 0:
                    continue
                else:
                    snpPOSDict[int(p.split("_")[0])] = [p.split("_")[1], p.split("_")[2]]

            # tmpPosList=[]
            # for g in geneDict.keys():
            #     strainGeneNumDict[g] = 0
            #     llimit = geneDict[g][0]
            #     rlimit = geneDict[g][1]
            #     for qq in snpPOSDict.keys():
            #         if qq not in tmpPosList:
            #             if llimit <= qq <= rlimit:
            #                 tmpPosList.append(qq)
            #                 strainGeneNumDict[g] +=1
            #             else:
            #                 continue
            #         else:
            #             continue
            # del tmpPosList
            # del g
            # del llimit
            # del rlimit
            # del qq
            snpPosList=snpPOSDict.keys()
            snpPosDf=pd.DataFrame(snpPosList,columns=["POS"])
            poss=snpPosDf['POS'].values
            for g in geneDict.keys():
                bool_array = np.zeros(len(snpPosDf), dtype=bool)
                llimit = geneDict[g][0]
                rlimit = geneDict[g][1]
                bool_array |= (poss >= llimit) & (poss <= rlimit)
                strainGeneNumDict[g] = np.sum(bool_array)


            randomSNPDict[ss]={}
            ctn=0
            for ct in range(1,randomCycleNum+1):
                ctn+=1
                np.random.seed(123 + ctn * 31)
                randomSNPDict[ss][ct]={}
                agn=0
                for ag, (start, end) in geneDict.items():
                    agn +=1
                    np.random.seed(123 + ctn * 31 + agn * 17)
                    snpPosList = list(np.random.choice(range(start, end + 1), size=strainGeneNumDict[ag], replace=False)) # 禁止在同一个位点上选择
                    snpBaseList = list(np.random.choice(["A", "T", "C", "G"], size=strainGeneNumDict[ag], replace=True)) # 允许选择相同突变后碱基
                    snpPosBaseDict=dict(zip(snpPosList,snpBaseList))
                    randomSNPDict[ss][ct][ag]=snpPosBaseDict

        del ss
        resultDict={}
        for cy in  tqdm(range(1,randomCycleNum+1)):
            baseRateBYgeneDict = {}
            baseRateBYnongeneDict = {}
            for ss in strainlist:
                for gg in geneDict.keys():
                    snpPOSDict={}
                    for zz in randomSNPDict[ss][cy][gg].keys():
                        snpPOSDict[zz] = randomSNPDict[ss][cy][gg][zz]
                    if gg not in baseRateBYgeneDict.keys():
                        baseRateBYgeneDict[gg] = [[0, 0, 0, 0], [0, 0, 0, 0]]
                    numA = 0
                    numT = 0
                    numC = 0
                    numG = 0
                    refnumA = 0
                    refnumT = 0
                    refnumC = 0
                    refnumG = 0
                    gllimit = geneDict[gg][0]
                    grlimit = geneDict[gg][1]
                    geneLength = grlimit - gllimit + 1

                    for ii in range(gllimit, grlimit + 1):
                        if ii not in snpPOSDict.keys():
                            if allLongNewRefseqdict[ii] == "A":
                                numA += 1
                                refnumA += 1
                            elif allLongNewRefseqdict[ii] == "T":
                                numT += 1
                                refnumT += 1
                            elif allLongNewRefseqdict[ii] == "C":
                                numC += 1
                                refnumC += 1
                            elif allLongNewRefseqdict[ii] == "G":
                                numG += 1
                                refnumG += 1
                            else:
                                print(allLongNewRefseqdict[ii])
                        else:
                            if snpPOSDict[ii] == "A":
                                numA += 1
                            elif snpPOSDict[ii] == "T":
                                numT += 1
                            elif snpPOSDict[ii] == "C":
                                numC += 1
                            elif snpPOSDict[ii] == "G":
                                numG += 1
                            else:
                                print(snpPOSDict[ii])

                            if allLongNewRefseqdict[ii] == "A":
                                refnumA += 1
                            elif allLongNewRefseqdict[ii] == "T":
                                refnumT += 1
                            elif allLongNewRefseqdict[ii] == "C":
                                refnumC += 1
                            elif allLongNewRefseqdict[ii] == "G":
                                refnumG += 1
                            else:
                                print(allLongNewRefseqdict[ii])

                    baseRateBYgeneDict[gg][0][0] += (refnumA / geneLength) / sampleNum
                    baseRateBYgeneDict[gg][0][1] += (refnumT / geneLength) / sampleNum
                    baseRateBYgeneDict[gg][0][2] += (refnumC / geneLength) / sampleNum
                    baseRateBYgeneDict[gg][0][3] += (refnumG / geneLength) / sampleNum
                    baseRateBYgeneDict[gg][1][0] += (numA / geneLength) / sampleNum
                    baseRateBYgeneDict[gg][1][1] += (numT / geneLength) / sampleNum
                    baseRateBYgeneDict[gg][1][2] += (numC / geneLength) / sampleNum
                    baseRateBYgeneDict[gg][1][3] += (numG / geneLength) / sampleNum

            klResultDict_gene = {}
            hResultDict_gene = {}
            dfResultDict_gene = {}
            for gene in baseRateBYgeneDict.keys():
                pA = baseRateBYgeneDict[gene][0][0]
                pT = baseRateBYgeneDict[gene][0][1]
                pC = baseRateBYgeneDict[gene][0][2]
                pG = baseRateBYgeneDict[gene][0][3]
                qA = baseRateBYgeneDict[gene][1][0]
                qT = baseRateBYgeneDict[gene][1][1]
                qC = baseRateBYgeneDict[gene][1][2]
                qG = baseRateBYgeneDict[gene][1][3]

                # print(pA+pT+pC+pG)
                klresult = KL(pA, pT, pC, pG, qA, qT, qC, qG)
                klResultDict_gene[gene] = klresult
                hresult = H(pA, pT, pC, pG, qA, qT, qC, qG)
                hResultDict_gene[gene] = hresult
                dfresult = DF(pA, pT, pC, pG, qA, qT, qC, qG)
                dfResultDict_gene[gene] = dfresult
                del klresult
                del hresult
                del dfresult
            resultDict[cy] =[klResultDict_gene,hResultDict_gene,dfResultDict_gene]


        resultDict2={}
        for c in sorted(resultDict.keys()):
            for gg in resultDict[c][0]:
                outLine="cycle_" + str(c)+"\t"+str(resultDict[c][0][gg])+"\t"+str(resultDict[c][1][gg])+"\t"+str(resultDict[c][2][gg])
                if gg not in resultDict2.keys():
                    resultDict2[gg]=[outLine]
                else:
                    resultDict2[gg].append(outLine)
                del outLine
            del gg
        del c

        for gg in resultDict2.keys():
            with open(outputFolder2+gg+"_ranSNP_cycle_"+str(randomCycleNum)+".txt","w") as output:
                output.write("Cycle"+"\t"+"KLvalue"+"\t"+"Hvalue"+"\t"+"DFvalue"+"\n")
                for ll in resultDict2[gg]:
                    output.write(ll+"\n")


print("finished!!!!!")
































