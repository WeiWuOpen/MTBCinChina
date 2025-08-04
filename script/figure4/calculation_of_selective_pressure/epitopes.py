

cladeNameList=["4.2cnONLYchina"]
epiTypeList=["TcellEpitopes"]


cycleNum =1000 # 模拟的次数
randomseed= 123


import os
import numpy as np
from scipy import stats
import math
from tqdm import tqdm

inputTamuraRateFile = r"...\script_and_testData\testData\allstrains_tamura_model_parameters_mean.txt"
inputcodonPOSsnpRateFile = r"...\script_and_testData\testData\allstrains_mean_mutation_rate_per_codon_position.txt"
inputTi_TvMeanRateFile = r"...\script_and_testData\testData\allstrains_transition_transversion_ratio_mean.txt"
inputCodonBiasFile = r"...\script_and_testData\testData\H37Rv_codon_usage_bias.txt"

inputAllGeneFile = r"...\script_and_testData\testData\\H37Rv.annotation_all.txt"
inputuncoverFolder = r"...\script_and_testData\testData\uncoverBed\\"
inputrefseqFile = r"...\script_and_testData\testData\H37Rv_r.fasta"


tamuraRateDict={}
with open(inputTamuraRateFile,"r") as input:
    for ll in input:
        llx=ll.strip().split()
        if llx[0] !="Rate1":
            tamuraRateDict["rate1"] = float(llx[0])
            tamuraRateDict["rate2"] = float(llx[1])


codonPosRateDict={}
with open(inputcodonPOSsnpRateFile,"r") as input:
    for l in input:
        lx=l.strip().split()
        if lx[0] != "FirstBaseSNPRate":
            codonPosRateDict["1st"] = float(lx[0])
            codonPosRateDict["2nd"] = float(lx[1])
            codonPosRateDict["3rd"] = float(lx[2])

titvDict = {}
with open(inputTi_TvMeanRateFile, "r") as input:
    for ll in input:
        llx = ll.strip().split()
        if llx[0] != "TransitionRate":
            titvDict["ti"] = float(llx[0])
            titvDict["tv"] = float(llx[1])
            titvDict["tvtype1"] = float(llx[2])  # AT & CG
            titvDict["tvtype2"] = float(llx[3])  # AC & TG
# 转换
transitionList = ["AG", "GA", "CT", "TC"]

# 颠换
transversionList = ["AC", "AT", "TA", "TG", "GC", "GT", "CA", "CG"]
transversionList_type1 = ["AT", "TA", "GC", "CG"]
transversionList_type2 = ["AC", "TG", "GT", "CA"]

# Taruma模型中对应两种速率的替换模式
rate1List = ["AG", "TG", "CG", "AC", "TC", "GC"]
rate2List = ["TA", "CA", "GA", "AT", "CT", "GT"]

allRateDict = {}
rateα = titvDict["ti"]
rateβ1 = titvDict["tvtype1"]
rateβ2 = titvDict["tvtype2"]
rateΘ1 = tamuraRateDict["rate1"]
rateΘ2 = tamuraRateDict["rate2"]
allRateDict["A"] = 2 * rateβ1 * rateΘ2 + 2 * rateβ2 * rateΘ1 + rateα * rateΘ1
allRateDict["T"] = 2 * rateβ1 * rateΘ2 + 2 * rateβ2 * rateΘ1 + rateα * rateΘ1
allRateDict["C"] = 2 * rateβ1 * rateΘ1 + 2 * rateβ2 * rateΘ2 + rateα * rateΘ2
allRateDict["G"] = 2 * rateβ1 * rateΘ1 + 2 * rateβ2 * rateΘ2 + rateα * rateΘ2
snpRateDict = {}
snpRateDict["AT"] = 2 * rateβ1 * rateΘ2
snpRateDict["AC"] = 2 * rateβ2 * rateΘ1
snpRateDict["AG"] = rateα * rateΘ1
snpRateDict["TA"] = 2 * rateβ1 * rateΘ2
snpRateDict["TC"] = rateα * rateΘ1
snpRateDict["TG"] = 2 * rateβ2 * rateΘ1
snpRateDict["CA"] = 2 * rateβ2 * rateΘ2
snpRateDict["CT"] = rateα * rateΘ2
snpRateDict["CG"] = 2 * rateβ1 * rateΘ1
snpRateDict["GA"] = rateα * rateΘ2
snpRateDict["GT"] = 2 * rateβ2 * rateΘ2
snpRateDict["GC"] = 2 * rateβ1 * rateΘ1




refseqdict = {}
with open(inputrefseqFile, "r") as input:
    for l in input:
        if l.strip()[0] != ">":
            refseq = list(l.strip())
for i in range(1, len(refseq) + 1):
    refseqdict[int(i)] = refseq[i - 1]


## 这个其实是用来翻译密码子的
condonBiasDict = {}
condonBiaslist = []
with open(inputCodonBiasFile, "r") as input:
    for l in input:
        lx = l.strip().split()
        if lx != [] and lx[0] != "Coding":
            condonBiaslist.append(l.strip().replace("U", "T"))

for ll in condonBiaslist:
    llx = ll.strip().split()
    condonBiasDict[llx[0]] = [llx[1], llx[3]]
    condonBiasDict[llx[6]] = [llx[7], llx[9]]
    condonBiasDict[llx[12]] = [llx[13], llx[15]]
    condonBiasDict[llx[18]] = [llx[19], llx[21]]



allgeneDict = {}
with open(inputAllGeneFile, "r") as input:
    for l in input:
        lx = l.strip().split()
        if lx[2] =="CDS":
            if lx[0][-1] == "c":
                allgeneDict[lx[0]] = lx[3] + "_" + lx[4] + "_" + "c"
            else:
                allgeneDict[lx[0]] = lx[3] + "_" + lx[4] + "_" + "f"


##############################################################################################################

for cladeName in cladeNameList:
    print(cladeName + " is working #######")

    inputlistFile = r"...\script_and_testData\testData\sampleList\\" + cladeName + ".txt"
    inputvcfFolder = r"...\script_and_testData\testData\vcf_4.2_CRMA\\"
    inputrefseqSNPFile_new = r"...\script_and_testData\testData\\"+cladeName.replace("ONLYchina","").replace("NOchina","")+"_CRMA.fasta"

    if cladeName[0] == "2":
        inputrefseqPOSFile_new = r"...\script_and_testData\testData\allL2_withGeolocInfo_refaMTBC_withprd.pos.txt"
    elif cladeName[0] == "4":
        inputrefseqPOSFile_new = r"...\script_and_testData\testData\allL4_withGeolocInfo_refaMTBC_withprd.pos.txt"
    else:
        print("wrong in inputrefseqPOSFile_new!!!!!!")



    strainlist = []
    with open(inputlistFile, "r") as input:
        for l in input:
            strainlist.append(l.strip())
    sampleNum = len(strainlist)
    # print(sampleNum)


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



    vcfdict = {}
    beddict = {}
    for i in strainlist:
        vcfdict[i] = []
        beddict[i] = []
        with open(inputvcfFolder + i + ".recode.ancestralREFseq.vcf", "r") as inputvcf, open(
                inputuncoverFolder + i + ".uncover.bed", "r") as inputbed:
            for line in inputbed:
                linex = line.strip().split()
                beddict[i].append(linex[1] + "_" + linex[2])
            for l in inputvcf:
                if l.strip()[0] != "#":
                    vcfdict[i].append(l.strip().split()[1] + "_" + l.strip().split()[3] + "_" + l.strip().split()[4])



    ########## 全基因组上全部CDS的全部待选ns，s位点数目与待选optimal,nonoptimal位点数目 ####################################
    print("##### start whole genome genes candidate ns s et al. ########## ")
    allNSsiteNum = 0
    allSsiteNum = 0
    alloptimalsiteNum = 0
    allnonoptimalsiteNum = 0
    allnormalsiteNum = 0
    allCDSLength = 0
    ggg = 0
    for gg in allgeneDict.keys():
        ggg += 1
        # print(str(ggg)+"/"+str(len(allgeneDict.keys())))
        for pp in range(int(allgeneDict[gg].split("_")[0]), int(allgeneDict[gg].split("_")[1]) + 1, 3):
            allCDSLength += 3
            if allgeneDict[gg].split("_")[2] == "f":
                codon = allLongNewRefseqdict[pp] + allLongNewRefseqdict[pp + 1] + allLongNewRefseqdict[pp + 2]
            elif allgeneDict[gg].split("_")[2] == "c":
                codon = anti_allLongNewRefseqdict[pp + 2] + anti_allLongNewRefseqdict[pp + 1] + \
                        anti_allLongNewRefseqdict[pp]
            else:
                print("WRONG IN allgeneDict!!!!!!!")

            for yi in range(1, 4):
                aminoDict = {}
                if yi == 1:
                    allrate = allRateDict[codon[0]]
                    for k in condonBiasDict.keys():
                        if codon[1:] == k[1:] and codon[0] != k[0]:
                            snpRate = snpRateDict[codon[0] + k[0]]
                            aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                            del snpRate
                        else:
                            continue
                    del allrate
                elif yi == 2:
                    allrate = allRateDict[codon[1]]
                    for k in condonBiasDict.keys():
                        if codon[0] == k[0] and codon[-1] == k[-1] and codon[1] != k[1]:
                            snpRate = snpRateDict[codon[1] + k[1]]
                            aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                            del snpRate
                        else:
                            continue
                    del allrate
                elif yi == 3:
                    allrate = allRateDict[codon[2]]
                    for k in condonBiasDict.keys():
                        if codon[:2] == k[:2] and codon[2] != k[2]:
                            snpRate = snpRateDict[codon[2] + k[2]]
                            aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                            del snpRate
                        else:
                            continue
                    del allrate
                else:
                    print("wrong in codonPos !!!!!")

                # print(aminoDict)
                sxx = 0
                nsxx = 0
                nopxx = 0
                opxx = 0
                normalxx = 0
                for iii in aminoDict.keys():
                    if aminoDict[iii][0] == condonBiasDict[codon][0]:
                        # print(aminoDict[iii][0])
                        sxx += 1 * aminoDict[iii][1]
                        if yi == 1:
                            if baseDict[iii[0]] != codon[0]:  # 由于 AT， CG互换的比例接近，所以不视为发生了两种密码子的更换
                                if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
                                    opxx += 1 * aminoDict[iii][1]
                                elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
                                    nopxx += 1 * aminoDict[iii][1]
                                else:
                                    print("wrong in bias !!!!!!!")
                            else:
                                normalxx += 1 * aminoDict[iii][1]
                        elif yi == 2:
                            if baseDict[iii[1]] != codon[1]:
                                if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
                                    opxx += 1 * aminoDict[iii][1]
                                elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
                                    nopxx += 1 * aminoDict[iii][1]
                                else:
                                    print("wrong in bias !!!!!!!")
                            else:
                                normalxx += 1 * aminoDict[iii][1]
                        elif yi == 3:
                            if baseDict[iii[2]] != codon[2]:
                                if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
                                    opxx += 1 * aminoDict[iii][1]
                                elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
                                    nopxx += 1 * aminoDict[iii][1]
                                else:
                                    print("wrong in bias !!!!!!!")
                            else:
                                normalxx += 1 * aminoDict[iii][1]
                        else:
                            print("wrong in codonpos !!!!!!!!")


                    else:
                        nsxx += 1 * aminoDict[iii][1]

                alloptimalsiteNum += opxx
                allnonoptimalsiteNum += nopxx
                allSsiteNum += sxx
                allNSsiteNum += nsxx
                allnormalsiteNum += normalxx

                del yi
                del opxx
                del nopxx
                del sxx
                del nsxx
                del normalxx
                del aminoDict

        # del codon

    # print("all gene: ")
    # print(allCDSLength)
    # print(allNSsiteNum)
    # print(allSsiteNum)
    # print(alloptimalsiteNum)
    # print(allnonoptimalsiteNum)
    # print(allnormalsiteNum)
    # print(alloptimalsiteNum+allnormalsiteNum+allnonoptimalsiteNum)
    # print("#############################")

    #######################################################################################################################

    ##############################################################################################

    ####### 全基因组上全部CDS的各种SNP ###############################################################
    allrealALLNumDict = {}
    allrealNsNumDict = {}
    allrealSNumDict = {}
    allrealOptimalNumDict = {}
    allrealNonoptimalNumDict = {}
    for ss in strainlist:
        allrealALLNumDict[ss] = {}
        allrealNsNumDict[ss] = 0
        allrealSNumDict[ss] = 0
        allrealOptimalNumDict[ss] = 0
        allrealNonoptimalNumDict[ss] = 0
        allrealAllsnpNum = 0
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
                codonpos = 0
                for g in allgeneDict.keys():
                    if int(allgeneDict[g].split("_")[0]) <= int(p.split("_")[0]) <= int(
                            allgeneDict[g].split("_")[1]):
                        if g not in allrealALLNumDict[ss].keys():
                            allrealALLNumDict[ss][g] = 1
                        else:
                            allrealALLNumDict[ss][g] += 1
                        if allgeneDict[g].split("_")[2] == "f":
                            if (int(p.split("_")[0]) - int(allgeneDict[g].split("_")[0]) + 1) % 3 == 1:
                                codon = p.split("_")[2] + allLongNewRefseqdict[int(p.split("_")[0]) + 1] + \
                                        allLongNewRefseqdict[int(p.split("_")[0]) + 2]
                                codonpos = 1
                                refcodon = p.split("_")[1] + allLongNewRefseqdict[int(p.split("_")[0]) + 1] + \
                                           allLongNewRefseqdict[int(p.split("_")[0]) + 2]

                            elif (int(p.split("_")[0]) - int(allgeneDict[g].split("_")[0]) + 1) % 3 == 2:
                                codon = allLongNewRefseqdict[int(p.split("_")[0]) - 1] + p.split("_")[2] + \
                                        allLongNewRefseqdict[int(p.split("_")[0]) + 1]
                                codonpos = 2
                                refcodon = allLongNewRefseqdict[int(p.split("_")[0]) - 1] + p.split("_")[1] + \
                                           allLongNewRefseqdict[int(p.split("_")[0]) + 1]

                            else:
                                codon = allLongNewRefseqdict[int(p.split("_")[0]) - 2] + allLongNewRefseqdict[
                                    int(p.split("_")[0]) - 1] + p.split("_")[2]
                                codonpos = 3
                                refcodon = allLongNewRefseqdict[int(p.split("_")[0]) - 2] + allLongNewRefseqdict[
                                    int(p.split("_")[0]) - 1] + p.split("_")[1]

                        else:
                            if (int(allgeneDict[g].split("_")[1]) - int(p.split("_")[0]) + 1) % 3 == 1:
                                codon = baseDict[p.split("_")[2]] + anti_allLongNewRefseqdict[
                                    int(p.split("_")[0]) - 1] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 2]
                                codonpos = 1
                                refcodon = baseDict[p.split("_")[1]] + anti_allLongNewRefseqdict[
                                    int(p.split("_")[0]) - 1] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 2]

                            elif (int(allgeneDict[g].split("_")[1]) - int(p.split("_")[0]) + 1) % 3 == 2:
                                codon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1] + baseDict[
                                    p.split("_")[2]] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 1]
                                codonpos = 2
                                refcodon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1] + baseDict[
                                    p.split("_")[1]] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 1]

                            else:
                                codon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 2] + \
                                        anti_allLongNewRefseqdict[
                                            int(p.split("_")[0]) + 1] + baseDict[p.split("_")[2]]
                                codonpos = 3
                                refcodon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 2] + \
                                           anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1] + baseDict[
                                               p.split("_")[1]]

                        break

                if codonpos == 0:
                    continue
                else:
                    allrealAllsnpNum += 1
                    # print(allrealAllsnpNum)
                    if condonBiasDict[codon][0] != condonBiasDict[refcodon][0]:
                        allrealNsNumDict[ss] += 1
                    else:
                        allrealSNumDict[ss] += 1
                        if codonpos == 1:
                            if baseDict[codon[0]] != refcodon[0]:
                                if float(condonBiasDict[codon][1]) > float(condonBiasDict[refcodon][1]):
                                    allrealOptimalNumDict[ss] += 1
                                elif float(condonBiasDict[codon][1]) < float(condonBiasDict[refcodon][1]):
                                    allrealNonoptimalNumDict[ss] += 1
                                else:
                                    print("wrong in bias !!!!!!!")
                        elif codonpos == 2:
                            if baseDict[codon[1]] != refcodon[1]:
                                if float(condonBiasDict[codon][1]) > float(condonBiasDict[refcodon][1]):
                                    allrealOptimalNumDict[ss] += 1
                                elif float(condonBiasDict[codon][1]) < float(condonBiasDict[refcodon][1]):
                                    allrealNonoptimalNumDict[ss] += 1
                                else:
                                    print("wrong in bias !!!!!!!")
                        elif codonpos == 3:
                            if baseDict[codon[2]] != refcodon[2]:
                                if float(condonBiasDict[codon][1]) > float(condonBiasDict[refcodon][1]):
                                    allrealOptimalNumDict[ss] += 1
                                elif float(condonBiasDict[codon][1]) < float(condonBiasDict[refcodon][1]):
                                    allrealNonoptimalNumDict[ss] += 1
                                else:
                                    print("wrong in bias !!!!!!!")
                        else:
                            print("wrong in codonpos !!!!!!!!")

        # del codon
        # del refcodon

    allrealSumNs = 0
    allrealSumS = 0
    allrealSumNop = 0
    allrealSumOp = 0
    allrealSumNSnum = 0
    allrealSumSnum = 0
    allrealSumNopnum = 0
    allrealSumOpnum = 0
    for allreal in allrealNsNumDict.keys():
        allrealSumNs += float(allrealNsNumDict[allreal]) / allNSsiteNum
        allrealSumNSnum += float(allrealNsNumDict[allreal])
    for allreal in allrealSNumDict.keys():
        allrealSumS += float(allrealSNumDict[allreal]) / allSsiteNum
        allrealSumSnum += float(allrealSNumDict[allreal])
    for allreal in allrealNonoptimalNumDict.keys():
        allrealSumNop += float(allrealNonoptimalNumDict[allreal]) / allnonoptimalsiteNum
        allrealSumNopnum += float(allrealNonoptimalNumDict[allreal])
    for allreal in allrealOptimalNumDict.keys():
        allrealSumOp += float(allrealOptimalNumDict[allreal]) / alloptimalsiteNum
        allrealSumOpnum += float(allrealOptimalNumDict[allreal])

    # print("########################")
    # print(cladeName)
    # print(allrealSumNSnum)
    # print(allrealSumSnum)
    # print(allrealSumNopnum)
    # print(allrealSumOpnum)
    # print("#########################")

    for epitype in epiTypeList:
            print(epitype+" starts ########")

            inputepitopefile = r"...\script_and_testData\testData\\"+epitype+".txt"

            outputFolder=r"...\script_and_testData\testData\output\calculation_of_selective_pressure\epitopes\\epitopes_cycle"+str(cycleNum)+"_"+epitype+"\\"
            outputFolder2 =r"...\script_and_testData\testData\output\calculation_of_selective_pressure\epitopes\\non-epitope-regions_cycle"+str(cycleNum)+"_"+epitype+"\\"


            outputfile9 = outputFolder+ cladeName + "_offsetDgree_ns.txt"
            outputfile10 = outputFolder+ cladeName + "_offsetDgree_s.txt"
            outputfile910 = outputFolder+ cladeName + "_offsetDgree_ns_dividedbyS.txt"
            outputfile11 = outputFolder+ cladeName + "_offsetDgree_nop_dividedbyS.txt"
            outputfile12 = outputFolder+ cladeName + "_offsetDgree_op_dividedbyS.txt"
            outputfile1112 = outputFolder+ cladeName + "_offsetDgree_nop_dividedbyOp.txt"

            outputfile13 = outputFolder + cladeName + "_offsetScore_ns.txt"
            outputfile14 = outputFolder+ cladeName + "_offsetScore_s.txt"
            outputfile15 = outputFolder+ cladeName + "_offsetScore_nop_dividedbyS.txt"
            outputfile16 = outputFolder+ cladeName + "_offsetScore_op_dividedbyS.txt"
            outputfile17 = outputFolder + cladeName + "_offsetScore_nop_dividedbyop.txt"
            outputfile18 = outputFolder + cladeName + "_offsetScore_ns_dividedbyS.txt"


            if not os.path.exists(outputFolder):
                os.mkdir(outputFolder)


            #################

            outputfile29 = outputFolder2 + cladeName + "_offsetDgree_ns.txt"
            outputfile210 = outputFolder2 + cladeName + "_offsetDgree_s.txt"
            outputfile2910 = outputFolder2 + cladeName + "_offsetDgree_ns_dividedbyS.txt"
            outputfile211 = outputFolder2 + cladeName + "_offsetDgree_nop_dividedbyS.txt"
            outputfile212 = outputFolder2 + cladeName + "_offsetDgree_op_dividedbyS.txt"
            outputfile21112 = outputFolder2 + cladeName + "_offsetDgree_nop_dividedbyOp.txt"

            outputfile213 = outputFolder2 + cladeName + "_offsetScore_ns.txt"
            outputfile214 = outputFolder2 + cladeName + "_offsetScore_s.txt"
            outputfile215 = outputFolder2 + cladeName + "_offsetScore_nop_dividedbyS.txt"
            outputfile216 = outputFolder2 + cladeName + "_offsetScore_op_dividedbyS.txt"
            outputfile217 = outputFolder2 + cladeName + "_offsetScore_nop_dividedbyop.txt"
            outputfile218 = outputFolder2 + cladeName + "_offsetScore_ns_dividedbyS.txt"

            if not os.path.exists(outputFolder2):
                os.mkdir(outputFolder2)

            #########################################################################################################


            epitopeDict = {}
            geneDict={}
            with open(inputepitopefile, "r") as input:
                for ll in input:
                    llx = ll.strip().split("\t")
                    if llx[0] != "Gene Pos-name":
                        if llx[0][-1] == "c":
                            epitopeDict[llx[0] + "_" + llx[3]] = str(int(llx[4]) * 3 - 2 + int(llx[1]) - 1) + "_" + \
                                                                    str(int(llx[5]) * 3 + int(llx[1]) - 1) + "_" + "c"

                            geneDict[llx[0]] = llx[1] +"_"+llx[2]+ "_" +"c" # 利用字典的key不能重复的特点自动去重
                        else:
                            epitopeDict[llx[0] + "_" + llx[3]] = str(int(llx[4]) * 3 - 2 + int(llx[1]) - 1) + "_" + \
                                                                    str(int(llx[5]) * 3 + int(llx[1]) - 1) + "_" + "f"

                            geneDict[llx[0]] = llx[1] +"_"+llx[2]+ "_" +"f"






            #######################################################################################################################
            ########## 全部抗原表位上的全部待选ns，s位点数目与待选optimal,nonoptimal位点数目 ############################################
            print("##### start all epitope candidate ns s et al. ########## ")
            epitopeNSsiteNum=0
            epitopeSsiteNum=0
            epitopeoptimalsiteNum=0
            epitopenonoptimalsiteNum=0
            epitopenormalsiteNum=0
            allEpitopeLength=0
            for gg in epitopeDict.keys():
                for pp in range(int(epitopeDict[gg].split("_")[0]), int(epitopeDict[gg].split("_")[1]) + 1, 3):
                    allEpitopeLength+=3
                    if epitopeDict[gg].split("_")[2] =="f":
                        codon = allLongNewRefseqdict[pp]+allLongNewRefseqdict[pp+1]+allLongNewRefseqdict[pp+2]
                    elif epitopeDict[gg].split("_")[2] =="c":
                        codon = anti_allLongNewRefseqdict[pp + 2] + anti_allLongNewRefseqdict[pp + 1] + anti_allLongNewRefseqdict[pp]
                    else:
                        print("WRONG IN epitopeDict!!!!!!!")

                    for yi in range(1, 4):
                        aminoDict = {}
                        if yi == 1:
                            allrate = allRateDict[codon[0]]
                            for k in condonBiasDict.keys():
                                if codon[1:] == k[1:] and codon[0] != k[0]:
                                    snpRate = snpRateDict[codon[0] + k[0]]
                                    aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                    del snpRate
                                else:
                                    continue
                            del allrate
                        elif yi == 2:
                            allrate = allRateDict[codon[1]]
                            for k in condonBiasDict.keys():
                                if codon[0] == k[0] and codon[-1] == k[-1] and codon[1] != k[1]:
                                    snpRate = snpRateDict[codon[1] + k[1]]
                                    aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                    del snpRate
                                else:
                                    continue
                            del allrate
                        elif yi == 3:
                            allrate = allRateDict[codon[2]]
                            for k in condonBiasDict.keys():
                                if codon[:2] == k[:2] and codon[2] != k[2]:
                                    snpRate = snpRateDict[codon[2] + k[2]]
                                    aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                    del snpRate
                                else:
                                    continue
                            del allrate
                        else:
                            print("wrong in codonPos !!!!!")

                        # print(aminoDict)
                        sxx = 0
                        nsxx = 0
                        nopxx = 0
                        opxx = 0
                        normalxx = 0
                        for iii in aminoDict.keys():
                            if aminoDict[iii][0] == condonBiasDict[codon][0]:
                                # print(aminoDict[iii][0])
                                sxx += 1 * aminoDict[iii][1]
                                if yi == 1:
                                    if baseDict[iii[0]] != codon[0]:  # 由于 AT， CG互换的比例接近，所以不视为发生了两种密码子的更换
                                        if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
                                            opxx += 1 * aminoDict[iii][1]
                                        elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
                                            nopxx += 1 * aminoDict[iii][1]
                                        else:
                                            print("wrong in bias !!!!!!!")
                                    else:
                                        normalxx += 1 * aminoDict[iii][1]
                                elif yi == 2:
                                    if baseDict[iii[1]] != codon[1]:
                                        if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
                                            opxx += 1 * aminoDict[iii][1]
                                        elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
                                            nopxx += 1 * aminoDict[iii][1]
                                        else:
                                            print("wrong in bias !!!!!!!")
                                    else:
                                        normalxx += 1 * aminoDict[iii][1]
                                elif yi == 3:
                                    if baseDict[iii[2]] != codon[2]:
                                        if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
                                            opxx += 1 * aminoDict[iii][1]
                                        elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
                                            nopxx += 1 * aminoDict[iii][1]
                                        else:
                                            print("wrong in bias !!!!!!!")
                                    else:
                                        normalxx += 1 * aminoDict[iii][1]
                                else:
                                    print("wrong in codonpos !!!!!!!!")


                            else:
                                nsxx += 1 * aminoDict[iii][1]

                        epitopeoptimalsiteNum += opxx
                        epitopenonoptimalsiteNum += nopxx
                        epitopeSsiteNum += sxx
                        epitopeNSsiteNum += nsxx
                        epitopenormalsiteNum += normalxx

                        del yi
                        del opxx
                        del nopxx
                        del sxx
                        del nsxx
                        del normalxx
                        del aminoDict

                # del codon



            # print(allEpitopeLength)
            # print(epitopeNSsiteNum)
            # print(epitopeSsiteNum)
            # print(epitopeoptimalsiteNum)
            # print(epitopenonoptimalsiteNum)
            # print(epitopenormalsiteNum)
            # print(epitopeoptimalsiteNum+epitopenormalsiteNum+epitopenonoptimalsiteNum)

            ###########################################################################################




            ##############################################################################################

            #######################################################################################################################
            ########## 全部抗原表位所在基因上的全部待选ns，s位点数目与待选optimal,nonoptimal位点数目 ############################################
            print("##### start genes with epitope candidate ns s et al. ########## ")
            geneNSsiteNum=0
            geneSsiteNum=0
            geneoptimalsiteNum=0
            genenonoptimalsiteNum=0
            genenormalsiteNum=0
            allgeneLength=0
            for gg in geneDict.keys():
                for pp in range(int(geneDict[gg].split("_")[0]), int(geneDict[gg].split("_")[1]) + 1, 3):
                    allgeneLength+=3
                    if geneDict[gg].split("_")[2] =="f":
                        codon = allLongNewRefseqdict[pp]+allLongNewRefseqdict[pp+1]+allLongNewRefseqdict[pp+2]
                    elif geneDict[gg].split("_")[2] =="c":
                        codon = anti_allLongNewRefseqdict[pp + 2] + anti_allLongNewRefseqdict[pp + 1] + anti_allLongNewRefseqdict[pp]
                    else:
                        print("WRONG IN geneDict!!!!!!!")

                    for yi in range(1, 4):
                        aminoDict = {}
                        if yi == 1:
                            allrate = allRateDict[codon[0]]
                            for k in condonBiasDict.keys():
                                if codon[1:] == k[1:] and codon[0] != k[0]:
                                    snpRate = snpRateDict[codon[0] + k[0]]
                                    aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                    del snpRate
                                else:
                                    continue
                            del allrate
                        elif yi == 2:
                            allrate = allRateDict[codon[1]]
                            for k in condonBiasDict.keys():
                                if codon[0] == k[0] and codon[-1] == k[-1] and codon[1] != k[1]:
                                    snpRate = snpRateDict[codon[1] + k[1]]
                                    aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                    del snpRate
                                else:
                                    continue
                            del allrate
                        elif yi == 3:
                            allrate = allRateDict[codon[2]]
                            for k in condonBiasDict.keys():
                                if codon[:2] == k[:2] and codon[2] != k[2]:
                                    snpRate = snpRateDict[codon[2] + k[2]]
                                    aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                    del snpRate
                                else:
                                    continue
                            del allrate
                        else:
                            print("wrong in codonPos !!!!!")

                        # print(aminoDict)
                        sxx = 0
                        nsxx = 0
                        nopxx = 0
                        opxx = 0
                        normalxx = 0
                        for iii in aminoDict.keys():
                            if aminoDict[iii][0] == condonBiasDict[codon][0]:
                                # print(aminoDict[iii][0])
                                sxx += 1 * aminoDict[iii][1]
                                if yi == 1:
                                    if baseDict[iii[0]] != codon[0]:  # 由于 AT， CG互换的比例接近，所以不视为发生了两种密码子的更换
                                        if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
                                            opxx += 1 * aminoDict[iii][1]
                                        elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
                                            nopxx += 1 * aminoDict[iii][1]
                                        else:
                                            print("wrong in bias !!!!!!!")
                                    else:
                                        normalxx += 1 * aminoDict[iii][1]
                                elif yi == 2:
                                    if baseDict[iii[1]] != codon[1]:
                                        if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
                                            opxx += 1 * aminoDict[iii][1]
                                        elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
                                            nopxx += 1 * aminoDict[iii][1]
                                        else:
                                            print("wrong in bias !!!!!!!")
                                    else:
                                        normalxx += 1 * aminoDict[iii][1]
                                elif yi == 3:
                                    if baseDict[iii[2]] != codon[2]:
                                        if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
                                            opxx += 1 * aminoDict[iii][1]
                                        elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
                                            nopxx += 1 * aminoDict[iii][1]
                                        else:
                                            print("wrong in bias !!!!!!!")
                                    else:
                                        normalxx += 1 * aminoDict[iii][1]
                                else:
                                    print("wrong in codonpos !!!!!!!!")


                            else:
                                nsxx += 1 * aminoDict[iii][1]

                        geneoptimalsiteNum += opxx
                        genenonoptimalsiteNum += nopxx
                        geneSsiteNum += sxx
                        geneNSsiteNum += nsxx
                        genenormalsiteNum += normalxx

                        del yi
                        del opxx
                        del nopxx
                        del sxx
                        del nsxx
                        del normalxx
                        del aminoDict

                # del codon



            # print(allgeneLength)
            # print(geneNSsiteNum)
            # print(geneSsiteNum)
            # print(geneoptimalsiteNum)
            # print(genenonoptimalsiteNum)
            # print(genenormalsiteNum)
            # print(geneoptimalsiteNum+genenormalsiteNum+genenonoptimalsiteNum)

            ###########################################################################################



            ####################################################################################################

            ############################################################################################
            # 表位上的真实突变
            realALLNumDict = {}
            realNsNumDict = {}
            realSNumDict = {}
            realOptimalNumDict = {}
            realNonoptimalNumDict = {}
            for ss in strainlist:
                snpnum=0
                realALLNumDict[ss] = {}
                realNsNumDict[ss] = 0
                realSNumDict[ss] = 0
                realOptimalNumDict[ss] = 0
                realNonoptimalNumDict[ss] = 0
                realAllsnpNum=0
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
                        codonpos = 0
                        for g in epitopeDict.keys():
                            if g not in realALLNumDict[ss].keys():
                                realALLNumDict[ss][g]=0
                            if int(epitopeDict[g].split("_")[0]) <= int(p.split("_")[0]) <= int(epitopeDict[g].split("_")[1]):
                                snpnum+=1
                                realALLNumDict[ss][g] += 1
                                if epitopeDict[g].split("_")[2] == "f":
                                    if (int(p.split("_")[0]) - int(epitopeDict[g].split("_")[0]) + 1) % 3 == 1:
                                        codon = p.split("_")[2] + allLongNewRefseqdict[int(p.split("_")[0]) + 1] + \
                                                allLongNewRefseqdict[int(p.split("_")[0]) + 2]
                                        codonpos = 1
                                        refcodon = p.split("_")[1] + allLongNewRefseqdict[int(p.split("_")[0]) + 1] + \
                                                   allLongNewRefseqdict[int(p.split("_")[0]) + 2]

                                    elif (int(p.split("_")[0]) - int(epitopeDict[g].split("_")[0]) + 1) % 3 == 2:
                                        codon = allLongNewRefseqdict[int(p.split("_")[0]) - 1] + p.split("_")[2] + \
                                                allLongNewRefseqdict[int(p.split("_")[0]) + 1]
                                        codonpos = 2
                                        refcodon = allLongNewRefseqdict[int(p.split("_")[0]) - 1] + p.split("_")[1] + \
                                                   allLongNewRefseqdict[int(p.split("_")[0]) + 1]

                                    else:
                                        codon = allLongNewRefseqdict[int(p.split("_")[0]) - 2] + allLongNewRefseqdict[
                                            int(p.split("_")[0]) - 1] + p.split("_")[2]
                                        codonpos = 3
                                        refcodon = allLongNewRefseqdict[int(p.split("_")[0]) - 2] + allLongNewRefseqdict[
                                            int(p.split("_")[0]) - 1] + p.split("_")[1]

                                else:
                                    if (int(epitopeDict[g].split("_")[1]) - int(p.split("_")[0]) + 1) % 3 == 1:
                                        codon = baseDict[p.split("_")[2]] + anti_allLongNewRefseqdict[
                                            int(p.split("_")[0]) - 1] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 2]
                                        codonpos = 1
                                        refcodon = baseDict[p.split("_")[1]] + anti_allLongNewRefseqdict[
                                            int(p.split("_")[0]) - 1] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 2]

                                    elif (int(epitopeDict[g].split("_")[1]) - int(p.split("_")[0]) + 1) % 3 == 2:
                                        codon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1] + baseDict[
                                            p.split("_")[2]] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 1]
                                        codonpos = 2
                                        refcodon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1] + baseDict[
                                            p.split("_")[1]] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 1]

                                    else:
                                        codon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 2] + anti_allLongNewRefseqdict[
                                            int(p.split("_")[0]) + 1] + baseDict[p.split("_")[2]]
                                        codonpos = 3
                                        refcodon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 2] + \
                                                   anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1] + baseDict[
                                                       p.split("_")[1]]

                                break

                        if codonpos == 0:
                            continue
                        else:
                            realAllsnpNum += 1
                            # print(realAllsnpNum)
                            if condonBiasDict[codon][0] != condonBiasDict[refcodon][0]:
                                realNsNumDict[ss] += 1
                            else:
                                realSNumDict[ss] += 1
                                if codonpos == 1:
                                    if baseDict[codon[0]] != refcodon[0]:
                                        if float(condonBiasDict[codon][1]) > float(condonBiasDict[refcodon][1]):
                                            realOptimalNumDict[ss] += 1
                                        elif float(condonBiasDict[codon][1]) < float(condonBiasDict[refcodon][1]):
                                            realNonoptimalNumDict[ss] += 1
                                        else:
                                            print("wrong in bias !!!!!!!")
                                elif codonpos == 2:
                                    if baseDict[codon[1]] != refcodon[1]:
                                        if float(condonBiasDict[codon][1]) > float(condonBiasDict[refcodon][1]):
                                            realOptimalNumDict[ss] += 1
                                        elif float(condonBiasDict[codon][1]) < float(condonBiasDict[refcodon][1]):
                                            realNonoptimalNumDict[ss] += 1
                                        else:
                                            print("wrong in bias !!!!!!!")
                                elif codonpos == 3:
                                    if baseDict[codon[2]] != refcodon[2]:
                                        if float(condonBiasDict[codon][1]) > float(condonBiasDict[refcodon][1]):
                                            realOptimalNumDict[ss] += 1
                                        elif float(condonBiasDict[codon][1]) < float(condonBiasDict[refcodon][1]):
                                            realNonoptimalNumDict[ss] += 1
                                        else:
                                            print("wrong in bias !!!!!!!")
                                else:
                                    print("wrong in codonpos !!!!!!!!")

                # del codon
                # del refcodon

                if snpnum==0:
                    print("sample NO snp in epitope: " + ss)




            realSumNs = 1E-6  # 加个极小的值，避免分母为0
            realSumS = 1E-6
            realSumNop = 1E-6
            realSumOp = 1E-6
            realSumNSnum = 0.1
            realSumSnum = 0.1
            realSumNopnum = 0.1
            realSumOpnum = 0.1
            for real in realNsNumDict.keys():
                realSumNs += float(realNsNumDict[real] / epitopeNSsiteNum)
                realSumNSnum += float(realNsNumDict[real])
            for real in realSNumDict.keys():
                realSumS += float(realSNumDict[real] / epitopeSsiteNum)
                realSumSnum += float(realSNumDict[real])
            for real in realNonoptimalNumDict.keys():
                realSumNop += float(realNonoptimalNumDict[real] / epitopenonoptimalsiteNum)
                realSumNopnum += float(realNonoptimalNumDict[real])
            for real in realOptimalNumDict.keys():
                realSumOp += float(realOptimalNumDict[real] / epitopeoptimalsiteNum)
                realSumOpnum += float(realOptimalNumDict[real])


            epitoperealSumNs = realSumNs
            epitoperealSumS = realSumS
            epitoperealSumNop = realSumNop
            epitoperealSumOp = realSumOp
            epitoperealSumNSnum = realSumNSnum
            epitoperealSumSnum = realSumSnum
            epitoperealSumNopnum=  realSumNopnum
            epitoperealSumOpnum= realSumOpnum

            epitoperealALLNumDict = realALLNumDict
            epitoperealNsNumDict = realNsNumDict
            epitoperealSNumDict = realSNumDict
            epitoperealOptimalNumDict = realOptimalNumDict
            epitoperealNonoptimalNumDict = realNonoptimalNumDict
            # print(epitoperealALLNumDict)

            del realSumNs
            del realSumS
            del realSumNop
            del realSumOp
            del realSumNSnum
            del realSumSnum
            del realSumNopnum
            del realSumOpnum
            del realALLNumDict
            del realNsNumDict
            del realSNumDict
            del realOptimalNumDict
            del realNonoptimalNumDict



            epitoperealSumNSS = (epitoperealSumNSnum + epitoperealSumSnum) / (
                        (epitopeNSsiteNum + epitopeSsiteNum) * len(strainlist))
            epitoperealSumnopop = (epitoperealSumOpnum + epitoperealSumNopnum) / (
                        (epitopenonoptimalsiteNum + epitopeoptimalsiteNum) * len(strainlist))

            allnssSNPrate = (allrealSumNSnum + allrealSumSnum) / ((allNSsiteNum + allSsiteNum) * len(strainlist))
            allnopopSNPrate = (allrealSumOpnum + allrealSumNopnum) / (
                        (allnonoptimalsiteNum + alloptimalsiteNum) * len(strainlist))

            epitoperALL = epitoperealSumNSS / allnssSNPrate
            epitoperopnop = epitoperealSumnopop / allnopopSNPrate

            print("epitoperALL: " + str(epitoperALL))
            print("epitoperopnop: " + str(epitoperopnop))

            ############################################################################################
            # 表位所在基因非表位区域上的真实突变
            realALLNumDict = {}
            realNsNumDict = {}
            realSNumDict = {}
            realOptimalNumDict = {}
            realNonoptimalNumDict = {}
            for ss in strainlist:
                realALLNumDict[ss] = {}
                realNsNumDict[ss] = 0
                realSNumDict[ss] = 0
                realOptimalNumDict[ss] = 0
                realNonoptimalNumDict[ss] = 0
                realAllsnpNum = 0
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
                        codonpos = 0
                        for g in geneDict.keys():
                            if g not in realALLNumDict[ss].keys():
                                realALLNumDict[ss][g] =0
                            if int(geneDict[g].split("_")[0]) <= int(p.split("_")[0]) <= int(geneDict[g].split("_")[1]):
                                realALLNumDict[ss][g] += 1
                                if geneDict[g].split("_")[2] == "f":
                                    if (int(p.split("_")[0]) - int(geneDict[g].split("_")[0]) + 1) % 3 == 1:
                                        codon = p.split("_")[2] + allLongNewRefseqdict[int(p.split("_")[0]) + 1] + \
                                                allLongNewRefseqdict[int(p.split("_")[0]) + 2]
                                        codonpos = 1
                                        refcodon = p.split("_")[1] + allLongNewRefseqdict[int(p.split("_")[0]) + 1] + \
                                                   allLongNewRefseqdict[int(p.split("_")[0]) + 2]

                                    elif (int(p.split("_")[0]) - int(geneDict[g].split("_")[0]) + 1) % 3 == 2:
                                        codon = allLongNewRefseqdict[int(p.split("_")[0]) - 1] + p.split("_")[2] + \
                                                allLongNewRefseqdict[int(p.split("_")[0]) + 1]
                                        codonpos = 2
                                        refcodon = allLongNewRefseqdict[int(p.split("_")[0]) - 1] + p.split("_")[1] + \
                                                   allLongNewRefseqdict[int(p.split("_")[0]) + 1]

                                    else:
                                        codon = allLongNewRefseqdict[int(p.split("_")[0]) - 2] + allLongNewRefseqdict[
                                            int(p.split("_")[0]) - 1] + p.split("_")[2]
                                        codonpos = 3
                                        refcodon = allLongNewRefseqdict[int(p.split("_")[0]) - 2] + allLongNewRefseqdict[
                                            int(p.split("_")[0]) - 1] + p.split("_")[1]

                                else:
                                    if (int(geneDict[g].split("_")[1]) - int(p.split("_")[0]) + 1) % 3 == 1:
                                        codon = baseDict[p.split("_")[2]] + anti_allLongNewRefseqdict[
                                            int(p.split("_")[0]) - 1] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 2]
                                        codonpos = 1
                                        refcodon = baseDict[p.split("_")[1]] + anti_allLongNewRefseqdict[
                                            int(p.split("_")[0]) - 1] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 2]

                                    elif (int(geneDict[g].split("_")[1]) - int(p.split("_")[0]) + 1) % 3 == 2:
                                        codon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1] + baseDict[
                                            p.split("_")[2]] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 1]
                                        codonpos = 2
                                        refcodon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1] + baseDict[
                                            p.split("_")[1]] + anti_allLongNewRefseqdict[int(p.split("_")[0]) - 1]

                                    else:
                                        codon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 2] + \
                                                anti_allLongNewRefseqdict[
                                                    int(p.split("_")[0]) + 1] + baseDict[p.split("_")[2]]
                                        codonpos = 3
                                        refcodon = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 2] + \
                                                   anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1] + baseDict[
                                                       p.split("_")[1]]

                                break

                        if codonpos == 0:
                            continue
                        else:
                            realAllsnpNum += 1
                            # print(realAllsnpNum)
                            if condonBiasDict[codon][0] != condonBiasDict[refcodon][0]:
                                realNsNumDict[ss] += 1
                            else:
                                realSNumDict[ss] += 1
                                if codonpos == 1:
                                    if baseDict[codon[0]] != refcodon[0]:
                                        if float(condonBiasDict[codon][1]) > float(condonBiasDict[refcodon][1]):
                                            realOptimalNumDict[ss] += 1
                                        elif float(condonBiasDict[codon][1]) < float(condonBiasDict[refcodon][1]):
                                            realNonoptimalNumDict[ss] += 1
                                        else:
                                            print("wrong in bias !!!!!!!")
                                elif codonpos == 2:
                                    if baseDict[codon[1]] != refcodon[1]:
                                        if float(condonBiasDict[codon][1]) > float(condonBiasDict[refcodon][1]):
                                            realOptimalNumDict[ss] += 1
                                        elif float(condonBiasDict[codon][1]) < float(condonBiasDict[refcodon][1]):
                                            realNonoptimalNumDict[ss] += 1
                                        else:
                                            print("wrong in bias !!!!!!!")
                                elif codonpos == 3:
                                    if baseDict[codon[2]] != refcodon[2]:
                                        if float(condonBiasDict[codon][1]) > float(condonBiasDict[refcodon][1]):
                                            realOptimalNumDict[ss] += 1
                                        elif float(condonBiasDict[codon][1]) < float(condonBiasDict[refcodon][1]):
                                            realNonoptimalNumDict[ss] += 1
                                        else:
                                            print("wrong in bias !!!!!!!")
                                else:
                                    print("wrong in codonpos !!!!!!!!")

                # del codon
                # del refcodon

            realSumNs = 1E-6  # 加个极小的值，避免分母为0
            realSumS = 1E-6
            realSumNop = 1E-6
            realSumOp = 1E-6
            realSumNSnum = 0.1
            realSumSnum = 0.1
            realSumNopnum = 0.1
            realSumOpnum = 0.1
            for real in realNsNumDict.keys():
                realSumNs += ( float(realNsNumDict[real]) - float(epitoperealNsNumDict[real]) ) / ( geneNSsiteNum -epitopeNSsiteNum )
                realSumNSnum += float(realNsNumDict[real])
            for real in realSNumDict.keys():
                realSumS +=  ( float(realSNumDict[real]) -  float(epitoperealSNumDict[real]) ) / ( geneSsiteNum -epitopeSsiteNum )
                realSumSnum += float(realSNumDict[real])
            for real in realNonoptimalNumDict.keys():
                realSumNop += ( float(realNonoptimalNumDict[real]) -  float(epitoperealNonoptimalNumDict[real]) )  / ( genenonoptimalsiteNum - epitopenonoptimalsiteNum)
                realSumNopnum += float(realNonoptimalNumDict[real])
            for real in realOptimalNumDict.keys():
                realSumOp += ( float(realOptimalNumDict[real]) -  float(epitoperealOptimalNumDict[real]) )  / ( geneoptimalsiteNum - epitopeoptimalsiteNum)
                realSumOpnum += float(realOptimalNumDict[real])

            geneOUTEPIrealSumNs = realSumNs
            geneOUTEPIrealSumS = realSumS
            geneOUTEPIrealSumNop = realSumNop
            geneOUTEPIrealSumOp = realSumOp
            geneOUTEPIrealSumNSnum = realSumNSnum - epitoperealSumNSnum
            geneOUTEPIrealSumSnum = realSumSnum - epitoperealSumSnum
            geneOUTEPIrealSumNopnum = realSumNopnum - epitoperealSumNopnum
            geneOUTEPIrealSumOpnum = realSumOpnum - epitoperealSumOpnum

            generealALLNumDict = realALLNumDict

            geneOUTEPIrealALLNumDict={}
            for skk in generealALLNumDict.keys():
                geneOUTEPIrealALLNumDict[skk]={}
                for ggg in generealALLNumDict[skk]:
                    geneOUTEPIrealALLNumDict[skk][ggg] = [generealALLNumDict[skk][ggg],[]]
                    for eee in epitoperealALLNumDict[skk]:
                        if eee.split("_")[0] != ggg:
                            continue
                        else:
                            geneOUTEPIrealALLNumDict[skk][ggg][0] = geneOUTEPIrealALLNumDict[skk][ggg][0] - epitoperealALLNumDict[skk][eee]
                            geneOUTEPIrealALLNumDict[skk][ggg][1].append(epitopeDict[eee])


            del realSumNs
            del realSumS
            del realSumNop
            del realSumOp
            del realSumNSnum
            del realSumSnum
            del realSumNopnum
            del realSumOpnum
            del realALLNumDict
            del realNsNumDict
            del realSNumDict
            del realOptimalNumDict
            del realNonoptimalNumDict


            geneOUTEPIrealSumNSS = (geneOUTEPIrealSumNSnum + geneOUTEPIrealSumSnum) / (
                        ((geneNSsiteNum - epitopeNSsiteNum) + (geneSsiteNum - epitopeSsiteNum)) * len(strainlist))
            geneOUTEPIrealSumnopop = (geneOUTEPIrealSumOpnum + geneOUTEPIrealSumNopnum) / (((
                                                                                                        genenonoptimalsiteNum - epitopenonoptimalsiteNum) + (
                                                                                                        geneoptimalsiteNum - epitopeoptimalsiteNum)) * len(
                strainlist))

            allnssSNPrate = (allrealSumNSnum + allrealSumSnum) / ((allNSsiteNum + allSsiteNum) * len(strainlist))
            allnopopSNPrate = (allrealSumOpnum + allrealSumNopnum) / (
                        (allnonoptimalsiteNum + alloptimalsiteNum) * len(strainlist))

            geneOUTEPIrALL = geneOUTEPIrealSumNSS / allnssSNPrate
            geneOUTEPIropnop = geneOUTEPIrealSumnopop / allnopopSNPrate

            print("geneOUTEPIrALL: " + str(geneOUTEPIrALL))
            print("geneOUTEPIropnop: " + str(geneOUTEPIropnop))



            # print(geneOUTEPIrealALLNumDict)

            #####################  模拟  #####################################################################
            #################################################################################################

            print("####### binomial test is working")
            cycleNum = cycleNum  # 模拟的次数

            ######################### 抗原表位上的模拟 ###################################################
            ranSumPnsDict = {}
            ranSumPsDict = {}
            ranSumPnopDict = {}
            ranSumPopDict = {}

            print("epitope is starting ################")
            for cy in tqdm(range(1, cycleNum + 1),desc="抗原表位模拟"):
                # for cccc in range(0, cycleNum + 1, int(cycleNum/10)):
                #     if cy == cccc:
                #         print("epitope: "+str(cy) + "/" + str(cycleNum))

                # 随机在所有表位的所有位点上进行和真实的allSNPnum数量相同的突变
                moniNsRateDict = {}
                moniSRateDict = {}
                moniNORateDict = {}
                moniORateDict = {}

                samplen = 0
                for sample in epitoperealALLNumDict.keys():
                    samplen += 1
                    moniNsRateDict[sample] = []
                    moniSRateDict[sample] = []
                    moniNORateDict[sample] = []
                    moniORateDict[sample] = []

                    # 在观测到真实突变的表位上生成等数量随机位置的突变
                    rgn = 0
                    for rg in epitoperealALLNumDict[sample].keys():
                        rgn += 1
                        fposList = []
                        sposList = []
                        tposList = []
                        if epitoperealALLNumDict[sample][rg]==0: # 没有发生突变
                            continue
                        else:
                            # 将三个密码子位置分开,按照突变率分开抽取
                            if epitopeDict[rg].split("_")[2] == "f":
                                for yp in range(int(epitopeDict[rg].split("_")[0]), int(epitopeDict[rg].split("_")[1]) + 1):
                                    if yp % 3 == 1:
                                        fposList.append(yp)
                                    elif yp % 3 == 2:
                                        sposList.append(yp)
                                    else:
                                        tposList.append(yp)
                            else:
                                for yp in range(int(epitopeDict[rg].split("_")[1]), int(epitopeDict[rg].split("_")[0]) - 1,-1):
                                    if yp % 3 == 0:
                                        fposList.append(yp)
                                    elif yp % 3 == 2:
                                        sposList.append(yp)
                                    else:
                                        tposList.append(yp)

                            if len(fposList) != len(sposList) or len(fposList) != len(tposList) or len(sposList) != len(tposList):
                                print("wrong!!!!!!! in f&s&tposList !!!!!")
                                print(str(int(epitopeDict[rg].split("_")[0]))+"_"+str(int(epitopeDict[rg].split("_")[1])))

                            for pp in range(1, epitoperealALLNumDict[sample][rg] + 1):
                                np.random.seed(randomseed*31+19*cycleNum + cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 97)
                                choicePos1=np.random.binomial(1,codonPosRateDict["1st"],1)
                                if choicePos1 ==1:
                                    choicePOSx =1
                                    del choicePos1
                                else:
                                    np.random.seed(randomseed*31+19*cycleNum + cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 197)
                                    secondRate=codonPosRateDict["2nd"] / (codonPosRateDict["2nd"]+codonPosRateDict["3rd"])
                                    choicePos2=np.random.binomial(1,secondRate,1)
                                    if choicePos2 ==1:
                                        choicePOSx =2
                                        del choicePos2
                                    else:
                                        choicePOSx =3
                                        del choicePos2

                                if choicePOSx ==1:
                                    np.random.seed(randomseed*31+19*cycleNum +cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 317)
                                    ranPos = np.random.choice(fposList,1)[0]
                                elif choicePOSx ==2:
                                    np.random.seed(randomseed*31+19*cycleNum +cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 317)
                                    ranPos = np.random.choice(sposList, 1)[0]
                                elif choicePOSx ==3:
                                    np.random.seed(randomseed*31+19*cycleNum +cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 317)
                                    ranPos = np.random.choice(tposList, 1)[0]


                                refcodon = ""
                                codonpos = 0
                                if epitopeDict[rg].split("_")[2] == "f":
                                    if (ranPos - int(epitopeDict[rg].split("_")[0]) + 1) % 3 == 1:
                                        codonpos = 1
                                        refcodon = allLongNewRefseqdict[ranPos] + allLongNewRefseqdict[ranPos + 1] + \
                                                   allLongNewRefseqdict[ranPos + 2]

                                    elif (ranPos - int(epitopeDict[rg].split("_")[0]) + 1) % 3 == 2:
                                        codonpos = 2
                                        refcodon = allLongNewRefseqdict[ranPos - 1] + allLongNewRefseqdict[ranPos] + \
                                                   allLongNewRefseqdict[ranPos + 1]

                                    else:
                                        codonpos = 3
                                        refcodon = allLongNewRefseqdict[ranPos - 2] + allLongNewRefseqdict[ranPos - 1] + \
                                                   allLongNewRefseqdict[ranPos]

                                else:
                                    if (int(epitopeDict[rg].split("_")[1]) - ranPos + 1) % 3 == 1:
                                        codonpos = 1
                                        refcodon = anti_allLongNewRefseqdict[ranPos] + anti_allLongNewRefseqdict[ranPos - 1] + \
                                                   anti_allLongNewRefseqdict[ranPos - 2]

                                    elif (int(epitopeDict[rg].split("_")[1]) - ranPos + 1) % 3 == 2:
                                        codonpos = 2
                                        refcodon = anti_allLongNewRefseqdict[ranPos + 1] + anti_allLongNewRefseqdict[ranPos] + \
                                                   anti_allLongNewRefseqdict[ranPos - 1]

                                    else:
                                        codonpos = 3
                                        refcodon = anti_allLongNewRefseqdict[ranPos + 2] + anti_allLongNewRefseqdict[
                                            ranPos + 1] + anti_allLongNewRefseqdict[ranPos]

                                if codonpos == 0:
                                    print("wrong in random codonpos !!!!!!")
                                else:

                                    allRateDict={}
                                    rateα = titvDict["ti"]
                                    rateβ1 = titvDict["tvtype1"]
                                    rateβ2 = titvDict["tvtype2"]
                                    rateΘ1 = tamuraRateDict["rate1"]
                                    rateΘ2 = tamuraRateDict["rate2"]
                                    allRateDict["A"] = 2 * rateβ1 * rateΘ2 + 2 * rateβ2 * rateΘ1 + rateα * rateΘ1
                                    allRateDict["T"] = 2 * rateβ1 * rateΘ2 + 2 * rateβ2 * rateΘ1 + rateα * rateΘ1
                                    allRateDict["C"] = 2 * rateβ1 * rateΘ1 + 2 * rateβ2 * rateΘ2 + rateα * rateΘ2
                                    allRateDict["G"] = 2 * rateβ1 * rateΘ1 + 2 * rateβ2 * rateΘ2 + rateα * rateΘ2
                                    snpRateDict={}
                                    snpRateDict["AT"] = 2 * rateβ1 * rateΘ2
                                    snpRateDict["AC"] = 2 * rateβ2 * rateΘ1
                                    snpRateDict["AG"] = rateα * rateΘ1
                                    snpRateDict["TA"] = 2 * rateβ1 * rateΘ2
                                    snpRateDict["TC"] = rateα * rateΘ1
                                    snpRateDict["TG"] = 2 * rateβ2 * rateΘ1
                                    snpRateDict["CA"] = 2 * rateβ2 * rateΘ2
                                    snpRateDict["CT"] = rateα * rateΘ2
                                    snpRateDict["CG"] = 2 * rateβ1 * rateΘ1
                                    snpRateDict["GA"] = rateα * rateΘ2
                                    snpRateDict["GT"] = 2 * rateβ2 * rateΘ2
                                    snpRateDict["GC"] = 2 * rateβ1 * rateΘ1


                                    aminoDict = {}
                                    if codonpos == 1:
                                        allrate = allRateDict[refcodon[0]]
                                        for k in condonBiasDict.keys():
                                            if refcodon[1:] == k[1:] and refcodon[0] != k[0]:
                                                snpRate=snpRateDict[refcodon[0] + k[0]]
                                                aminoDict[k] = [condonBiasDict[k][0], snpRate/allrate]
                                                del snpRate
                                            else:
                                                continue
                                        del allrate
                                    elif codonpos == 2:
                                        allrate = allRateDict[refcodon[1]]
                                        for k in condonBiasDict.keys():
                                            if refcodon[0] == k[0] and refcodon[-1] == k[-1] and refcodon[1] != k[1]:
                                                snpRate = snpRateDict[refcodon[1] + k[1]]
                                                aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                                del snpRate
                                            else:
                                                continue
                                        del allrate
                                    elif codonpos == 3:
                                        allrate = allRateDict[refcodon[2]]
                                        for k in condonBiasDict.keys():
                                            if refcodon[:2] == k[:2]  and refcodon[2] != k[2]:
                                                snpRate = snpRateDict[refcodon[2] + k[2]]
                                                aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                                del snpRate
                                            else:
                                                continue
                                        del allrate
                                    else:
                                        print("wrong in codonpos !!!!!")

                                    # print(aminoDict)
                                    sxx = 0
                                    nsxx = 0
                                    nopxx = 0
                                    opxx = 0
                                    for iii in aminoDict.keys():
                                        if aminoDict[iii][0] == condonBiasDict[refcodon][0]:
                                            # print(aminoDict[iii][0])
                                            sxx += 1*aminoDict[iii][1]
                                            if codonpos == 1:
                                                if baseDict[iii[0]] != refcodon[0]:
                                                    if float(condonBiasDict[iii][1]) > float(condonBiasDict[refcodon][1]):
                                                        opxx += 1*aminoDict[iii][1]
                                                    elif float(condonBiasDict[iii][1]) < float(condonBiasDict[refcodon][1]):
                                                        nopxx += 1*aminoDict[iii][1]
                                                    else:
                                                        print("wrong in bias !!!!!!!")
                                                else:
                                                    continue
                                            elif codonpos == 2:
                                                if baseDict[iii[1]] != refcodon[1]:
                                                    if float(condonBiasDict[iii][1]) > float(condonBiasDict[refcodon][1]):
                                                        opxx += 1 * aminoDict[iii][1]
                                                    elif float(condonBiasDict[iii][1]) < float(condonBiasDict[refcodon][1]):
                                                        nopxx += 1 * 1 * aminoDict[iii][1]
                                                    else:
                                                        print("wrong in bias !!!!!!!")
                                                else:
                                                    continue
                                            elif codonpos == 3:
                                                if baseDict[iii[2]] != refcodon[2]:
                                                    if float(condonBiasDict[iii][1]) > float(condonBiasDict[refcodon][1]):
                                                        opxx += 1 * aminoDict[iii][1]
                                                    elif float(condonBiasDict[iii][1]) < float(condonBiasDict[refcodon][1]):
                                                        nopxx += 1 * aminoDict[iii][1]
                                                    else:
                                                        print("wrong in bias !!!!!!!")
                                                else:
                                                    continue
                                            else:
                                                print("wrong in codonpos !!!!!!!!")


                                        else:
                                            nsxx += 1 * aminoDict[iii][1]

                                    rateNS = nsxx
                                    rateS = sxx
                                    rateNOP = nopxx
                                    rateOP = opxx

                                    if round(rateS + rateNS,8) != 1:
                                        print("wrong!!!!! rateS rateNS")
                                        print(rateS)
                                        print(rateNS)
                                    moniNsRateDict[sample].append(rateNS)
                                    moniSRateDict[sample].append(rateS)
                                    moniNORateDict[sample].append(rateNOP)
                                    moniORateDict[sample].append(rateOP)



                sumPns = 1E-6  # 加个极小值避免0
                sumPs = 1E-6
                sumPnop = 1E-6
                sumPop = 1E-6


                ggn = 0
                for gg in moniNsRateDict.keys():
                    ggn += 1

                    ppn = 0
                    if len(moniNsRateDict[gg]) !=0:
                        for pp in moniNsRateDict[gg]:
                            ppn += 1
                            np.random.seed(randomseed*31+19*cycleNum+cy * 13 + ggn * 111 + ppn * 41 + 37)
                            pns = np.random.binomial(1, pp, 1)
                            sumPns += float(pns / epitopeNSsiteNum)
                        # del pp
                        # del ppn
                        # del pns

                    ccn = 0
                    if len( moniSRateDict[gg]) !=0:
                        for cc in moniSRateDict[gg]:
                            ccn += 1
                            np.random.seed(randomseed*31+19*cycleNum+cy * 13 + ggn * 111 + ccn * 41 + 37)
                            ps = np.random.binomial(1, cc, 1)
                            sumPs += float(ps / epitopeSsiteNum)
                        # del cc
                        # del ccn
                        # del ps

                    aan = 0
                    if len( moniNORateDict[gg]) !=0:
                        for aa in moniNORateDict[gg]:
                            aan += 1
                            np.random.seed(randomseed*31+19*cycleNum+cy * 13 + ggn * 111 + aan * 41 + 37)
                            pnop = np.random.binomial(1, aa, 1)
                            sumPnop += float(pnop / epitopenonoptimalsiteNum)
                        # del aa
                        # del aan
                        # del pnop

                    bbn = 0
                    if len( moniORateDict[gg]) !=0:
                        for bb in moniORateDict[gg]:
                            bbn += 1
                            np.random.seed(randomseed*31+19*cycleNum+cy * 13 + ggn * 111 + bbn * 41 + 37)
                            pop = np.random.binomial(1, bb, 1)
                            sumPop += float(pop / epitopeoptimalsiteNum)
                        # del bb
                        # del bbn
                        # del pop


                ranSumPnsDict[cy] = float(sumPns)
                ranSumPsDict[cy] = float(sumPs)
                ranSumPnopDict[cy] = float(sumPnop)
                ranSumPopDict[cy] = float(sumPop)


                del sumPns
                del sumPs
                del sumPnop
                del sumPop
                del moniORateDict
                del moniNsRateDict
                del moniNORateDict
                del moniSRateDict

            # s或op为0的cycle加个极小的值 # 避免除数为0
            for ccy in ranSumPsDict.keys():
                if ranSumPsDict[ccy] == 0:
                    ranSumPsDict[ccy] = 1E-6
            for ccy in ranSumPopDict.keys():
                if ranSumPopDict[ccy] == 0:
                    ranSumPopDict[ccy] = 1E-6


            ############ 计算模拟值的95%CI,由于样本足够大，可以视为正态分布 # 废
            # 计算总体标准差
            ranSumPnsList = []
            for cc in ranSumPnsDict.keys():
                ranSumPnsList.append(ranSumPnsDict[cc])
            meanNs = np.mean(ranSumPnsList)
            stdNs = np.std(ranSumPnsList)

            ranSumPsList = []
            for cc in ranSumPsDict.keys():
                ranSumPsList.append(ranSumPsDict[cc])
            meanS = np.mean(ranSumPsList)
            stdS = np.std(ranSumPsList)

            ranSumPns_SumPsList = []
            for rrr in range(0, len(ranSumPnsList)):
                ranSumPns_SumPsList.append(ranSumPnsList[rrr] / ranSumPsList[rrr])
            meanNS_S = np.mean(ranSumPns_SumPsList)
            stdNS_S = np.std(ranSumPns_SumPsList)

            ranSumPnopList = []
            ranSumPnop_SumPsList = []
            for cc in ranSumPnopDict.keys():
                ranSumPnopList.append(ranSumPnopDict[cc])
            for rrr in range(0, len(ranSumPnopList)):
                ranSumPnop_SumPsList.append(ranSumPnopList[rrr] / ranSumPsList[rrr])
            meanNOP_S = np.mean(ranSumPnop_SumPsList)
            stdNOP_S = np.std(ranSumPnop_SumPsList)

            ranSumPopList = []
            ranSumPop_SumPsList = []
            for cc in ranSumPopDict.keys():
                ranSumPopList.append(ranSumPopDict[cc])
            for rrr in range(0, len(ranSumPopList)):
                ranSumPop_SumPsList.append(ranSumPopList[rrr] / ranSumPsList[rrr])
            meanOP_S = np.mean(ranSumPop_SumPsList)
            stdOP_S = np.std(ranSumPop_SumPsList)

            ranSumPnop_SumPopList = []
            for rrr in range(0, len(ranSumPnopList)):
                ranSumPnop_SumPopList.append(ranSumPnopList[rrr] / ranSumPopList[rrr])
            meanNOP_OP = np.mean(ranSumPnop_SumPopList)
            stdNOP_OP = np.std(ranSumPnop_SumPopList)


            # 计算偏移度
            rNSlist = []
            ranMeanNS = np.mean(ranSumPnsList)
            for uu in ranSumPnsList:
                rNS = (epitoperealSumNs - uu) / ranMeanNS
                rNSlist.append(rNS)

            rSlist = []
            ranMeanS = np.mean(ranSumPsList)
            for uu in ranSumPsList:
                rS = (epitoperealSumS - uu) / ranMeanS
                rSlist.append(rS)

            rNS_Slist = []
            ranMeanNS_S = np.mean(ranSumPns_SumPsList)
            for uu in ranSumPns_SumPsList:
                rNS_S = (epitoperealSumNs/epitoperealSumS - uu) / ranMeanNS_S
                rNS_Slist.append(rNS_S)


            rNOP_Slist = []
            ranMeanNOP_S = np.mean(ranSumPnop_SumPsList)
            for uu in ranSumPnop_SumPsList:
                rNOP_S = (epitoperealSumNop/epitoperealSumS - uu) / ranMeanNOP_S
                rNOP_Slist.append(rNOP_S)

            rOP_Slist = []
            ranMeanOP_S = np.mean(ranSumPop_SumPsList)
            for uu in ranSumPop_SumPsList:
                rOP_S = (epitoperealSumOp/epitoperealSumS - uu) / ranMeanOP_S
                rOP_Slist.append(rOP_S)

            rNOP_OPlist = []
            ranMeanNOP_OP = np.mean(ranSumPnop_SumPopList)
            for uu in ranSumPnop_SumPopList:
                rNOP_OP = (epitoperealSumNop / epitoperealSumOp - uu) / ranMeanNOP_OP
                rNOP_OPlist.append(rNOP_OP)


            # 计算偏移指数
            e = math.e

            rindexNSlist=[]
            for va in rNSlist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b = epitoperALL
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) >-1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / (epitoperALL)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexNSlist.append(ri)
            del va
            del ri
            del a
            del b
            del c


            rindexSlist=[]
            for va in rSlist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b = epitoperALL
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / (epitoperALL)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexSlist.append(ri)
            del va
            del ri
            del a
            del b
            del c


            rindexNS_Slist = []
            for va in rNS_Slist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b = epitoperALL
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / (epitoperALL)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexNS_Slist.append(ri)
            del va
            del ri
            del a
            del b
            del c


            rindexNOP_Slist = []
            for va in rNOP_Slist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b = epitoperopnop
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / (epitoperopnop)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexNOP_Slist.append(ri)
            del va
            del ri
            del a
            del b
            del c


            rindexOP_Slist = []
            for va in rOP_Slist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b =  epitoperopnop
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / ( epitoperopnop)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexOP_Slist.append(ri)
            del va
            del ri
            del a
            del b
            del c


            rindexNOP_OPlist = []
            for va in rNOP_OPlist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b =  epitoperopnop
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / ( epitoperopnop)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexNOP_OPlist.append(ri)
            del va
            del ri
            del a
            del b
            del c





            with open(outputfile9, "w") as output:
                output.write("OffsetDegree_ns_" + cladeName + "\n")
                for bb in rNSlist:
                    output.write(str(bb) + "\n")

            with open(outputfile10, "w") as output:
                output.write("OffsetDegree_s_" + cladeName + "\n")
                for bb in rSlist:
                    output.write(str(bb) + "\n")


            with open(outputfile910, "w") as output:
                output.write("OffsetDegree_ns/s_" + cladeName + "\n")
                for bb in rNS_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile11, "w") as output:
                output.write("OffsetDegree_nop/s_" + cladeName + "\n")
                for bb in rNOP_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile12, "w") as output:
                output.write("OffsetDegree_op/s_" + cladeName + "\n")
                for bb in rOP_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile1112, "w") as output:
                output.write("OffsetDegree_nop/op_" + cladeName + "\n")
                for bb in rNOP_OPlist:
                    output.write(str(bb) + "\n")


            with open(outputfile13, "w") as output:
                output.write("OffsetIndex_ns_" + cladeName + "\n")
                for bb in rindexNSlist:
                    output.write(str(bb) + "\n")

            with open(outputfile14, "w") as output:
                output.write("OffsetIndex_ns_" + cladeName + "\n")
                for bb in rindexSlist:
                    output.write(str(bb) + "\n")


            with open(outputfile15, "w") as output:
                output.write("OffsetIndex_nop/s_" + cladeName + "\n")
                for bb in rindexNOP_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile16, "w") as output:
                output.write("OffsetIndex_op/s_" + cladeName + "\n")
                for bb in rindexOP_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile17, "w") as output:
                output.write("OffsetIndex_nop/op_" + cladeName + "\n")
                for bb in rindexNOP_OPlist:
                    output.write(str(bb) + "\n")

            with open(outputfile18, "w") as output:
                output.write("OffsetIndex_ns/s_" + cladeName + "\n")
                for bb in rindexNS_Slist:
                    output.write(str(bb) + "\n")

            del ranSumPnsDict
            del ranSumPsDict
            del ranSumPnopDict
            del ranSumPopDict



            ######################### 抗原表位所在基因上的非表位区域上的模拟 ###################################################
            ranSumPnsDict = {}
            ranSumPsDict = {}
            ranSumPnopDict = {}
            ranSumPopDict = {}
            print("geneOUTEPI is starting ################")
            for cy in tqdm(range(1, cycleNum + 1),desc="表位所在基因非表位区域模拟"):
                # 随机在所有表位的所有位点上进行和真实的allSNPnum数量相同的突变
                moniNsRateDict = {}
                moniSRateDict = {}
                moniNORateDict = {}
                moniORateDict = {}

                samplen = 0
                for sample in geneOUTEPIrealALLNumDict.keys():
                    samplen += 1
                    moniNsRateDict[sample] = []
                    moniSRateDict[sample] = []
                    moniNORateDict[sample] = []
                    moniORateDict[sample] = []

                    # 在观测到真实突变的基因非表位区域上生成等数量随机位置的突变
                    rgn = 0
                    for rg in geneOUTEPIrealALLNumDict[sample].keys():
                        rgn += 1
                        fposList = []
                        sposList = []
                        tposList = []
                        if  geneOUTEPIrealALLNumDict[sample][rg][0] ==0:
                            continue
                        else:
                            # 将三个密码子位置分开,按照突变率分开抽取
                            if len(geneOUTEPIrealALLNumDict[sample][rg][1]) == 0:
                                print("WRONG!!!! no epitope in gene!!!!!!")
                            else:
                                tmpEPIposList=[]
                                for ccc in geneOUTEPIrealALLNumDict[sample][rg][1]:
                                    tmpEPIposList.append([int(ccc.split("_")[0]),int(ccc.split("_")[1])])

                                if geneDict[rg].split("_")[2] == "f":
                                    for yp in range(int(geneDict[rg].split("_")[0]), int(geneDict[rg].split("_")[1]) + 1):
                                        ypnn=0
                                        for tp in tmpEPIposList:
                                            if tp[0] <= yp <= tp[1]:
                                                ypnn +=1
                                                break
                                            else:
                                                continue
                                        if ypnn >0:
                                            continue
                                        else:
                                            if yp % 3 == 1:
                                                fposList.append(yp)
                                            elif yp % 3 == 2:
                                                sposList.append(yp)
                                            else:
                                                tposList.append(yp)
                                        del ypnn
                                        del yp

                                else:
                                    for yp in range(int(geneDict[rg].split("_")[1]), int(geneDict[rg].split("_")[0]) - 1, -1):
                                        ypnn = 0
                                        for tp in tmpEPIposList:
                                            if tp[0] <= yp <= tp[1]:
                                                ypnn += 1
                                                break
                                            else:
                                                continue
                                        if ypnn > 0:
                                            continue
                                        else:
                                            if yp % 3 == 0:
                                                fposList.append(yp)
                                            elif yp % 3 == 2:
                                                sposList.append(yp)
                                            else:
                                                tposList.append(yp)
                                        del ypnn
                                        del yp

                            if len(fposList) != len(sposList) or len(fposList) != len(tposList) or len(sposList) != len(tposList):
                                print("wrong!!!!!!! in f&s&tposList !!!!!")
                                print(str(int(geneDict[rg].split("_")[0])) + "_" + str(int(geneDict[rg].split("_")[1])))

                            for pp in range(1, geneOUTEPIrealALLNumDict[sample][rg][0] + 1):
                                np.random.seed(randomseed*31+19*cycleNum+cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 97)
                                choicePos1 = np.random.binomial(1, codonPosRateDict["1st"], 1)
                                if choicePos1 == 1:
                                    choicePOSx = 1
                                    del choicePos1
                                else:
                                    np.random.seed(randomseed*31+19*cycleNum+cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 197)
                                    secondRate = codonPosRateDict["2nd"] / (codonPosRateDict["2nd"] + codonPosRateDict["3rd"])
                                    choicePos2 = np.random.binomial(1, secondRate, 1)
                                    if choicePos2 == 1:
                                        choicePOSx = 2
                                        del choicePos2
                                    else:
                                        choicePOSx = 3
                                        del choicePos2

                                if choicePOSx == 1:
                                    np.random.seed(randomseed*31+19*cycleNum+cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 317)
                                    ranPos = np.random.choice(fposList, 1)[0]
                                elif choicePOSx == 2:
                                    np.random.seed(randomseed*31+19*cycleNum+cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 317)
                                    ranPos = np.random.choice(sposList, 1)[0]
                                elif choicePOSx == 3:
                                    np.random.seed(randomseed*31+19*cycleNum+cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 317)
                                    ranPos = np.random.choice(tposList, 1)[0]

                                refcodon = ""
                                codonpos = 0
                                if geneDict[rg].split("_")[2] == "f":
                                    if (ranPos - int(geneDict[rg].split("_")[0]) + 1) % 3 == 1:
                                        codonpos = 1
                                        refcodon = allLongNewRefseqdict[ranPos] + allLongNewRefseqdict[ranPos + 1] + \
                                                   allLongNewRefseqdict[ranPos + 2]

                                    elif (ranPos - int(geneDict[rg].split("_")[0]) + 1) % 3 == 2:
                                        codonpos = 2
                                        refcodon = allLongNewRefseqdict[ranPos - 1] + allLongNewRefseqdict[ranPos] + \
                                                   allLongNewRefseqdict[ranPos + 1]

                                    else:
                                        codonpos = 3
                                        refcodon = allLongNewRefseqdict[ranPos - 2] + allLongNewRefseqdict[ranPos - 1] + \
                                                   allLongNewRefseqdict[ranPos]

                                else:
                                    if (int(geneDict[rg].split("_")[1]) - ranPos + 1) % 3 == 1:
                                        codonpos = 1
                                        refcodon = anti_allLongNewRefseqdict[ranPos] + anti_allLongNewRefseqdict[ranPos - 1] + \
                                                   anti_allLongNewRefseqdict[ranPos - 2]

                                    elif (int(geneDict[rg].split("_")[1]) - ranPos + 1) % 3 == 2:
                                        codonpos = 2
                                        refcodon = anti_allLongNewRefseqdict[ranPos + 1] + anti_allLongNewRefseqdict[ranPos] + \
                                                   anti_allLongNewRefseqdict[ranPos - 1]

                                    else:
                                        codonpos = 3
                                        refcodon = anti_allLongNewRefseqdict[ranPos + 2] + anti_allLongNewRefseqdict[
                                            ranPos + 1] + anti_allLongNewRefseqdict[ranPos]

                                if codonpos == 0:
                                    print("wrong in random codonpos !!!!!!")
                                else:

                                    allRateDict = {}
                                    rateα = titvDict["ti"]
                                    rateβ1 = titvDict["tvtype1"]
                                    rateβ2 = titvDict["tvtype2"]
                                    rateΘ1 = tamuraRateDict["rate1"]
                                    rateΘ2 = tamuraRateDict["rate2"]
                                    allRateDict["A"] = 2 * rateβ1 * rateΘ2 + 2 * rateβ2 * rateΘ1 + rateα * rateΘ1
                                    allRateDict["T"] = 2 * rateβ1 * rateΘ2 + 2 * rateβ2 * rateΘ1 + rateα * rateΘ1
                                    allRateDict["C"] = 2 * rateβ1 * rateΘ1 + 2 * rateβ2 * rateΘ2 + rateα * rateΘ2
                                    allRateDict["G"] = 2 * rateβ1 * rateΘ1 + 2 * rateβ2 * rateΘ2 + rateα * rateΘ2
                                    snpRateDict = {}
                                    snpRateDict["AT"] = 2 * rateβ1 * rateΘ2
                                    snpRateDict["AC"] = 2 * rateβ2 * rateΘ1
                                    snpRateDict["AG"] = rateα * rateΘ1
                                    snpRateDict["TA"] = 2 * rateβ1 * rateΘ2
                                    snpRateDict["TC"] = rateα * rateΘ1
                                    snpRateDict["TG"] = 2 * rateβ2 * rateΘ1
                                    snpRateDict["CA"] = 2 * rateβ2 * rateΘ2
                                    snpRateDict["CT"] = rateα * rateΘ2
                                    snpRateDict["CG"] = 2 * rateβ1 * rateΘ1
                                    snpRateDict["GA"] = rateα * rateΘ2
                                    snpRateDict["GT"] = 2 * rateβ2 * rateΘ2
                                    snpRateDict["GC"] = 2 * rateβ1 * rateΘ1

                                    aminoDict = {}
                                    if codonpos == 1:
                                        allrate = allRateDict[refcodon[0]]
                                        for k in condonBiasDict.keys():
                                            if refcodon[1:] == k[1:] and refcodon[0] != k[0]:
                                                snpRate = snpRateDict[refcodon[0] + k[0]]
                                                aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                                del snpRate
                                            else:
                                                continue
                                        del allrate
                                    elif codonpos == 2:
                                        allrate = allRateDict[refcodon[1]]
                                        for k in condonBiasDict.keys():
                                            if refcodon[0] == k[0] and refcodon[-1] == k[-1] and refcodon[1] != k[1]:
                                                snpRate = snpRateDict[refcodon[1] + k[1]]
                                                aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                                del snpRate
                                            else:
                                                continue
                                        del allrate
                                    elif codonpos == 3:
                                        allrate = allRateDict[refcodon[2]]
                                        for k in condonBiasDict.keys():
                                            if refcodon[:2] == k[:2] and refcodon[2] != k[2]:
                                                snpRate = snpRateDict[refcodon[2] + k[2]]
                                                aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
                                                del snpRate
                                            else:
                                                continue
                                        del allrate
                                    else:
                                        print("wrong in codonpos !!!!!")

                                    # print(aminoDict)
                                    sxx = 0
                                    nsxx = 0
                                    nopxx = 0
                                    opxx = 0
                                    for iii in aminoDict.keys():
                                        if aminoDict[iii][0] == condonBiasDict[refcodon][0]:
                                            # print(aminoDict[iii][0])
                                            sxx += 1 * aminoDict[iii][1]
                                            if codonpos == 1:
                                                if baseDict[iii[0]] != refcodon[0]:
                                                    if float(condonBiasDict[iii][1]) > float(condonBiasDict[refcodon][1]):
                                                        opxx += 1 * aminoDict[iii][1]
                                                    elif float(condonBiasDict[iii][1]) < float(condonBiasDict[refcodon][1]):
                                                        nopxx += 1 * aminoDict[iii][1]
                                                    else:
                                                        print("wrong in bias !!!!!!!")
                                                else:
                                                    continue
                                            elif codonpos == 2:
                                                if baseDict[iii[1]] != refcodon[1]:
                                                    if float(condonBiasDict[iii][1]) > float(condonBiasDict[refcodon][1]):
                                                        opxx += 1 * aminoDict[iii][1]
                                                    elif float(condonBiasDict[iii][1]) < float(condonBiasDict[refcodon][1]):
                                                        nopxx += 1 * 1 * aminoDict[iii][1]
                                                    else:
                                                        print("wrong in bias !!!!!!!")
                                                else:
                                                    continue
                                            elif codonpos == 3:
                                                if baseDict[iii[2]] != refcodon[2]:
                                                    if float(condonBiasDict[iii][1]) > float(condonBiasDict[refcodon][1]):
                                                        opxx += 1 * aminoDict[iii][1]
                                                    elif float(condonBiasDict[iii][1]) < float(condonBiasDict[refcodon][1]):
                                                        nopxx += 1 * aminoDict[iii][1]
                                                    else:
                                                        print("wrong in bias !!!!!!!")
                                                else:
                                                    continue
                                            else:
                                                print("wrong in codonpos !!!!!!!!")


                                        else:
                                            nsxx += 1 * aminoDict[iii][1]

                                    rateNS = nsxx
                                    rateS = sxx
                                    rateNOP = nopxx
                                    rateOP = opxx

                                    if round(rateS + rateNS, 8) != 1:
                                        print("wrong!!!!! rateS rateNS")
                                        print(rateS)
                                        print(rateNS)
                                    moniNsRateDict[sample].append(rateNS)
                                    moniSRateDict[sample].append(rateS)
                                    moniNORateDict[sample].append(rateNOP)
                                    moniORateDict[sample].append(rateOP)

                sumPns = 1E-6  # 加个极小值避免0
                sumPs = 1E-6
                sumPnop = 1E-6
                sumPop = 1E-6

                ggn = 0
                for gg in moniNsRateDict.keys():
                    ggn += 1

                    ppn = 0
                    for pp in moniNsRateDict[gg]:
                        ppn += 1
                        np.random.seed(randomseed*31+19*cycleNum+cy * 13 + ggn * 111 + ppn * 41 + 37)
                        pns = np.random.binomial(1, pp, 1)
                        sumPns += float(pns / (geneNSsiteNum - epitopeNSsiteNum))


                    ccn = 0
                    for cc in moniSRateDict[gg]:
                        ccn += 1
                        np.random.seed(randomseed*31+19*cycleNum+cy * 13 + ggn * 111 + ccn * 41 + 37)
                        ps = np.random.binomial(1, cc, 1)
                        sumPs += float(ps / (geneSsiteNum - epitopeSsiteNum))


                    aan = 0
                    for aa in moniNORateDict[gg]:
                        aan += 1
                        np.random.seed(randomseed*31+19*cycleNum+cy * 13 + ggn * 111 + aan * 41 + 37)
                        pnop = np.random.binomial(1, aa, 1)
                        sumPnop += float(pnop / (genenonoptimalsiteNum- epitopenonoptimalsiteNum))


                    bbn = 0
                    for bb in moniORateDict[gg]:
                        bbn += 1
                        np.random.seed(randomseed*31+19*cycleNum+cy * 13 + ggn * 111 + bbn * 41 + 37)
                        pop = np.random.binomial(1, bb, 1)
                        sumPop += float(pop / (geneoptimalsiteNum- epitopeoptimalsiteNum))




                ranSumPnsDict[cy] = float(sumPns)
                ranSumPsDict[cy] = float(sumPs)
                ranSumPnopDict[cy] = float(sumPnop)
                ranSumPopDict[cy] = float(sumPop)

                del sumPns
                del sumPs
                del sumPnop
                del sumPop
                del moniNsRateDict
                del moniSRateDict
                del moniNORateDict
                del moniORateDict



            # s或op为0的cycle加个极小的值 # 避免除数为0
            for ccy in ranSumPsDict.keys():
                if ranSumPsDict[ccy] == 0:
                    ranSumPsDict[ccy] = 1E-6
            for ccy in ranSumPopDict.keys():
                if ranSumPopDict[ccy] == 0:
                    ranSumPopDict[ccy] = 1E-6

            ############ 计算模拟值的95%CI,由于样本足够大，可以视为正态分布 # 废
            # 计算总体标准差
            ranSumPnsList = []
            for cc in ranSumPnsDict.keys():
                ranSumPnsList.append(ranSumPnsDict[cc])
            meanNs = np.mean(ranSumPnsList)
            stdNs = np.std(ranSumPnsList)

            ranSumPsList = []
            for cc in ranSumPsDict.keys():
                ranSumPsList.append(ranSumPsDict[cc])
            meanS = np.mean(ranSumPsList)
            stdS = np.std(ranSumPsList)

            ranSumPns_SumPsList = []
            for rrr in range(0, len(ranSumPnsList)):
                ranSumPns_SumPsList.append(ranSumPnsList[rrr] / ranSumPsList[rrr])
            meanNS_S = np.mean(ranSumPns_SumPsList)
            stdNS_S = np.std(ranSumPns_SumPsList)

            ranSumPnopList = []
            ranSumPnop_SumPsList = []
            for cc in ranSumPnopDict.keys():
                ranSumPnopList.append(ranSumPnopDict[cc])
            for rrr in range(0, len(ranSumPnopList)):
                ranSumPnop_SumPsList.append(ranSumPnopList[rrr] / ranSumPsList[rrr])
            meanNOP_S = np.mean(ranSumPnop_SumPsList)
            stdNOP_S = np.std(ranSumPnop_SumPsList)

            ranSumPopList = []
            ranSumPop_SumPsList = []
            for cc in ranSumPopDict.keys():
                ranSumPopList.append(ranSumPopDict[cc])
            for rrr in range(0, len(ranSumPopList)):
                ranSumPop_SumPsList.append(ranSumPopList[rrr] / ranSumPsList[rrr])
            meanOP_S = np.mean(ranSumPop_SumPsList)
            stdOP_S = np.std(ranSumPop_SumPsList)

            ranSumPnop_SumPopList = []
            for rrr in range(0, len(ranSumPnopList)):
                ranSumPnop_SumPopList.append(ranSumPnopList[rrr] / ranSumPopList[rrr])
            meanNOP_OP = np.mean(ranSumPnop_SumPopList)
            stdNOP_OP = np.std(ranSumPnop_SumPopList)


            # 计算偏移度
            rNSlist = []
            ranMeanNS = np.mean(ranSumPnsList)
            for uu in ranSumPnsList:
                rNS = (geneOUTEPIrealSumNs - uu) / ranMeanNS
                rNSlist.append(rNS)

            rSlist = []
            ranMeanS = np.mean(ranSumPsList)
            for uu in ranSumPsList:
                rS = (geneOUTEPIrealSumS - uu) / ranMeanS
                rSlist.append(rS)

            rNS_Slist = []
            ranMeanNS_S = np.mean(ranSumPns_SumPsList)
            for uu in ranSumPns_SumPsList:
                rNS_S = ((geneOUTEPIrealSumNs / geneOUTEPIrealSumS) - uu) / ranMeanNS_S
                rNS_Slist.append(rNS_S)

            rNOP_Slist = []
            ranMeanNOP_S = np.mean(ranSumPnop_SumPsList)
            for uu in ranSumPnop_SumPsList:
                rNOP_S = ((geneOUTEPIrealSumNop / geneOUTEPIrealSumS) - uu) / ranMeanNOP_S
                rNOP_Slist.append(rNOP_S)

            rOP_Slist = []
            ranMeanOP_S = np.mean(ranSumPop_SumPsList)
            for uu in ranSumPop_SumPsList:
                rOP_S = ((geneOUTEPIrealSumOp / geneOUTEPIrealSumS) - uu) / ranMeanOP_S
                rOP_Slist.append(rOP_S)

            rNOP_OPlist = []
            ranMeanNOP_OP = np.mean(ranSumPnop_SumPopList)
            for uu in ranSumPnop_SumPopList:
                rNOP_OP = ((geneOUTEPIrealSumNop / geneOUTEPIrealSumOp) - uu) / ranMeanNOP_OP
                rNOP_OPlist.append(rNOP_OP)

            # 计算偏移指数
            e = math.e

            rindexNSlist = []
            for va in rNSlist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b = geneOUTEPIrALL
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / (geneOUTEPIrALL)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexNSlist.append(ri)
            del va
            del ri
            del a
            del b
            del c

            rindexSlist = []
            for va in rSlist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b = geneOUTEPIrALL
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / (geneOUTEPIrALL)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexSlist.append(ri)
            del va
            del ri
            del a
            del b
            del c

            rindexNS_Slist = []
            for va in rNS_Slist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b = geneOUTEPIrALL
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / (geneOUTEPIrALL)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexNS_Slist.append(ri)
            del va
            del ri
            del a
            del b
            del c

            rindexNOP_Slist = []
            for va in rNOP_Slist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b = geneOUTEPIropnop
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / (geneOUTEPIropnop)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexNOP_Slist.append(ri)
            del va
            del ri
            del a
            del b
            del c

            rindexOP_Slist = []
            for va in rOP_Slist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b = geneOUTEPIropnop
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / (geneOUTEPIropnop)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexOP_Slist.append(ri)
            del va
            del ri
            del a
            del b
            del c

            rindexNOP_OPlist = []
            for va in rNOP_OPlist:
                if va == 0:
                    ri = 0
                elif va > 0:
                    a = float(np.tanh(va)) + 1
                    b = geneOUTEPIropnop
                    c = - math.log(b + 1, e) * math.log(va + 1, e)
                    ri = 1 - (a ** c)
                elif va < 0:
                    if float(np.tanh(va)) > -1:
                        a = 1 / (float(np.tanh(va)) + 1)
                    else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                        a = 1E10
                    b = 1 / (geneOUTEPIropnop)
                    c = - math.log(b + 1, e) * math.log(-va + 1, e)
                    ri = (a ** c) - 1
                rindexNOP_OPlist.append(ri)
            del va
            del ri
            del a
            del b
            del c




            with open(outputfile29, "w") as output:
                output.write("OffsetDegree_ns_" + cladeName + "\n")
                for bb in rNSlist:
                    output.write(str(bb) + "\n")

            with open(outputfile210, "w") as output:
                output.write("OffsetDegree_s_" + cladeName + "\n")
                for bb in rSlist:
                    output.write(str(bb) + "\n")

            with open(outputfile2910, "w") as output:
                output.write("OffsetDegree_ns/s_" + cladeName + "\n")
                for bb in rNS_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile211, "w") as output:
                output.write("OffsetDegree_nop/s_" + cladeName + "\n")
                for bb in rNOP_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile212, "w") as output:
                output.write("OffsetDegree_op/s_" + cladeName + "\n")
                for bb in rOP_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile21112, "w") as output:
                output.write("OffsetDegree_nop/op_" + cladeName + "\n")
                for bb in rNOP_OPlist:
                    output.write(str(bb) + "\n")

            with open(outputfile213, "w") as output:
                output.write("OffsetIndex_ns_" + cladeName + "\n")
                for bb in rindexNSlist:
                    output.write(str(bb) + "\n")

            with open(outputfile214, "w") as output:
                output.write("OffsetIndex_ns_" + cladeName + "\n")
                for bb in rindexSlist:
                    output.write(str(bb) + "\n")

            with open(outputfile215, "w") as output:
                output.write("OffsetIndex_nop/s_" + cladeName + "\n")
                for bb in rindexNOP_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile216, "w") as output:
                output.write("OffsetIndex_op/s_" + cladeName + "\n")
                for bb in rindexOP_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile217, "w") as output:
                output.write("OffsetIndex_nop/op_" + cladeName + "\n")
                for bb in rindexNOP_OPlist:
                    output.write(str(bb) + "\n")

            with open(outputfile218, "w") as output:
                output.write("OffsetIndex_ns/s_" + cladeName + "\n")
                for bb in rindexNS_Slist:
                    output.write(str(bb) + "\n")



print("finished!!!")













