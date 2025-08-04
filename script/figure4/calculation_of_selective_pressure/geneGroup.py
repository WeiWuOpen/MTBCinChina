


cladeNameList=["4.2cnONLYchina"]
geneGroupList = ["MIKGs","essentialGenes"]


# 模拟的次数
cycleNum=1000
epsilon = 1e-10  # 设定一个非常小的数作为误差范围,用来代替浮点数里的0，因为浮点数位数的问题，直接==0或==0.0无法判定


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


for cladeName in cladeNameList:
    inputlistFile = r"...\script_and_testData\testData\sampleList\\" + cladeName + ".txt"
    inputvcfFolder = r"...\script_and_testData\testData\vcf_4.2_CRMA\\"
    inputrefseqSNPFile_new = r"...\script_and_testData\testData\\"+cladeName.replace("ONLYchina","").replace("NOchina","")+"_CRMA.fasta"

    if cladeName[0] == "2":
        inputrefseqPOSFile_new = r"...\script_and_testData\testData\allL2_withGeolocInfo_refaMTBC_withprd.pos.txt"
    elif cladeName[0] == "4":
        inputrefseqPOSFile_new = r"...\script_and_testData\testData\allL4_withGeolocInfo_refaMTBC_withprd.pos.txt"
    else:
        print("wrong in inputrefseqPOSFile_new!!!!!!")


    tamuraRateDict = {}
    with open(inputTamuraRateFile, "r") as input:
        for ll in input:
            llx = ll.strip().split()
            if llx[0] != "Rate1":
                tamuraRateDict["rate1"] = float(llx[0])
                tamuraRateDict["rate2"] = float(llx[1])

    codonPosRateDict = {}
    with open(inputcodonPOSsnpRateFile, "r") as input:
        for l in input:
            lx = l.strip().split()
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

    strainlist = []
    with open(inputlistFile, "r") as input:
        for l in input:
            strainlist.append(l.strip())
    # 避免重复样本
    strainlist = np.unique(strainlist)

    sampleNum = len(strainlist)

    allgeneDict = {}
    with open(inputAllGeneFile, "r") as input:
        for l in input:
            lx = l.strip().split()
            if lx[2] == "CDS":
                if lx[5] == "+":
                    allgeneDict[lx[0]] = lx[3] + "_" + lx[4] + "_" + "f"
                else:
                    allgeneDict[lx[0]] = lx[3] + "_" + lx[4] + "_" + "c"

    vcfdict = {}
    beddict = {}
    allvcfSNPnum = 0
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
                    allvcfSNPnum += 1

    allvcfSNPrate = allvcfSNPnum / len(refseq)




    ########## 全基因组上全部CDS的全部待选ns，s位点数目与待选optimal,nonoptimal位点数目 ####################################
    # print("############ start whole genome genes candidate ns s et al. ##############")
    # allNSsiteNum = 0
    # allSsiteNum = 0
    # alloptimalsiteNum = 0
    # allnonoptimalsiteNum = 0
    # allnormalsiteNum = 0
    # allCDSLength = 0
    # # ggg=0
    # for gg in tqdm(allgeneDict.keys()):
    #     # ggg+=1
    #     # print(str(ggg)+"/"+str(len(allgeneDict.keys())))
    #     for pp in range(int(allgeneDict[gg].split("_")[0]), int(allgeneDict[gg].split("_")[1]) + 1, 3):
    #         allCDSLength += 3
    #         if allgeneDict[gg].split("_")[2] == "f":
    #             codon = allLongNewRefseqdict[pp] + allLongNewRefseqdict[pp + 1] + allLongNewRefseqdict[pp + 2]
    #         elif allgeneDict[gg].split("_")[2] == "c":
    #             codon = anti_allLongNewRefseqdict[pp + 2] + anti_allLongNewRefseqdict[pp + 1] + \
    #                     anti_allLongNewRefseqdict[pp]
    #         else:
    #             print("WRONG IN allgeneDict!!!!!!!")
    #
    #         for yi in range(1, 4):
    #             aminoDict = {}
    #             if yi == 1:
    #                 allrate = allRateDict[codon[0]]
    #                 for k in condonBiasDict.keys():
    #                     if codon[1:] == k[1:] and codon[0] != k[0]:
    #                         snpRate = snpRateDict[codon[0] + k[0]]
    #                         aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
    #                         del snpRate
    #                     else:
    #                         continue
    #                 del allrate
    #             elif yi == 2:
    #                 allrate = allRateDict[codon[1]]
    #                 for k in condonBiasDict.keys():
    #                     if codon[0] == k[0] and codon[-1] == k[-1] and codon[1] != k[1]:
    #                         snpRate = snpRateDict[codon[1] + k[1]]
    #                         aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
    #                         del snpRate
    #                     else:
    #                         continue
    #                 del allrate
    #             elif yi == 3:
    #                 allrate = allRateDict[codon[2]]
    #                 for k in condonBiasDict.keys():
    #                     if codon[:2] == k[:2] and codon[2] != k[2]:
    #                         snpRate = snpRateDict[codon[2] + k[2]]
    #                         aminoDict[k] = [condonBiasDict[k][0], snpRate / allrate]
    #                         del snpRate
    #                     else:
    #                         continue
    #                 del allrate
    #             else:
    #                 print("wrong in codonPos !!!!!")
    #
    #             # print(aminoDict)
    #             sxx = 0
    #             nsxx = 0
    #             nopxx = 0
    #             opxx = 0
    #             normalxx = 0
    #             for iii in aminoDict.keys():
    #                 if aminoDict[iii][0] == condonBiasDict[codon][0]:
    #                     # print(aminoDict[iii][0])
    #                     sxx += 1 * aminoDict[iii][1]
    #                     if yi == 1:
    #                         if baseDict[iii[0]] != codon[0]:  # 由于 AT， CG互换的比例接近，所以不视为发生了两种密码子的更换
    #                             if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
    #                                 opxx += 1 * aminoDict[iii][1]
    #                             elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
    #                                 nopxx += 1 * aminoDict[iii][1]
    #                             else:
    #                                 print("wrong in bias !!!!!!!")
    #                         else:
    #                             normalxx += 1 * aminoDict[iii][1]
    #                     elif yi == 2:
    #                         if baseDict[iii[1]] != codon[1]:
    #                             if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
    #                                 opxx += 1 * aminoDict[iii][1]
    #                             elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
    #                                 nopxx += 1 * aminoDict[iii][1]
    #                             else:
    #                                 print("wrong in bias !!!!!!!")
    #                         else:
    #                             normalxx += 1 * aminoDict[iii][1]
    #                     elif yi == 3:
    #                         if baseDict[iii[2]] != codon[2]:
    #                             if float(condonBiasDict[iii][1]) > float(condonBiasDict[codon][1]):
    #                                 opxx += 1 * aminoDict[iii][1]
    #                             elif float(condonBiasDict[iii][1]) < float(condonBiasDict[codon][1]):
    #                                 nopxx += 1 * aminoDict[iii][1]
    #                             else:
    #                                 print("wrong in bias !!!!!!!")
    #                         else:
    #                             normalxx += 1 * aminoDict[iii][1]
    #                     else:
    #                         print("wrong in codonpos !!!!!!!!")
    #
    #
    #                 else:
    #                     nsxx += 1 * aminoDict[iii][1]
    #
    #             alloptimalsiteNum += opxx
    #             allnonoptimalsiteNum += nopxx
    #             allSsiteNum += sxx
    #             allNSsiteNum += nsxx
    #             allnormalsiteNum += normalxx
    #
    #             del yi
    #             del opxx
    #             del nopxx
    #             del sxx
    #             del nsxx
    #             del normalxx
    #             del aminoDict
    #
    #     del codon

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
    # 只要全基因组的CDS没有变化，上面的值是固定的，没必要每次都算
    allCDSLength = 3958911
    allNSsiteNum = 2843652.0343105034
    allSsiteNum = 1115258.9657158025
    alloptimalsiteNum = 178470.48630486924
    allnonoptimalsiteNum = 776737.504495235
    allnormalsiteNum = 160050.97490044663
    ######################################################################################





    ######################################################################################

    ####### 全基因组上全部CDS的各种SNP ####################################################
    print("Start 全基因组上全部CDS的各种SNP")
    allrealALLNumDict = {}
    allrealNsNumDict = {}
    allrealSNumDict = {}
    allrealOptimalNumDict = {}
    allrealNonoptimalNumDict = {}
    for ss in tqdm(strainlist):
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
                    if int(allgeneDict[g].split("_")[0]) <= int(p.split("_")[0]) <= int(allgeneDict[g].split("_")[1]):
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

        del codon
        del refcodon

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

    print("########################")
    print(cladeName)
    print(allrealSumNSnum)
    print(allrealSumSnum)
    print(allrealSumNopnum)
    print(allrealSumOpnum)
    print("#########################")
    ####################################################################################################

    for geneGroupName in geneGroupList:
        try:
            print(cladeName+"_"+geneGroupName + " is working #######")

            inputGenefile = r"...\script_and_testData\testData\geneGroup\\" + geneGroupName + ".txt"

            outputFolder2 = r"...\script_and_testData\testData\output\calculation_of_selective_pressure\geneGroup\\" + geneGroupName + "\\"
            outputFolder= outputFolder2 + "refcladeMRCA_"+str(cycleNum)+"\\"


            outputfile910 = outputFolder+ cladeName + "_offsetDgree_ns_dividedbyS.txt"
            outputfile1112 = outputFolder+ cladeName + "_offsetDgree_nop_dividedbyOp.txt"
            outputfile17 = outputFolder + cladeName + "_offsetScore_nop_dividedbyop.txt"
            outputfile18 = outputFolder + cladeName + "_offsetScore_ns_dividedbyS.txt"
            #########################################################################################################

            if not os.path.exists(outputFolder2):
                os.mkdir(outputFolder2)
            if not os.path.exists(outputFolder):
                os.mkdir(outputFolder)



            geneDict = {}
            with open(inputGenefile, "r") as input:
                for l in input:
                    if l.strip() !="":
                        lx = l.strip().split("\t")
                        if lx[4] == "+": # 正向基因
                            geneDict[lx[0]] = lx[2] + "_" + lx[3] + "_" + "f"
                        else:
                            geneDict[lx[0]] = lx[2] + "_" + lx[3] + "_" + "c"



            ########## 目标CDS上的全部待选ns，s位点数目与待选optimal,nonoptimal位点数目 ####################################
            NSsiteNum=0
            SsiteNum=0
            optimalsiteNum=0
            nonoptimalsiteNum=0
            normalsiteNum=0
            allGeneLength=0
            for gg in geneDict.keys():
                for pp in range(int(geneDict[gg].split("_")[0]), int(geneDict[gg].split("_")[1]) + 1, 3):
                    allGeneLength+=3
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
                        normalxx =0
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
                                        normalxx +=  1 * aminoDict[iii][1]
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


                        optimalsiteNum += opxx
                        nonoptimalsiteNum += nopxx
                        SsiteNum += sxx
                        NSsiteNum += nsxx
                        normalsiteNum += normalxx

                        del yi
                        del opxx
                        del nopxx
                        del sxx
                        del nsxx
                        del normalxx
                        del aminoDict

                del codon

            print("target gene group: ")
            print(allGeneLength)
            print(NSsiteNum)
            print(SsiteNum)
            print(optimalsiteNum)
            print(nonoptimalsiteNum)
            print(optimalsiteNum+normalsiteNum+nonoptimalsiteNum)
            print("#######################")






            ################################################################################################################
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
                        for g in geneDict.keys():
                            if int(geneDict[g].split("_")[0]) <= int(p.split("_")[0]) <= int(geneDict[g].split("_")[1]):
                                if g not in realALLNumDict[ss].keys():
                                    realALLNumDict[ss][g] = 1
                                else:
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



            realSumNs = 0
            realSumS = 0
            realSumNop = 0
            realSumOp = 0
            realSumNSnum = 0
            realSumSnum = 0
            realSumNopnum=0
            realSumOpnum=0
            for real in realNsNumDict.keys():
                realSumNs += float(realNsNumDict[real]) / NSsiteNum
                realSumNSnum += float(realNsNumDict[real])
            for real in realSNumDict.keys():
                realSumS += float(realSNumDict[real]) / SsiteNum
                realSumSnum += float(realSNumDict[real])
            for real in realNonoptimalNumDict.keys():
                realSumNop += float(realNonoptimalNumDict[real]) / nonoptimalsiteNum
                realSumNopnum += float(realNonoptimalNumDict[real])
            for real in realOptimalNumDict.keys():
                realSumOp += float(realOptimalNumDict[real]) / optimalsiteNum
                realSumOpnum += float(realOptimalNumDict[real])

            realSumNSS = (realSumNSnum + realSumSnum) / ((NSsiteNum + SsiteNum) * len(strainlist))
            realSumnopop = (realSumOpnum + realSumNopnum) / ((nonoptimalsiteNum + optimalsiteNum) * len(strainlist))

            allnssSNPrate = (allrealSumNSnum + allrealSumSnum) / ((allNSsiteNum + allSsiteNum) * len(strainlist))
            allnopopSNPrate = (allrealSumOpnum + allrealSumNopnum) / ((allnonoptimalsiteNum + alloptimalsiteNum) * len(strainlist))

            rALL = realSumNSS / allnssSNPrate
            ropnop = realSumnopop / allnopopSNPrate
            print(rALL)
            print(ropnop)



            #####################  模拟  #####################################################################
            print("####### binomial test is working")
            cycleNum = cycleNum  # 模拟的次数
            ranSumPnsDict = {}
            ranSumPsDict = {}
            ranSumPnopDict = {}
            ranSumPopDict = {}

            ranSumPnsnormalDict = {}


            for cy in tqdm(range(1, cycleNum + 1),desc="binomial test cycle"):
                # for cccc in range(0, cycleNum + 1, int(cycleNum/10)):
                #     if cy == cccc:
                #         print(str(cy) + "/" + str(cycleNum))

                # 随机在所有基因的所有位点上进行和真实的allSNPnum数量相同的突变
                moniNsRateDict = {}
                moniSRateDict = {}
                moniNORateDict = {}
                moniORateDict = {}

                moniNsNormalRateDict = {}


                samplen = 0
                for sample in realALLNumDict.keys():
                    samplen += 1
                    moniNsRateDict[sample] = []
                    moniSRateDict[sample] = []
                    moniNORateDict[sample] = []
                    moniORateDict[sample] = []

                    moniNsNormalRateDict[sample] = []

                    # 在观测到真实突变的基因上生成等数量随机位置的突变
                    rgn = 0
                    for rg in realALLNumDict[sample].keys():
                        rgn += 1
                        fposList = []
                        sposList = []
                        tposList = []
                        # 将三个密码子位置分开,按照突变率分开抽取
                        if geneDict[rg].split("_")[2] == "f":
                            for yp in range(int(geneDict[rg].split("_")[0]), int(geneDict[rg].split("_")[1]) + 1):
                                if yp % 3 == 1:
                                    fposList.append(yp)
                                elif yp % 3 == 2:
                                    sposList.append(yp)
                                else:
                                    tposList.append(yp)
                        else:
                            for yp in range(int(geneDict[rg].split("_")[1]), int(geneDict[rg].split("_")[0]) - 1,-1):
                                if yp % 3 == 0:
                                    fposList.append(yp)
                                elif yp % 3 == 2:
                                    sposList.append(yp)
                                else:
                                    tposList.append(yp)

                        if len(fposList) != len(sposList) or len(fposList) != len(tposList) or len(sposList) != len(tposList):
                            print("wrong!!!!!!! in f&s&tposList !!!!!")
                            print(str(int(geneDict[rg].split("_")[0]))+"_"+str(int(geneDict[rg].split("_")[1])))

                        for pp in range(1, realALLNumDict[sample][rg] + 1):
                            np.random.seed(cycleNum * 19 +cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 97)
                            choicePos1=np.random.binomial(1,codonPosRateDict["1st"],1)
                            if choicePos1 ==1:
                                choicePOSx =1
                                del choicePos1
                            else:
                                np.random.seed(cycleNum * 19 + cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 197)
                                secondRate=codonPosRateDict["2nd"] / (codonPosRateDict["2nd"]+codonPosRateDict["3rd"])
                                choicePos2=np.random.binomial(1,secondRate,1)
                                if choicePos2 ==1:
                                    choicePOSx =2
                                    del choicePos2
                                else:
                                    choicePOSx =3
                                    del choicePos2

                            if choicePOSx ==1:
                                np.random.seed(cycleNum * 19 + cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 317)
                                ranPos = np.random.choice(fposList,1)[0]
                            elif choicePOSx ==2:
                                np.random.seed(cycleNum * 19 + cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 317)
                                ranPos = np.random.choice(sposList, 1)[0]
                            elif choicePOSx ==3:
                                np.random.seed(cycleNum * 19 + cy * 13 + samplen * 37 + rgn * 21 + pp * 17 + 317)
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


                sumPns = 0
                sumPs = 0
                sumPnop = 0
                sumPop = 0

                ggn = 0
                for gg in moniNsRateDict.keys():
                    ggn += 1

                    ppn = 0
                    for pp in moniNsRateDict[gg]:
                        ppn += 1
                        np.random.seed(cycleNum * 19 + cy * 13 + ggn * 111 + ppn * 41 + 37)
                        pns = np.random.binomial(1, pp, 1)
                        sumPns += pns / NSsiteNum
                    # del pp
                    # del ppn
                    # del pns
                    ccn = 0
                    for cc in moniSRateDict[gg]:
                        ccn += 1
                        np.random.seed(cycleNum * 19 + cy * 13 + ggn * 111 + ccn * 41 + 37)
                        ps = np.random.binomial(1, cc, 1)
                        sumPs += ps / SsiteNum
                    # del cc
                    # del ccn
                    # del ps
                    aan = 0
                    for aa in moniNORateDict[gg]:
                        aan += 1
                        np.random.seed(cycleNum * 19 + cy * 13 + ggn * 111 + aan * 41 + 37)
                        pnop = np.random.binomial(1, aa, 1)
                        sumPnop += pnop / nonoptimalsiteNum
                    # del aa
                    # del aan
                    # del pnop
                    bbn = 0
                    for bb in moniORateDict[gg]:
                        bbn += 1
                        np.random.seed(cycleNum * 19 + cy * 13 + ggn * 111 + bbn * 41 + 37)
                        pop = np.random.binomial(1, bb, 1)
                        sumPop += pop / optimalsiteNum
                    # del bb
                    # del bbn
                    # del pop


                ranSumPnsDict[cy] = float(sumPns)
                ranSumPsDict[cy] = float(sumPs)
                ranSumPnopDict[cy] = float(sumPnop)
                ranSumPopDict[cy] = float(sumPop)


            ############ 计算模拟值的95%CI,由于样本足够大，可以视为正态分布 # 废
            # 计算总体标准差
            ranSumPnsList = []
            for cc in ranSumPnsDict.keys():
                ranSumPnsList.append(ranSumPnsDict[cc])
            ranSumPsList = []
            for cc in ranSumPsDict.keys():
                ranSumPsList.append(ranSumPsDict[cc])
            ranSumPnopList = []
            for cc in ranSumPnopDict.keys():
                ranSumPnopList.append(ranSumPnopDict[cc])
            ranSumPopList = []
            for cc in ranSumPopDict.keys():
                ranSumPopList.append(ranSumPopDict[cc])


            ranSumPns_SumPsList=[]
            for rrr in range(0,len(ranSumPnsList)):
                if ranSumPnsList[rrr] < epsilon:
                    ranSumPns_SumPsList.append(0.0)
                else:
                    if ranSumPsList[rrr] < epsilon:
                        ranSumPsList[rrr] = 1e-4  # 加个极小值避免除数为0
                    ranSumPns_SumPsList.append(ranSumPnsList[rrr] / ranSumPsList[rrr])
            meanNS_S = np.mean(ranSumPns_SumPsList)
            stdNS_S = np.std(ranSumPns_SumPsList)





            ranSumPnop_SumPopList = []
            for rrr in range(0, len(ranSumPnopList)):
                if ranSumPnopList[rrr] < epsilon:
                    ranSumPnop_SumPopList.append(0)
                else:
                    if ranSumPopList[rrr] < epsilon:
                        ranSumPopList[rrr] = 1e-4  # 加个极小值避免除数为0
                    ranSumPnop_SumPopList.append(ranSumPnopList[rrr] / ranSumPopList[rrr])
            meanNOP_OP = np.mean(ranSumPnop_SumPopList)
            # nopopCi = stats.norm.interval(0.95, loc=meanNOP_OP, scale=stats.sem(ranSumPnop_SumPopList))
            stdNOP_OP = np.std(ranSumPnop_SumPopList)



            # 计算偏移度
            rNS_Slist = []
            ranMeanNS_S = float(np.mean(ranSumPns_SumPsList))
            for uu in ranSumPns_SumPsList:
                if realSumS < epsilon and realSumNs < epsilon:  # 没有突变发生
                    rNS_S = 0.0
                    rNS_Slist.append(rNS_S)
                else:
                    if realSumS < epsilon:
                        realSumS = 1e-4
                    if ranMeanNS_S < epsilon:
                        ranMeanNS_S = 0.1
                    rNS_S = (realSumNs / realSumS - uu) / ranMeanNS_S
                    rNS_Slist.append(rNS_S)

            rNOP_OPlist = []
            ranMeanNOP_OP = float(np.mean(ranSumPnop_SumPopList))
            for uu in ranSumPnop_SumPopList:
                if realSumS < epsilon:  # 没有同义突变发生
                    rNOP_OP = 0.0
                    rNOP_OPlist.append(rNOP_OP)
                else:
                    if realSumOp < epsilon:  # 没有向optimal密码子突变发生
                        realSumOp = 1e-4
                    if ranMeanNOP_OP < epsilon:
                        ranMeanNOP_OP = 0.1

                    rNOP_OP = (realSumNop / realSumOp - uu) / ranMeanNOP_OP
                    rNOP_OPlist.append(rNOP_OP)

            # 计算偏移指数
            e = math.e
            rindexNS_Slist = []
            for va in rNS_Slist:
                if rALL < epsilon:
                    ri = 0.0
                else:
                    if 0 <= va < epsilon:
                        ri = 0.0
                    elif va > 0:
                        a = float(np.tanh(va)) + 1
                        b = rALL
                        c = - math.log(b + 1, e) * math.log(va + 1, e)
                        ri = 1 - (a ** c)
                    elif va < 0:
                        if float(np.tanh(va)) > -1:
                            a = 1 / (float(np.tanh(va)) + 1)
                        else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                            a = 1E10
                        b = 1 / (rALL)
                        c = - math.log(b + 1, e) * math.log(-va + 1, e)
                        ri = (a ** c) - 1
                rindexNS_Slist.append(ri)

            rindexNOP_OPlist = []
            for va in rNOP_OPlist:
                if ropnop < epsilon:
                    ri = 0.0
                else:
                    if 0 <= va < epsilon:
                        ri = 0.0
                    elif va > 0:
                        a = float(np.tanh(va)) + 1
                        b = ropnop
                        c = - math.log(b + 1, e) * math.log(va + 1, e)
                        ri = 1 - (a ** c)
                    elif va < 0:
                        if float(np.tanh(va)) > -1:
                            a = 1 / (float(np.tanh(va)) + 1)
                        else:  # 如果va的值过于小（负数），tanh由于精确度不够会直接输出-1.0，所以直接将a设置为一个很大的值
                            a = 1E10
                        b = 1 / (ropnop)
                        c = - math.log(b + 1, e) * math.log(-va + 1, e)
                        ri = (a ** c) - 1
                rindexNOP_OPlist.append(ri)




            with open(outputfile910, "w") as output:
                output.write("OffsetDegree_ns/s_" + cladeName + "\n")
                for bb in rNS_Slist:
                    output.write(str(bb) + "\n")

            with open(outputfile1112, "w") as output:
                output.write("OffsetDegree_nop/op_" + cladeName + "\n")
                for bb in rNOP_OPlist:
                    output.write(str(bb) + "\n")

            with open(outputfile17, "w") as output:
                output.write("OffsetIndex_nop/op_" + cladeName + "\n")
                for bb in rindexNOP_OPlist:
                    output.write(str(bb) + "\n")

            with open(outputfile18, "w") as output:
                output.write("OffsetIndex_ns/s_" + cladeName + "\n")
                for bb in rindexNS_Slist:
                    output.write(str(bb) + "\n")

        except ZeroDivisionError:
            print(f"{geneGroupName} has zeroDivision")

        except Exception as e:
            print(f"{geneGroupName} has mistake:  {e}")



print("finished!!!")













