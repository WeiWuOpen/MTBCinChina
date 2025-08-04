

cladeNameList=["4.2"]
geneGroupList = ["MIKGs","essentialGenes","TcellEpitopes","all"]


treeModelType="bsl"
# 分为多少个取样区间(dnds和偏好性)
TimeSet=20

setLimitRate1=0.25
setLimitRate2=0.5


inputrefseqFile = r"...\script_and_testData\testData\H37Rv_r.fasta"
inputCodonBiasFile = r"...\script_and_testData\testData\H37Rv_codon_usage_bias.txt"
inputTi_TvMeanRateFile = r"...\script_and_testData\testData\allstrains_transition_transversion_ratio_mean.txt"
inputTamuraRateFile = r"...\script_and_testData\testData\allstrains_tamura_model_parameters_mean.txt"


import copy
import pandas as pd
import numpy as np
import os
from sklearn.linear_model import LinearRegression
import math
from scipy import stats
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm


for cladeName in tqdm(cladeNameList,desc=treeModelType+" Clade"):
    print(cladeName + " is working #######")

    inputrefseqSNPFile_new = r"...\script_and_testData\testData\\"+cladeName+"_CRMA.fasta"
    if cladeName[0]=="2":
        inputrefseqPOSFile_new = r"...\script_and_testData\testData\allL2_withGeolocInfo_refaMTBC_withprd.pos.txt"
    elif cladeName[0]=="4":
        inputrefseqPOSFile_new = r"...\script_and_testData\testData\allL4_withGeolocInfo_refaMTBC_withprd.pos.txt"
    else:
        print("wrong in inputrefseqPOSFile_new!!!!!!")

    inputFastaFile = r"...\script_and_testData\testData\treetimeResult_refCRMA\\" + cladeName + "_ancestral_results\\ancestral_sequences_r.fasta"
    inputTreeFile = r"...\script_and_testData\testData\treetimeResult_refCRMA\\" + cladeName + "_ancestral_results\\annotated_tree_del.nexus"
    inputBeastTreeFile = r"...\script_and_testData\testData\output\historical_mutation_estimates\\" + cladeName + treeModelType + ".combine.burnin20.repeat3.annotation._nodename.newick"
    inputSNPposFile = r"...\script_and_testData\testData\\" + cladeName + "_refaCRMA_withprd_NOout.pos.txt"



    for geneGroupName in geneGroupList:
        print(cladeName+"_"+geneGroupName + " is working #######")

        inputgenelistFile=r"...\script_and_testData\testData\geneGroup\\" + geneGroupName + ".txt"

        outputFolder4=r"...\script_and_testData\testData\output\historical_mutation_estimates\\"+geneGroupName+"\\"
        outputFolder=outputFolder4+"\\results_refcladeCRMA\\"
        outputFolder2 = outputFolder+treeModelType+"\\"
        outputFolder3 = outputFolder2+str(TimeSet)+"\\"


        if not os.path.exists(outputFolder4):
            os.mkdir(outputFolder4)
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder)
        if not os.path.exists(outputFolder2):
            os.mkdir(outputFolder2)
        if not os.path.exists(outputFolder3):
            os.mkdir(outputFolder3)

        outputFile=outputFolder3 +cladeName +"_"+treeModelType+"_ns_s_allRate_byYBP.txt"
        outputFile2=outputFolder3 + cladeName +"_"+treeModelType+"_dnds_byYBP.txt"
        outputFile3=outputFolder3 +cladeName +"_"+treeModelType+"_codon_preference_byYBP.txt"
        outputFile4=outputFolder3 + cladeName +"_"+treeModelType+"_dndo_byYBP.txt"



        snpPosDict={}
        snpNum=0
        with open(inputSNPposFile,"r") as input:
            for l in input:
                snpNum+=1
                snpPosDict[str(snpNum)] =l.strip()
        del snpNum


        fastaList=[]
        with open(inputFastaFile,"r") as input:
            for l in input:
                fastaList.append(l.strip())
        fastaDict={}
        for r in range(0,len(fastaList)):
            if fastaList[r][0]==">":
                if "/" not in fastaList[r]:
                    fastaDict[fastaList[r].replace(">","")] = fastaList[r+1]
                else: # 把样本名中的"/"后的部分去除
                    fastaDict[fastaList[r].replace(">", "").split("/")[0]] = fastaList[r + 1]



        refseqdict={}
        with open(inputrefseqFile,"r") as input:
            for l in input:
                if l.strip()[0] !=">":
                    refseq=list(l.strip())
        refseqlength=len(refseq)
        for i in range(1,len(refseq)+1):
            refseqdict[int(i)] = refseq[i-1]

        newRefseqdict={}
        newRefPoslist=[]
        with open(inputrefseqSNPFile_new,"r") as input, open(inputrefseqPOSFile_new,"r") as inputpos:
            for l in input:
                if l.strip()[0] != ">":
                    newrefseq = list(l.strip())
            for ll in inputpos:
                newRefPoslist.append(ll.strip())
        for i in range(0,len(newrefseq)):
            newRefseqdict[int(newRefPoslist[i])] = newrefseq[i]

        allLongNewRefseqdict={}
        for k in refseqdict.keys():
            if k in newRefseqdict.keys() and newRefseqdict[k] !="-":
                allLongNewRefseqdict[k] = newRefseqdict[k]
            else:
                allLongNewRefseqdict[k] = refseqdict[k]

        anti_allLongNewRefseqdict={}
        baseDict={"A":"T","T":"A","C":"G","G":"C"}
        for baseindex in allLongNewRefseqdict.keys():
            anti_allLongNewRefseqdict[baseindex] = baseDict[allLongNewRefseqdict[baseindex]]


        genedict={}
        with open(inputgenelistFile,"r") as input:
            for l in input:
                lx=l.strip().split("\t")
                if lx[4] =="+":
                    genedict[lx[0]] =[lx[2],lx[3],"+"]
                else:
                    genedict[lx[0]] = [lx[2], lx[3], "-"]
        print("Gene number: "+str(len(genedict)))



        condonBiasDict={}
        condonBiaslist=[]
        with open (inputCodonBiasFile,"r") as input:
            for l in input:
                lx=l.strip().split()
                if lx !=[] and lx[0] !="Coding":
                    condonBiaslist.append(l.strip().replace("U","T"))
        for ll in condonBiaslist:
            llx=ll.strip().split()
            condonBiasDict[llx[0]] = [llx[1],llx[3]]
            condonBiasDict[llx[6]] = [llx[7],llx[9]]
            condonBiasDict[llx[12]] = [llx[13],llx[15]]
            condonBiasDict[llx[18]] = [llx[19],llx[21]]


        tamuraRateDict={}
        with open(inputTamuraRateFile,"r") as input:
            for ll in input:
                llx=ll.strip().split()
                if llx[0] !="Rate1":
                    tamuraRateDict["rate1"] = float(llx[0])
                    tamuraRateDict["rate2"] = float(llx[1])


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




        #### 判断同义突变是否发生optimal 和 nonoptimal 密码子的替换 #################################################
        def isOptimalCodon(acodon, bcodon): # acodon 是突变前的，bcodon是突变后的
            resultxxx=""
            differentBaseA=""
            differentBaseB = ""
            if condonBiasDict[acodon][0] != condonBiasDict[bcodon][0]:
                print(acodon + "\t" + bcodon + " are not s!!!!")
            else:
                # 这里默认acodon和bcodon只有一个碱基的差异
                for iu in range(0,3):
                    if acodon[iu] != bcodon[iu]:
                        differentBaseA=acodon[iu]
                        differentBaseB=bcodon[iu]
                if  differentBaseA =="" or  differentBaseB=="":
                    print("wrong!!! in isOptimalCodon function!!!")
                if baseDict[differentBaseA] == differentBaseB:
                    resultxxx="normal"
                else:
                    if float(condonBiasDict[acodon][1]) > float(condonBiasDict[bcodon][1]):
                        resultxxx="nonoptimal"
                    else:
                        resultxxx="optimal"

            return resultxxx
        ###########################################################################################

        ########## 目标基因群上全部待选ns，s位点数目与待选optimal,nonoptimal位点数目 ####################################
        NSsiteNum=0
        SsiteNum=0
        optimalNum=0
        nonoptimalNum=0
        normalNum=0
        allGeneLength=0
        ggg=0
        for gg in genedict.keys():
            ggg+=1
            for pp in range(int(genedict[gg][0]), int(genedict[gg][1]) + 1, 3):
                allGeneLength += 3
                if genedict[gg][2] =="+":
                    codon = allLongNewRefseqdict[pp]+allLongNewRefseqdict[pp+1]+allLongNewRefseqdict[pp+2]
                elif genedict[gg][2] =="-":
                    codon = anti_allLongNewRefseqdict[pp + 2] + anti_allLongNewRefseqdict[pp + 1] + anti_allLongNewRefseqdict[pp]
                else:
                    print("WRONG IN genedict!!!!!!!")

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

                    optimalNum += opxx
                    nonoptimalNum += nopxx
                    SsiteNum += sxx
                    NSsiteNum += nsxx
                    normalNum += normalxx

                    del yi
                    del opxx
                    del nopxx
                    del sxx
                    del nsxx
                    del normalxx
                    del aminoDict

            del codon

        # print("all target gene: ")
        # print(allGeneLength)
        # print(NSsiteNum)
        # print(SsiteNum)
        # print(optimalNum)
        # print(nonoptimalNum)
        # print(normalNum)
        # print(optimalNum+normalNum+nonoptimalNum)
        # print("#############################")




        with open(inputTreeFile,"r") as input:
            for l in input:
                if l.strip().split()[0] =="Tree":
                    tree=l.strip().split()[1].replace("tree1=","").replace("'", "")  # 去除有些树的样本名外的引号

        with open(inputBeastTreeFile,"r") as input:
            btree=input.readline().strip().replace("'", "")  # 去除有些树的样本名外的引号



        ######################## 计算根节点的bp ##########################################################
        # 将数据最后一个样本点的所有枝长相加
        ################# 自定义判断判断是否为数字函数(isdigit函数对小数无效）
        def is_number(s):
            try:  # 如果能运⾏ float(s) 语句，返回 True（字符串 s 是浮点数）
                float(s)
                return True
            except ValueError:  # ValueError 为 Python 的⼀种标准异常，表⽰"传⼊⽆效的参数"
                pass  # 如果引发了 ValueError 这种异常，不做任何事情（pass：不做任何事情，⼀般⽤做占位语句）
            try:
                import unicodedata  # 处理 ASCII 码的包
                unicodedata.numeric(s)  # 把⼀个表⽰数字的字符串转换为浮点数返回的函数
                return True
            except (TypeError, ValueError):
                pass
                return False
        #####################
        # 去除tree中存在的node名
        indexList=[]
        for c in range(0,len(btree)):
            indexList.append(c)
        for cc in range(0,len(btree)):
            if btree[cc]==")" and btree[cc+1]=="N": # 这里node名统一是以N（Nodexxx)开头的，若有其他命名法，需要更改
                for uu in range(cc+1,len(btree)):
                    if btree[uu] !=":":
                        indexList.remove(uu)
                    else:
                        break
        btree2=""
        for ccc in indexList:
            btree2 += btree[ccc]
        btree2 += ";"

        rootBP=0
        splitTree=btree2.split("):")
        splitTree[-1] = splitTree[-1].replace(");", "")  # 去除掉最后一个字符串中的”);“符号

        for tt in range(len(splitTree) - 1, -1, -1):  # 从后往前进行遍历
            if (is_number(splitTree[tt])):  # 当数据只有数字时
                t = '%.4f' % float(splitTree[tt])  # 四舍五入保留四位小数
                rootBP += float(t)
                continue
            else:  # 最后一个枝长包含样本名的处理
                t = '%.4f' % float(splitTree[tt].split(":")[-1])
                rootBP += float(t)
                break
        rootBP = float('%.4f' % rootBP)  # 处理浮点数运算的精度丢失问题
        rootBP =  rootBP

        # 修正存在负数BP时候的rootBP
        treeCutx1 = btree2.split(":")
        treeCutx2 = [j for i in treeCutx1 for j in i.split(",")]
        treeCutx3 = [j for i in treeCutx2 for j in i.split(")")]
        treeCutx4 = [j for i in treeCutx3 for j in i.split("(")]
        absnegativeBP=0.0
        for tc in treeCutx4:
            if (is_number(tc)):  # 当数据只有数字时
                if float(tc) <0: # 寻找负数BP
                   if abs(float(tc)) > absnegativeBP:
                       absnegativeBP = abs(float(tc))  # 寻找绝对值最大的负数BP
        rootBP = rootBP + absnegativeBP

        print("rootBP :  {}".format(str(rootBP)))
        #######################################################################################



        # 找到所有的NODE
        nodeDict={}
        for i in range(0,len(tree)):
            if tree[i]=="(":
                sampleList=[]
                obindex=i
                obnum=1 # 前括号数量
                fbnum=0 # 后括号数量
                for q in range(i+1,len(tree)):
                    if tree[q]=="(":
                        obnum +=1
                    elif tree[q]==")":
                        fbnum +=1
                        if fbnum == obnum:
                            fbindex=q
                            break
                    else:
                        continue
                for d in range(fbindex+1,len(tree)):
                    if tree[d] == ":" or tree[d]==";":
                        nodeName=tree[fbindex+1:d]
                        break
                    else:
                        continue
                innerCut=tree[obindex+1:fbindex]
                innerCutx1=innerCut.split(":")
                innerCutx2 = [j for i in innerCutx1 for j in i.split(",")]
                innerCutx3 = [j for i in innerCutx2 for j in i.split(")")]
                innerCutx4 = [j for i in innerCutx3 for j in i.split("(")]
                for sp in innerCutx4:
                    if "RR" in sp or "SC" in sp:  # 样本命名中的开头的ERR、SRR、SC，若有其他命名方法，需要修改
                        sampleList.append(sp.split("/")[0])
                nodeDict[nodeName]=[sampleList]
                del sampleList
                del obnum
                del obindex
                del fbindex
                del fbnum
                del innerCut
                del innerCutx1
                del innerCutx2
                del innerCutx3
                del innerCutx4
                del q
                del d
                del sp
            else:
                continue


        beastNodeDict={}
        for i in range(0,len(btree)):
            if btree[i]=="(":
                sampleList=[]
                obindex=i
                obnum=1 # 前括号数量
                fbnum=0 # 后括号数量
                for q in range(i+1,len(btree)):
                    if btree[q]=="(":
                        obnum +=1
                    elif btree[q]==")":
                        fbnum +=1
                        if fbnum == obnum:
                            fbindex=q
                            break
                    else:
                        continue
                for d in range(fbindex+1,len(btree)):
                    if btree[d] == ":" or btree[d]==";":
                        nodeName=btree[fbindex+1:d]
                        colonIndex=d
                        break
                    else:
                        continue
                if btree[colonIndex] ==":":
                    for v in range(colonIndex+1,len(btree)):
                        if btree[v]=="," or btree[v]==")":
                            nodeLength=float(btree[colonIndex+1:v])
                            break
                        else:
                            continue
                elif btree[colonIndex] ==";": # 最外层的根节点
                    nodeLength=0.0
                    v=""
                else:
                    print("wrong!!! in colonIndex")


                innerCut=btree[obindex+1:fbindex]
                innerCutx1=innerCut.split(":")
                innerCutx2 = [j for i in innerCutx1 for j in i.split(",")]
                innerCutx3 = [j for i in innerCutx2 for j in i.split(")")]
                innerCutx4 = [j for i in innerCutx3 for j in i.split("(")]
                for sp in innerCutx4:
                    if "RR" in sp or "SC" in sp:  ######### 样本命名中的开头的ERR、SRR、SC，若有其他命名方法，需要修改
                        sampleList.append(sp.split("/")[0])
                beastNodeDict[nodeName]=[sampleList,nodeLength]
                del sampleList
                del obnum
                del obindex
                del fbindex
                del fbnum
                del innerCut
                del innerCutx1
                del innerCutx2
                del innerCutx3
                del innerCutx4
                del q
                del d
                del sp
                del nodeLength
                del nodeName
                del colonIndex
                del v
            else:
                continue

        ## 寻找父节点
        nodeDictf=copy.deepcopy(nodeDict)
        for n in nodeDictf.keys():
            fnodeName="root"
            sampleNum=999999999999999
            for b in nodeDict.keys():
                if b !=n:
                    is_subset=all(elem in nodeDict[b][0] for elem in nodeDictf[n][0]) # 判断后者是否是前者的子集
                    if is_subset:
                        if len(nodeDict[b][0]) < sampleNum: # 寻找包含样本最少的node，就是父节点
                            sampleNum =  len(nodeDict[b][0])
                            fnodeName = b
            nodeDictf[n].append(fnodeName)
            del fnodeName
            del b
            del sampleNum


        beastNodeDictf=copy.deepcopy(beastNodeDict)
        for n in beastNodeDictf.keys():
            fnodeName="root"
            sampleNum=999999999999999
            for b in beastNodeDict.keys():
                if b !=n:
                    is_subset=all(elem in beastNodeDict[b][0] for elem in beastNodeDictf[n][0]) # 判断后者是否是前者的子集
                    if is_subset:
                        if len(beastNodeDict[b][0]) < sampleNum: # 寻找包含样本最少的node，就是父节点
                            sampleNum =  len(beastNodeDict[b][0])
                            fnodeName = b
            beastNodeDictf[n].append(fnodeName)
            del fnodeName
            del b
            del sampleNum


        # 把nodelength转化成nodeBP
        beastNodeDictfr = copy.deepcopy(beastNodeDictf)
        for zz in beastNodeDictf.keys():
            dnodeBP=0
            if beastNodeDictf[zz][2] !="root":
                tmpNode=zz
                dnodeBP += beastNodeDictf[tmpNode][1]
                while True:
                    tmpNode = beastNodeDictf[tmpNode][2]
                    if tmpNode=="root":
                        break
                    else:
                        dnodeBP +=  beastNodeDictf[tmpNode][1]
            else:
                dnodeBP=0
            nodeNP=rootBP - dnodeBP
            beastNodeDictfr[zz][1] = nodeNP
            del dnodeBP
            del nodeNP


        ## match两种树的node名
        beastNodeDictfm=copy.deepcopy(beastNodeDictfr)
        for ee in beastNodeDictf.keys():
            for ww in nodeDictf.keys():
                set1=set(nodeDictf[ww][0])
                set2=set(beastNodeDictf[ee][0])
                intersectionSet= set1 & set2 # 两者的交集
                if len(intersectionSet) == len(set1) and len(intersectionSet) == len(set2): ## node包含的样本(不重复样本)完全相同
                    beastNodeDictfm[ee].append(ww)
                    break
                else:
                    continue
        del ee
        del ww
        del set1
        del set2
        del intersectionSet





        ### filter node  #############
        # 检查无法对应的node，若有，删除
        beastNodeDictfm_del1=copy.deepcopy(beastNodeDictfm)
        for pp in beastNodeDictfm.keys():
            if len(beastNodeDictfm[pp]) < 4:
                # print(pp+"\t"+str(len(beastNodeDictfm[pp])))
                del beastNodeDictfm_del1[pp]
        del pp
        print("#########################")


        # 剔除负数length的node
        negativeBPnodeList=[]
        beastNodeDictfm_del2=copy.deepcopy(beastNodeDictfm_del1)
        for ppp in beastNodeDictfm_del1.keys():
            if beastNodeDictf[ppp][1] <0:  # 注意使用的是beastNodeDictf，即判断nodeLength是负数，而不是nodeBP
                # print(ppp+"\t"+str(beastNodeDictf[ppp][1]))
                negativeBPnodeList.append(ppp)
                del beastNodeDictfm_del2[ppp]
        del ppp
        print("#######################")

        # 剔除负数length的子节点,只去除直接的一层
        beastNodeDictfm_del3 = copy.deepcopy(beastNodeDictfm_del2)
        for pppp in beastNodeDictfm_del2.keys():
            if beastNodeDictfm_del2[pppp][2] in negativeBPnodeList:
                # print(pppp)
                del beastNodeDictfm_del3[pppp]
        del pppp
        print("#######################")


        # 剔除与父节点时间距离过远的node
        beastNodeDictfm_del4 = copy.deepcopy(beastNodeDictfm_del3)
        for p5 in beastNodeDictfm_del3.keys():
            if beastNodeDictfm_del3[p5][2] !="root":
                nodeBP=beastNodeDictfm_del3[p5][1]
                fnodeBP=beastNodeDictfm[ beastNodeDictfm_del3[p5][2] ][1]
                if nodeBP/fnodeBP < setLimitRate1 : # 自设定倍数
                    # print(p5+"\t"+str(nodeBP)+"\t"+str(fnodeBP))
                    del beastNodeDictfm_del4[p5]
            else:
                nodeNP=""
                fnodeBP=""
                continue
        del p5
        del nodeNP
        del fnodeBP
        print("########################")

        print("after filter node: "+ str(len(beastNodeDictfm_del4))+"/"+str(len(beastNodeDictfm))+" rate:"+str(len(beastNodeDictfm_del4)/len(beastNodeDictfm)))


        snpDict={}
        for kk in beastNodeDictfm_del4.keys():
            if beastNodeDictfm_del4[kk][2] !="root":
                nodeNamet=beastNodeDictfm_del4[kk][-1]
                fnodeNamet=beastNodeDictfm[ beastNodeDictfm_del4[kk][2]  ][-1] # 注意这里使用的beastNodeDictfm，因为beastNodeDictfm_del4有可能把对应的父节点剔除了
                if nodeNamet in fastaDict.keys() and fnodeNamet in fastaDict.keys():
                    snpList=[]
                    for p in range(0,len(fastaDict[nodeNamet])):
                        if fastaDict[ nodeNamet ][p] != fastaDict[ fnodeNamet ][p]:
                            if fastaDict[ nodeNamet ][p] !="-" and fastaDict[ fnodeNamet ][p]  !="-":
                                if fastaDict[nodeNamet][p] != "N" and fastaDict[fnodeNamet][p] != "N":
                                    snpList.append(str(p+1)+"_"+ fastaDict[ fnodeNamet ][p]+"_"+ fastaDict[ nodeNamet ][p])
                    if len(snpList) > 0:
                        snpDict[kk] =snpList
                    del snpList
                else:
                    print(kk+" NOT match to fasta Node, del!!!!!!!")
        del kk
        del p
        del nodeNamet

        #####################################################
        ##### 从末端node到tip ################################
        finalNodeList = []
        for node in beastNodeDictfm_del4.keys():
            if len(beastNodeDictfm_del4[node][0]) == 2:  ## 当且仅当只包含两个样本的node是末端node（二叉树）
                finalNodeList.append(node)
        del node

        # 剔除与tip距离过远的末端node
        finalNodeList_2=[]
        sssn=0
        for ss in finalNodeList:
            nodeBP= int(beastNodeDictfm_del4[ss][1])
            if nodeBP / rootBP < setLimitRate2:  # 自设定倍数
                finalNodeList_2.append(ss)
            else:
                sssn+=1
        print(str(sssn)+"/"+str(len(finalNodeList))+" final node is filtered !!!!!!!!!!!!!!!!")


        snpDict_tip={}
        for mm in finalNodeList_2:
            mm1=mm+"_1"
            mm2=mm+"_2"
            nodeNamet = beastNodeDictfm_del4[mm][-1]
            tipName1=beastNodeDictfm_del4[mm][0][0]
            tipName2 = beastNodeDictfm_del4[mm][0][1]
            if nodeNamet in fastaDict.keys() and tipName1 in fastaDict.keys() and tipName2 in fastaDict.keys():
                snpList1 = []
                snpList2 = []
                for pp in range(0, len(fastaDict[nodeNamet])):
                    if fastaDict[nodeNamet][pp] != fastaDict[tipName1][pp]:
                        if fastaDict[nodeNamet][pp] != "-" and fastaDict[tipName1][pp] != "-":
                            if fastaDict[nodeNamet][pp] != "N" and fastaDict[tipName1][pp] != "N":
                                snpList1.append(str(pp + 1) + "_" + fastaDict[nodeNamet][pp] + "_" + fastaDict[tipName1][pp])

                for pp in range(0, len(fastaDict[nodeNamet])):
                    if fastaDict[nodeNamet][pp] != fastaDict[tipName2][pp]:
                        if fastaDict[nodeNamet][pp] != "-" and fastaDict[tipName2][pp] != "-":
                            if fastaDict[nodeNamet][pp] != "N" and fastaDict[tipName2][pp] != "N":
                                snpList2.append(str(pp + 1) + "_" + fastaDict[nodeNamet][pp] + "_" + fastaDict[tipName2][pp])

                if len(snpList1) > 0:
                    snpDict_tip[mm1] = snpList1
                if len(snpList2) > 0:
                    snpDict_tip[mm2] = snpList2
                del snpList1
                del snpList2
            else:
                print(mm + " NOT match to fasta Node or tip, del!!!!!!!")

        #######################################################



        snpDict2={}
        for kk in snpDict.keys():
            snpDict2[kk]=[]
            for iii in snpDict[kk]:
                iiix=iii.split("_")
                snpDict2[kk].append(snpPosDict[iiix[0]]+"_"+iiix[1]+"_"+iiix[2])
        # del kk
        # del iii
        for kk in snpDict_tip.keys():
            snpDict2[kk]=[]
            for iii in snpDict_tip[kk]:
                iiix=iii.split("_")
                snpDict2[kk].append(snpPosDict[iiix[0]]+"_"+iiix[1]+"_"+iiix[2])
        # del kk
        # del iii


        snpDict3={}
        for kk in snpDict2.keys():
            snpDict3[kk]=[]
            for y in snpDict2[kk]:
                #### 判断是否在CDS内 ####
                geneList=[]
                for gg in sorted(genedict.keys()):
                    if int(genedict[gg][0])  <= int(y.split("_")[0]) <= int(genedict[gg][1]):
                        geneList.append(gg)
                        break
                    else:
                        continue
                if len(geneList) >0:
                    snpDict3[kk].append(y+"_"+geneList[0])
                else:
                    continue


        dndsDict={}
        biasDict={}
        for kk in snpDict3.keys():
            nsNum=0
            sNum=0
            biasDict[kk] = [0, 0, 0]
            for p in snpDict3[kk]:
                        if genedict[p.split("_")[-1]][-1]=="+":
                            if (int(p.split("_")[0])-int(genedict[p.split("_")[-1]][0])+1)%3 ==1:
                                aminoAcid=p.split("_")[2]+ allLongNewRefseqdict[int(p.split("_")[0]) + 1]+ allLongNewRefseqdict[int(p.split("_")[0]) + 2]
                                aminoAcidREF=p.split("_")[1]+ allLongNewRefseqdict[int(p.split("_")[0]) + 1]+ allLongNewRefseqdict[int(p.split("_")[0]) + 2]

                            elif (int(p.split("_")[0])-int(genedict[p.split("_")[-1]][0])+1)%3 ==2:
                                aminoAcid=allLongNewRefseqdict[int(p.split("_")[0]) - 1]+ p.split("_")[2]+ allLongNewRefseqdict[int(p.split("_")[0]) + 1]
                                aminoAcidREF=allLongNewRefseqdict[int(p.split("_")[0]) - 1] + p.split("_")[1]+ allLongNewRefseqdict[int(p.split("_")[0]) + 1]

                            else:
                                aminoAcid = allLongNewRefseqdict[int(p.split("_")[0]) - 2]+ allLongNewRefseqdict[int(p.split("_")[0]) - 1]+ p.split("_")[2]
                                aminoAcidREF = allLongNewRefseqdict[int(p.split("_")[0]) - 2]+ allLongNewRefseqdict[int(p.split("_")[0]) - 1]+ p.split("_")[1]

                        elif genedict[p.split("_")[-1]][-1]=="-":
                            if (int(genedict[p.split("_")[-1]][1]) - int(p.split("_")[0]) + 1) < 0:
                                print("gene end pos and snp pos is wrong !!!!!!!!!!")
                            else:
                                if (int(genedict[p.split("_")[-1]][1]) -int(p.split("_")[0]) + 1) % 3 == 1:
                                    aminoAcid = baseDict[p.split("_")[2]]+ anti_allLongNewRefseqdict[int(p.split("_")[0]) - 1]+ anti_allLongNewRefseqdict[int(p.split("_")[0]) - 2]
                                    aminoAcidREF = baseDict[p.split("_")[1]]+ anti_allLongNewRefseqdict[int(p.split("_")[0]) - 1]+ anti_allLongNewRefseqdict[int(p.split("_")[0]) - 2]

                                elif (int(genedict[p.split("_")[-1]][1]) -int(p.split("_")[0]) + 1) % 3 == 2:
                                    aminoAcid = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1]+ baseDict[p.split("_")[2]]+ anti_allLongNewRefseqdict[int(p.split("_")[0]) - 1]
                                    aminoAcidREF = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1]+ baseDict[p.split("_")[1]]+ anti_allLongNewRefseqdict[int(p.split("_")[0]) - 1]

                                elif (int(genedict[p.split("_")[-1]][1]) - int(p.split("_")[0]) + 1) % 3 == 0:
                                    aminoAcid = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 2]+ anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1]+ baseDict[p.split("_")[2]]
                                    aminoAcidREF = anti_allLongNewRefseqdict[int(p.split("_")[0]) + 2]+ anti_allLongNewRefseqdict[int(p.split("_")[0]) + 1]+ baseDict[p.split("_")[1]]
                                else:
                                    print("codon is wrong !!!!!!!!!")
                        else:
                            print("wrong!!!!")

                        if condonBiasDict[aminoAcid][0] == condonBiasDict[aminoAcidREF][0]:
                            sNum+=1
                            if isOptimalCodon(aminoAcidREF,aminoAcid) =="nonoptimal":
                                biasDict[kk][0] +=1
                            elif isOptimalCodon(aminoAcidREF,aminoAcid) =="optimal":
                                biasDict[kk][1] +=1
                            elif isOptimalCodon(aminoAcidREF,aminoAcid) =="normal":
                                biasDict[kk][2] +=1
                            else:
                                print("Something wrong in codon preference!!!!")
                        else:
                            nsNum += 1
            dndsDict[kk]=[nsNum,sNum]



        ###################################################################################
        # Rate ns s all
        resultDict={}
        for cc in dndsDict.keys():
            if "_" not in cc:
                nodeBPz=beastNodeDictfm[cc][1]
                fnodeBPz=beastNodeDictfm[ beastNodeDictfm[cc][2] ][1]
                resultDict[cc] = [fnodeBPz,nodeBPz,dndsDict[cc][0],dndsDict[cc][1]]
            else: # 这是末端node到tip
                nodeBPz=0
                fnodeBPz=beastNodeDictfm[cc.split("_")[0]][1]
                resultDict[cc] = [fnodeBPz,nodeBPz,dndsDict[cc][0],dndsDict[cc][1]]

        del cc
        del nodeBPz
        del fnodeBPz


        resultDict2={}
        for rr in resultDict.keys():
            if resultDict[rr][1] != resultDict[rr][0]:
                nsRate= (resultDict[rr][2] /allGeneLength)
                sRate= (resultDict[rr][3] / allGeneLength)
                allRate= (resultDict[rr][2]+resultDict[rr][3]) / allGeneLength
                resultDict2[rr] = [resultDict[rr][0],resultDict[rr][1],nsRate,sRate,allRate]
        del rr

        byYearResultDict_2={}
        num=0
        for nn in resultDict2.keys():
            num+=1
            for t in range( int(resultDict2[nn][1]) , int(resultDict2[nn][0])+1 ):
                if t not in byYearResultDict_2.keys():
                    byYearResultDict_2[t]=[[resultDict2[nn][2]/(int(resultDict2[nn][0])-int(resultDict2[nn][1])+1)],\
                                         [resultDict2[nn][3]/(int(resultDict2[nn][0])-int(resultDict2[nn][1])+1)],\
                                         [resultDict2[nn][4]/(int(resultDict2[nn][0])-int(resultDict2[nn][1])+1)] ]
                else:
                    byYearResultDict_2[t][0] .append( resultDict2[nn][2] / (int(resultDict2[nn][0]) - int(resultDict2[nn][1]) + 1) )
                    byYearResultDict_2[t][1] .append( resultDict2[nn][3] / (int(resultDict2[nn][0]) - int(resultDict2[nn][1]) + 1) )
                    byYearResultDict_2[t][2] .append( resultDict2[nn][4] / (int(resultDict2[nn][0]) - int(resultDict2[nn][1]) + 1) )


        gg = 0
        byYearResultDict = {}
        for kk in byYearResultDict_2.keys():
            if len(byYearResultDict_2[kk][0]) == 0 or len(byYearResultDict_2[kk][1]) == 0 or len(byYearResultDict_2[kk][2]) == 0:
                gg += 1
            else:
                byYearResultDict[kk] = byYearResultDict_2[kk]
        del kk

        print("no sample (bp) in ns s all , remove!!  rate: " + str(gg) + " / " + str(len(byYearResultDict_2)))



        ## 利用百分位数方法去除异常值 #########
        # byYearResultDict3 = {}
        # for kk in byYearResultDict.keys():
        #     byYearResultDict3[kk] = [[], [], []]
        #     outlinerIndexList = []
        #     for rr in range(0, len(byYearResultDict[kk])):
        #         iqr = np.percentile(byYearResultDict[kk][rr], 75) - np.percentile(byYearResultDict[kk][rr], 25)
        #         uplimit = np.percentile(byYearResultDict[kk][rr], 75) + 1.5 * iqr
        #         downlimit = np.percentile(byYearResultDict[kk][rr], 25) - 1.5 * iqr
        #         iiIndex = -1
        #         for ii in byYearResultDict[kk][rr]:
        #             iiIndex += 1
        #             if ii < downlimit or ii > uplimit:
        #                 outlinerIndexList.append(iiIndex)
        #     outlinerIndexList = np.unique(outlinerIndexList)
        #     for rr in range(0, len(byYearResultDict[kk])):
        #         ppIndex = -1
        #         for pp in byYearResultDict[kk][rr]:
        #             ppIndex += 1
        #             if ppIndex not in outlinerIndexList:
        #                 byYearResultDict3[kk][rr].append(pp)
        # del kk
        # del ii
        # del pp
        # del rr
        # del iiIndex
        # del ppIndex
        # del outlinerIndexList

        byYearResultDict3=byYearResultDict

        byYearResultDict3_2 = {}
        # 计算每个年上的均值和总体标准差
        for bbb in byYearResultDict3.keys():
            byYearResultDict3_2[bbb] = [[np.mean(byYearResultDict3[bbb][0]), np.std(byYearResultDict3[bbb][0])], \
                                      [np.mean(byYearResultDict3[bbb][1]), np.std(byYearResultDict3[bbb][1])], \
                                      [np.mean(byYearResultDict3[bbb][2]), np.std(byYearResultDict3[bbb][2])], ]


        byYearResultDict4 = {}
        ttt = 0
        for ccc in byYearResultDict3_2.keys():
            byYearResultDict4[ccc] = byYearResultDict3_2[ccc]


        #############################################
        ## 偏好性
        biasDict2={}
        for kkkk in biasDict.keys():
            if "_" not in kkkk:
                nodeBPy = beastNodeDictfm[kkkk][1]
                fnodeBPy = beastNodeDictfm[ beastNodeDictfm[kkkk][2] ][1]
                biasDict2[kkkk] = [fnodeBPy,nodeBPy, biasDict[kkkk][0],biasDict[kkkk][1],biasDict[kkkk][2]]
            else: # 这是末端node到tip
                nodeBPy=0
                fnodeBPy = beastNodeDictfm[kkkk.split("_")[0]][1]
                biasDict2[kkkk] = [fnodeBPy, nodeBPy, biasDict[kkkk][0], biasDict[kkkk][1], biasDict[kkkk][2]]
        del kkkk
        del nodeBPy
        del fnodeBPy


        biasDict3={}
        for rr in biasDict2.keys():
            if biasDict2[rr][1] != biasDict2[rr][0]:
                nonoptimalRate= biasDict2[rr][2] / allGeneLength
                optimalRate= biasDict2[rr][3] / allGeneLength
                normalRate=  biasDict2[rr][4] / allGeneLength
                biasDict3[rr] = [biasDict2[rr][0],biasDict2[rr][1],nonoptimalRate,optimalRate,normalRate]


        byYearBiasDict_2={}
        for nn in biasDict3.keys():
            for t in range( int(biasDict3[nn][1]) , int(biasDict3[nn][0])+1):
                if t not in byYearBiasDict_2.keys():
                    byYearBiasDict_2[t]=[[biasDict3[nn][2]/(int(biasDict3[nn][0])-int(biasDict3[nn][1])+1)],\
                                         [biasDict3[nn][3]/(int(biasDict3[nn][0])-int(biasDict3[nn][1])+1)],\
                                         [biasDict3[nn][4]/(int(biasDict3[nn][0])-int(biasDict3[nn][1])+1)] ]
                else:
                    byYearBiasDict_2[t][0] .append( biasDict3[nn][2] / (int(biasDict3[nn][0])-int(biasDict3[nn][1])+    1) )
                    byYearBiasDict_2[t][1] .append( biasDict3[nn][3] / (int(biasDict3[nn][0]) - int(biasDict3[nn][1]) + 1) )
                    byYearBiasDict_2[t][2] .append( biasDict3[nn][4] / (int(biasDict3[nn][0]) - int(biasDict3[nn][1]) + 1) )


        gg=0
        byYearBiasDict={}
        for kk in byYearBiasDict_2.keys():
            if len(byYearBiasDict_2[kk][0]) ==0 or len(byYearBiasDict_2[kk][1])==0 or len(byYearBiasDict_2[kk][2])==0:
                gg+=1
            else:
                byYearBiasDict[kk] = byYearBiasDict_2[kk]

        print("no sample (ybp) in bias , remove!!  rate: "+str(gg)+" / "+str(len(byYearBiasDict_2)))




        ## 利用百分位数方法去除异常值 #########
        # byYearBiasDict3 = {}
        # for kk in byYearBiasDict.keys():
        #     byYearBiasDict3[kk] = [[], [], []]
        #     outlinerIndexList = []
        #     for rr in range(0, len(byYearBiasDict[kk])):
        #         iqr = np.percentile(byYearBiasDict[kk][rr], 75) - np.percentile(byYearBiasDict[kk][rr], 25)
        #         uplimit = np.percentile(byYearBiasDict[kk][rr], 75) + 1.5 * iqr
        #         downlimit = np.percentile(byYearBiasDict[kk][rr], 25) - 1.5 * iqr
        #         iiIndex = -1
        #         for ii in byYearBiasDict[kk][rr]:
        #             iiIndex += 1
        #             if ii < downlimit or ii > uplimit:
        #                 outlinerIndexList.append(iiIndex)
        #     outlinerIndexList = np.unique(outlinerIndexList)
        #     for rr in range(0, len(byYearBiasDict[kk])):
        #         ppIndex = -1
        #         for pp in byYearBiasDict[kk][rr]:
        #             ppIndex += 1
        #             if ppIndex not in outlinerIndexList:
        #                 byYearBiasDict3[kk][rr].append(pp)
        # del kk
        # del ii
        # del pp
        # del rr
        # del iiIndex
        # del ppIndex
        # del outlinerIndexList

        byYearBiasDict3=byYearBiasDict

        byYearBiasDict3_2={}
        # 计算均值和总体标准差
        for bbb in byYearBiasDict3.keys():
            byYearBiasDict3_2[bbb]=[[np.mean(byYearBiasDict3[bbb][0]),np.std(byYearBiasDict3[bbb][0])],\
                                     [np.mean(byYearBiasDict3[bbb][1]),np.std(byYearBiasDict3[bbb][1])],\
                                     [np.mean(byYearBiasDict3[bbb][2]),np.std(byYearBiasDict3[bbb][2])],]

        byYearBiasDict4={}
        ttt=0
        for ccc in byYearBiasDict3_2.keys():
            byYearBiasDict4[ccc] = byYearBiasDict3_2[ccc]




        ############ dnds和dndo ################################
        dndsResultDict={}
        dndoResultDict={}
        timecc= rootBP / TimeSet
        for ti in range(0,TimeSet):
            ## dnds #######################
            nsMeanList=[]
            sMeanList=[]
            for kk in byYearResultDict4.keys():
                if ti*timecc <= kk < (ti+1)*timecc:
                    # 换分母，从基因全长换成潜在的ns/s数目
                    nsMeanList.append(byYearResultDict4[kk][0][0] *allGeneLength/NSsiteNum )
                    sMeanList.append(byYearResultDict4[kk][1][0] *allGeneLength/SsiteNum )


            if len(nsMeanList) == 0 and len(sMeanList) == 0:
                print(str(ti * timecc) + "\t" + str((ti + 1) * timecc) + " this time period is empty,remove (dnds) !!!")
            else:
                ## 设定时间段内所有NS和S各自的总和来算dn/ds
                nsmeanSum = sum(nsMeanList)
                smeanSum = sum(sMeanList)
                allSum=nsmeanSum+smeanSum
                if allSum ==0: # 没有发生任何突变，设置dnds=0
                    dndsMean = 0
                elif smeanSum ==0 and nsmeanSum !=0 :  # 没有同义突变发生，设置dnds=2.5（便于画图）
                    dndsMean=2.5
                else:
                    dndsMean = nsmeanSum / smeanSum
                    if dndsMean >2: # 人为设置为2，便于画图
                        dndsMean=2
                for time in range(int(ti * timecc) + 1, int((ti + 1) * timecc)):
                    dndsResultDict[time] = dndsMean
                del dndsMean


            ## 偏好性 ###############################
            nonoptimalMeanList=[]
            optimalMeanList=[]
            for kk in byYearBiasDict4.keys():
                if ti*timecc <= kk < (ti+1)*timecc:
                    # 换分母，从基因全长替换成潜在的optimal/nonoptimal数目
                    nonoptimalMeanList.append(byYearBiasDict4[kk][0][0] *allGeneLength/nonoptimalNum)
                    optimalMeanList.append(byYearBiasDict4[kk][1][0] *allGeneLength/optimalNum)


            if len(nonoptimalMeanList) ==0 and len(optimalMeanList) ==0:
                print(str(ti*timecc)+"\t"+str((ti+1)*timecc)+" this time period is empty,remove (偏好性) !!!")
            else:
                nonoptimalmeanSum = sum(nonoptimalMeanList)
                optimalmeanSum = sum(optimalMeanList)
                if optimalmeanSum==0 and nonoptimalmeanSum==0:  # 没有发生这两种突变，设置dndo=0
                    dndoMean=0
                else:
                    if optimalmeanSum==0: # 没有发生向optimal密码子的突变
                       dndoMean=4.5  ## 人为设置为4.5，便于画图
                    else:
                        dndoMean = nonoptimalmeanSum / optimalmeanSum
                        if dndoMean >4: # 人为设置为4，便于画图
                            dndoMean=4
                for time in range(int(ti*timecc+1) , int((ti+1)*timecc)):
                    dndoResultDict[time] = dndoMean
                del dndoMean


        dndsResultDict2={}
        for kkk in dndsResultDict.keys():
            # 去除nan值
            if math.isnan(dndsResultDict[kkk]):
                print(str(kkk) +" ybp dnds is nan, to remove!!!")
            else:
                dndsResultDict2[kkk]=dndsResultDict[kkk]


        dndoResultDict2={}
        for kkk in dndoResultDict.keys():
            # 去除nan值
            if math.isnan(dndoResultDict[kkk]):
                print(str(kkk) +" ybp dndo is nan, to remove!!!")
            else:
                dndoResultDict2[kkk]=dndoResultDict[kkk]


        with open(outputFile,"w") as output:
            output.write(cladeName+"_YBP"+"\t"+cladeName+"_ns"+"\t"+cladeName+"_nsStd"+"\t"+cladeName+"_s"+"\t"+cladeName+"_sStd"\
                        +"\t"+cladeName+"_all"+"\t"+cladeName+"_allStd"+"\n")
            for kkk in sorted(byYearResultDict4.keys()):
                output.write(str(kkk)+"\t"+str(byYearResultDict4[kkk][0][0])+"\t"+str(byYearResultDict4[kkk][0][1])+"\t"+str(byYearResultDict4[kkk][1][0])+"\t"+str(byYearResultDict4[kkk][1][1])\
                             +"\t"+str(byYearResultDict4[kkk][2][0])+"\t"+str(byYearResultDict4[kkk][2][1])+"\n")

        with open(outputFile2,"w") as output:
            output.write(cladeName+"_meanYBP"+"\t"+cladeName+"_dnds"+"\n")
            for kkkk in sorted(dndsResultDict2.keys()):
                output.write(str(kkkk)+"\t"+str(dndsResultDict2[kkkk])+"\n")


        with open(outputFile3,"w") as output:
            output.write(cladeName+"_YBP"+"\t"+cladeName+"_nonoptimal"+"\t"+cladeName+"_nonoptimalStd"+"\t"+cladeName+"_optimal"+"\t"+cladeName+"_optimalStd"\
                         +"\t"+cladeName+"_normal"+"\t"+cladeName+"_normalStd"+"\n")
            for kkk in sorted(byYearBiasDict4.keys()):
                output.write(str(kkk)+"\t"+str(byYearBiasDict4[kkk][0][0])+"\t"+str(byYearBiasDict4[kkk][0][1])+"\t"+str(byYearBiasDict4[kkk][1][0])+"\t"+str(byYearBiasDict4[kkk][1][1])\
                             +"\t"+str(byYearBiasDict4[kkk][2][0])+"\t"+str(byYearBiasDict4[kkk][2][1])+"\n")


        with open(outputFile4,"w") as output:
            output.write(cladeName+"_meanYBP"+"\t"+cladeName+"_dndo"+"\n")
            for kkkk in sorted(dndoResultDict2.keys()):
                output.write(str(kkkk)+"\t"+str(dndoResultDict2[kkkk])+"\n")

print("finished!!!!")





































