

cladeNameList=["4.2cnONLYchina","4.2bNOchina"]

import math

from tqdm import tqdm


for cladeName in cladeNameList:
    print(cladeName + " is working #######")

    inputlistFile = r"...\script_and_testData\testData\sampleList\\" + cladeName + ".txt"
    inputrefseqFile = r"...\script_and_testData\testData\H37Rv_r.fasta"
    inputrefseqSNPFile_new = r"...\script_and_testData\testData\aMTBC_noCanetti_withprd.fasta"
    inputrefseqPOSFile_new = r"...\script_and_testData\testData\all_exceptAnimal_withCanetti_refH37Rv_withprd.pos.txt"
    inputAnnotationFile = r"...\script_and_testData\testData\H37Rv.annotation_all.txt"
    inputvcfFolder = r"...\script_and_testData\testData\vcf\\"
    inputuncoverFolder = r"...\script_and_testData\testData\uncoverBed\\"

    outputFile=r"...\script_and_testData\testData\output\information_entropy\\"+cladeName+"_aMTBC_DFandKLvalues.txt"

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


    strainlist = []
    with open(inputlistFile, "r") as input:
        for l in input:
            strainlist.append(l.strip())
    sampleNum=len(strainlist)

    geneDict = {}
    nongeneDict ={}
    with open(inputAnnotationFile, "r") as input:
        for l in input:
            lx = l.strip().split("\t")
            if lx[2] == "CDS":
                if lx[0][-1] == "c":
                    geneDict[lx[0]] = lx[3] + "_" + lx[4] + "_" +  lx[2]
                else:
                    geneDict[lx[0]] = lx[3] + "_" + lx[4] + "_" +  lx[2]
            else:
                if lx[2] != "repeat_region" and lx[4] != "5pUTR":
                    nongeneDict[lx[0]] =  lx[3] + "_" + lx[4] + "_" + lx[2]

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


    # # 交叉熵
    # def H(pA,pT,pC,pG,qA,qT,qC,qG): # 前者为参考
    #     # import math
    #     e=math.e
    #     result= - (pA*math.log(qA,e) +pT*math.log(qT,e) +pC*math.log(qC,e) +pG*math.log(qG,e))
    #     return result

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

        # if pA==0:
        #     result = (-(pT * math.log(qT, e) + pC * math.log(qC, e) + pG * math.log(qG, e))) - \
        #              (- (pT * math.log(pT, e) + pC * math.log(pC, e) + pG * math.log(pG, e)))
        # elif pT==0:
        #     result = (- (pA * math.log(qA, e)  + pC * math.log(qC, e) + pG * math.log(qG, e))) - \
        #              (- (pA * math.log(pA, e)  + pC * math.log(pC, e) + pG * math.log(pG, e)))
        # elif pC==0:
        #     result = (- (pA * math.log(qA, e) + pT * math.log(qT, e)  + pG * math.log(qG, e))) - \
        #              (- (pA * math.log(pA, e) + pT * math.log(pT, e)  + pG * math.log(pG, e)))
        # elif pG==0:
        #     result = (- (pA * math.log(qA, e) + pT * math.log(qT, e) + pC * math.log(qC, e))) - \
        #              (- (pA * math.log(pA, e) + pT * math.log(pT, e) + pC * math.log(pC, e)))
        # else:
        #     result= (- (pA*math.log(qA,e) +pT*math.log(qT,e) +pC*math.log(qC,e) +pG*math.log(qG,e))) - \
        #         (- (pA * math.log(pA, e) + pT * math.log(pT, e) + pC * math.log(pC, e) + pG * math.log(pG, e)))
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
        # if pA==0:
        #     result = - (pT * math.log(qT, e) + pC * math.log(qC, e) + pG * math.log(qG, e))
        # elif pT==0:
        #     result = - (pA * math.log(qA, e)  + pC * math.log(qC, e) + pG * math.log(qG, e))
        # elif pC==0:
        #     result = - (pA * math.log(qA, e) + pT * math.log(qT, e)  + pG * math.log(qG, e))
        # elif pG==0:
        #     result = - (pA * math.log(qA, e) + pT * math.log(qT, e) + pC * math.log(qC, e))
        # else:
        #     result = - (pA * math.log(qA, e) + pT * math.log(qT, e) + pC * math.log(qC, e) + pG * math.log(qG, e))
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


    # 计算4种碱基出现的频率
    baseRateBYgeneDict={}
    baseRateBYnongeneDict={}
    for ss in tqdm(strainlist):
        snpPOSDict={}
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
                snpPOSDict[int(p.split("_")[0])] = [p.split("_")[1],p.split("_")[2]]

        for gg in geneDict.keys():
            if gg not in baseRateBYgeneDict.keys():
                baseRateBYgeneDict[gg]=[[0,0,0,0],[0,0,0,0]]
            numA = 0
            numT = 0
            numC = 0
            numG = 0
            refnumA = 0
            refnumT = 0
            refnumC = 0
            refnumG = 0
            gllimit= int(geneDict[gg].split("_")[0])
            grlimit = int(geneDict[gg].split("_")[1])
            geneLength=grlimit-gllimit+1

            for ii in range(gllimit,grlimit+1):
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
                    if snpPOSDict[ii][1] == "A":
                        numA += 1
                    elif snpPOSDict[ii][1] == "T":
                        numT += 1
                    elif snpPOSDict[ii][1] == "C":
                        numC += 1
                    elif snpPOSDict[ii][1] == "G":
                        numG += 1
                    else:
                        print(snpPOSDict[ii][1])

                    if snpPOSDict[ii][0] == "A":
                        refnumA += 1
                    elif snpPOSDict[ii][0] == "T":
                        refnumT += 1
                    elif snpPOSDict[ii][0] == "C":
                        refnumC += 1
                    elif snpPOSDict[ii][0] == "G":
                        refnumG += 1
                    else:
                        print(snpPOSDict[ii][0])


            baseRateBYgeneDict[gg][0][0] += (refnumA/geneLength) / sampleNum
            baseRateBYgeneDict[gg][0][1] += (refnumT / geneLength) / sampleNum
            baseRateBYgeneDict[gg][0][2] += (refnumC / geneLength) / sampleNum
            baseRateBYgeneDict[gg][0][3] += (refnumG / geneLength) / sampleNum
            baseRateBYgeneDict[gg][1][0] += (numA / geneLength) / sampleNum
            baseRateBYgeneDict[gg][1][1] += (numT / geneLength) / sampleNum
            baseRateBYgeneDict[gg][1][2] += (numC / geneLength) / sampleNum
            baseRateBYgeneDict[gg][1][3] += (numG / geneLength) / sampleNum

        for gg in nongeneDict.keys():
            if gg not in baseRateBYnongeneDict.keys():
                baseRateBYnongeneDict[gg]=[[0,0,0,0],[0,0,0,0]]
            numA = 0
            numT = 0
            numC = 0
            numG = 0
            refnumA = 0
            refnumT = 0
            refnumC = 0
            refnumG = 0
            gllimit= int(nongeneDict[gg].split("_")[0])
            grlimit = int(nongeneDict[gg].split("_")[1])
            nongeneLength=grlimit-gllimit+1

            for ii in range(gllimit,grlimit+1):
                if ii not in snpPOSDict.keys():
                    if allLongNewRefseqdict[ii] == "A":   # 注意：不考虑基因顺反的原因是，在反义链上的基因交换A和T以及G和C之后不影响信息系统的组成
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
                    if snpPOSDict[ii][1] == "A":
                        numA += 1
                    elif snpPOSDict[ii][1] == "T":
                        numT += 1
                    elif snpPOSDict[ii][1] == "C":
                        numC += 1
                    elif snpPOSDict[ii][1] == "G":
                        numG += 1
                    else:
                        print(snpPOSDict[ii][1])

                    if snpPOSDict[ii][0] == "A":
                        refnumA += 1
                    elif snpPOSDict[ii][0] == "T":
                        refnumT += 1
                    elif snpPOSDict[ii][0] == "C":
                        refnumC += 1
                    elif snpPOSDict[ii][0] == "G":
                        refnumG += 1
                    else:
                        print(snpPOSDict[ii][0])


            baseRateBYnongeneDict[gg][0][0] += (refnumA / nongeneLength) / sampleNum
            baseRateBYnongeneDict[gg][0][1] += (refnumT / nongeneLength) / sampleNum
            baseRateBYnongeneDict[gg][0][2] += (refnumC / nongeneLength) / sampleNum
            baseRateBYnongeneDict[gg][0][3] += (refnumG / nongeneLength) / sampleNum
            baseRateBYnongeneDict[gg][1][0] += (numA / nongeneLength) / sampleNum
            baseRateBYnongeneDict[gg][1][1] += (numT / nongeneLength) / sampleNum
            baseRateBYnongeneDict[gg][1][2] += (numC / nongeneLength) / sampleNum
            baseRateBYnongeneDict[gg][1][3] += (numG / nongeneLength) / sampleNum


    klResultDict_gene={}
    hResultDict_gene={}
    dfResultDict_gene={}
    for gene in baseRateBYgeneDict.keys():
        pA = baseRateBYgeneDict[gene][0][0]
        pT = baseRateBYgeneDict[gene][0][1]
        pC = baseRateBYgeneDict[gene][0][2]
        pG = baseRateBYgeneDict[gene][0][3]
        qA = baseRateBYgeneDict[gene][1][0]
        qT = baseRateBYgeneDict[gene][1][1]
        qC = baseRateBYgeneDict[gene][1][2]
        qG = baseRateBYgeneDict[gene][1][3]
        klresult=KL(pA,pT,pC,pG,qA,qT,qC,qG)
        klResultDict_gene[gene]=klresult
        hresult = H(pA, pT, pC, pG, qA, qT, qC, qG)
        hResultDict_gene[gene] = hresult
        dfresult=DF(pA, pT, pC, pG, qA, qT, qC, qG)
        dfResultDict_gene[gene] = dfresult
        del klresult
        del hresult
        del dfresult

    klResultDict_nongene={}
    hResultDict_nongene={}
    dfResultDict_nongene={}
    for nongene in baseRateBYnongeneDict.keys():
        pA = baseRateBYnongeneDict[nongene][0][0]
        pT = baseRateBYnongeneDict[nongene][0][1]
        pC = baseRateBYnongeneDict[nongene][0][2]
        pG = baseRateBYnongeneDict[nongene][0][3]
        qA = baseRateBYnongeneDict[nongene][1][0]
        qT = baseRateBYnongeneDict[nongene][1][1]
        qC = baseRateBYnongeneDict[nongene][1][2]
        qG = baseRateBYnongeneDict[nongene][1][3]

        klresult=KL(pA,pT,pC,pG,qA,qT,qC,qG)
        klResultDict_nongene[nongene]=klresult
        hresult = H(pA, pT, pC, pG, qA, qT, qC, qG)
        hResultDict_nongene[nongene] = hresult
        dfresult = DF(pA, pT, pC, pG, qA, qT, qC, qG)
        dfResultDict_nongene[nongene] = dfresult
        del klresult
        del hresult
        del dfresult

    with open(outputFile,"w") as output:
        output.write("Gene"+"\t"+"Type"+"\t"+"StartPos"+"\t"+"EndPos"+"\t"+cladeName+"_KLvalue"+"\t"+cladeName+"_Hvalue"+"\t"+cladeName+"_DFvalue"+"\n")
        for kr in klResultDict_gene.keys():
            output.write(kr +"\t"+ geneDict[kr].split("_")[2] +"\t"+ geneDict[kr].split("_")[0] +"\t"+ geneDict[kr].split("_")[1] +"\t"+ str(klResultDict_gene[kr]) +"\t"+ str(hResultDict_gene[kr]) +"\t"+ str(dfResultDict_gene[kr])+"\n")
        del kr
        for kr in klResultDict_nongene.keys():
            output.write(kr +"\t"+ nongeneDict[kr].split("_")[2] +"\t"+ nongeneDict[kr].split("_")[0] +"\t"+nongeneDict[kr].split("_")[1] +"\t"+str(klResultDict_nongene[kr])+"\t"+ str(hResultDict_nongene[kr]) +"\t"+ str(dfResultDict_nongene[kr])+"\n")


print("finished!!!!!!!")
