

cladeNameList=["4.2"]


for cladeName in cladeNameList:
    print(cladeName + " is working #######")

    inputFolder=r"F:\forStudy\finalResults\TB_MS\script_and_testData\testData\treetimeResult_refCRMA\\"
    inputfile=inputFolder+cladeName+"_ancestral_results\\annotated_tree.nexus"
    outputfile=inputFolder+cladeName+"_ancestral_results\\annotated_tree_del.nexus"

    titleList=[]
    with open(inputfile,"r") as input:
        for l in input:
            if l.strip().split()[0] !="Tree" and l.strip() !="End;" and l.strip() !="Begin Trees;":
                titleList.append(l)
            elif l.strip().split()[0] =="Tree":
                tree=l.strip().split()[1].replace("tree1=","")

    tree=tree.replace("]","[")
    treeList=tree.split("[")
    newTree=""
    for i in range(0,len(treeList),2):
        newTree+=treeList[i]
    # print(newTree)

    with open(outputfile,"w") as output:
        for ii in titleList:
            output.write(ii.replace("\n","")+"\n")
        output.write("End;"+"\n")
        output.write("Begin Trees;" + "\n")
        output.write(" "+"Tree tree1=")
        output.write(newTree+"\n")
        output.write("End;")


print("finished!!!!!")