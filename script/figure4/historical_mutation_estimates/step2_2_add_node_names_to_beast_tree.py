

inputFolder=r"...\script_and_testData\testData\beastTree\\"
outputFolder=r"...\script_and_testData\testData\output\historical_mutation_estimates\\"



import os
from tqdm import tqdm

if not os.path.exists(outputFolder):
    os.mkdir(outputFolder)

inputFileList=os.listdir(inputFolder)


for f in tqdm(inputFileList):

    with open(inputFolder+f,"r") as input:
        for ll in input:
                tree=ll.strip()

    treeList=tree.split(")")


    newTree=""
    nodeNum=0
    for ii in range(0,len(treeList)):
        if ii != len(treeList)-1:
            newTree += treeList[ii]+")"+"NodeX"+str(nodeNum)
            nodeNum+=1
        else:
            newTree += treeList[ii]


    with open(outputFolder+f.replace("newick","")+"_nodename.newick","w") as output:
        output.write(newTree.replace("'","")) # 去除figtree给样本名添加的' '

print("finished!!!!")