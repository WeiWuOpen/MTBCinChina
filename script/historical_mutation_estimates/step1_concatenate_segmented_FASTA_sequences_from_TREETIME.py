

cladeNameList=["4.2"]


for cladeName in cladeNameList:

    inputfile =r"...\script_and_testData\testData\treetimeResult_refCRMA\\"+cladeName+"_ancestral_results\\ancestral_sequences.fasta"

    outputfile=r"...\script_and_testData\testData\treetimeResult_refCRMA\\"+cladeName+"_ancestral_results\\ancestral_sequences_r.fasta"

    linelist=[]
    with open(inputfile,"r") as input:
        for line in input:
            if line.strip() !="":
                linex=line.strip()
                linelist.append(linex)
    n=0
    with open(outputfile,"w") as output:
        for i in linelist:
            if i[0] == ">":
                n+=1
                if n ==1:
                    output.write(i + "\n")
                else:
                    output.write("\n"+ i + "\n")
            else:
                output.write(i)


print("finished!!!!!!")