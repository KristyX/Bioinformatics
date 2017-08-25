import re

with open('.../output_files/protein_features/meme_output/Enzyme.txt','r') as fp1:
    data1=fp1.readlines()
    #protein_name is hsa:column
    drug_name=[]
    for line in data1:
        if line.split()[1] in drug_name:
            continue
        else:
            drug_name.append(line.split()[1])
fp1.close()

with open ('.../output_files/drug_features/R_fp/4860Enzyme_fp.txt', 'r') as fp2:
    data2 = fp2.readlines()
fp2.close()
i = 5
csf_position = []
while i <= len(data2):
    csf_position.append([int(s) for s in data2[i].split() if s.isdigit()])
    i += 6
csf_matrix =[]
for row in range(len(drug_name)):
    csf_temp = []
    for column in range(4860): #PubChem dataset: 881
        csf_temp.append(0)
    csf_matrix.append(csf_temp)
for l in range(len(drug_name)):
    for c in range(len(csf_position[l])):
        n = int(csf_position[l][c])
        csf_matrix[l][n] = 1

with open('4860Enzyme_csfMatrix.txt', 'w') as out:
    for n in range (len(drug_name)):
        out.write(drug_name[n] +'\t')
        for m in range (4860): # the number of csf  PubChem
            out.write(str(csf_matrix[n][m]) + '\t')
        out.write('\n')
out.close()





