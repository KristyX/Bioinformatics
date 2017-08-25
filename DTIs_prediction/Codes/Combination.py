import numpy as np
import math
import random
import collections
with open('.../output_files/Nuclear_receptor/Nuclear_receptor.txt','r') as fp_dataset:   #Enzyme(c)
    interact_data = fp_dataset.readlines()
fp_dataset.close()
with open ('.../output_files/Nuclear_receptor/4860Nuclear_receptor_csfMatrix.txt','r') as fp_drug:
    drug_data = fp_drug.readlines()
fp_drug.close()
with open ('.../output_files/Nuclear_receptor/Nuclear_receptor100_is.txt','r') as fp_protein:
    protein_data = fp_protein.readlines()
fp_protein.close()
interact_Matrix = []
drug_Matrix = {}
protein_Matrix = {}
protein_names = []
drug_names = []
for m in range(len(interact_data)):
    interact_Matrix.append(interact_data[m].split())
for i in range(len(drug_data)):
    drug_Matrix[(drug_data[i].split('\t'))[0]] = (drug_data[i].split('\t'))[1:4861]
    drug_names.append((drug_data[i].split('\t'))[0])
for j in range(len(protein_data)):
    protein_Matrix[(protein_data[j].split('\t'))[0]] = (protein_data[j].split('\t'))[1: 101] #need to change (+1)
    protein_names.append((protein_data[j].split('\t'))[0])
drugTargets_dic = {}
drug_previous = []
proteinDrugs_dic = {}
protein_previous = ''
for d in range(len(interact_Matrix)):
    if interact_Matrix[d][1] not in drug_previous:
        drugTargets_dic[interact_Matrix[d][1]] = [interact_Matrix[d][0]]
        drug_previous.append(interact_Matrix[d][1])
    elif interact_Matrix[d][1] in drug_previous:
        drugTargets_dic[interact_Matrix[d][1]].append(interact_Matrix[d][0])
    if interact_Matrix[d][0] != protein_previous:
        proteinDrugs_dic[interact_Matrix[d][0]] = [interact_Matrix[d][1]]
        protein_previous = interact_Matrix[d][0]
    elif interact_Matrix[d][0] == protein_previous:
        proteinDrugs_dic[interact_Matrix[d][0]].append(interact_Matrix[d][1])
#dic could be drugTargets_dic or proteinDrugs_dic; matrix could be protein_Matrix or drug_Matrix; names could be protein_names or drug_names
def RNselection(dic, matrix, names):
    features = []
    meanx = []
    varx = []
    RN_samples = []
    for n in range(len(matrix.values()[0])):
        features.append([float(matrix.values()[0][n])])
    for m in range(len(matrix)-1):
        for n in range(len(matrix.values()[0])):
            #print m, n
            features[n].append(float(matrix.values()[m+1][n]))
    w_up = []
    for i in range(len(features)):
        if np.var(features[i]) == 0:
            w_up.append(0.0)
        else:
            w_up.append(math.sqrt(np.mean(features[i])) / float(np.var(features[i])))
        meanx.append(np.mean(features[i])) #mean of first second.. feature
        varx.append(np.var(features[i]))  #var of first second... feature
    w_down = sum(w_up)
    for k in dic: #k is key of dic, dic[k] is value
        positive_s = dic[k]
        negative_s = list(names)
        sum_positive_zai = 0
        prob_dic = {}
        #print negative_s
        for l in range(len(dic[k])):
            negative_s.remove(dic[k][l])
        for name in positive_s:
            for p in range(len(matrix[name])):
                if varx[p] != 0:
                    sum_positive_zai += abs((meanx[p] * (float(matrix[name][p]) - meanx[p])) / float(varx[p] * w_down))
                elif varx[p] == 0:
                    continue
        mean_positive_zai = sum_positive_zai / len(positive_s)
        for name in negative_s:
            negative_zai = 0
            for q in range(len(matrix[name])):
                if varx[q] != 0:
                    negative_zai += abs((meanx[q] * (float(matrix[name][q]) - meanx[q])) / float(varx[q] * w_down))
                elif varx[q] == 0:
                    continue
            prob = abs(negative_zai - mean_positive_zai)
            prob_dic[name] = prob
        prob_value = sorted(prob_dic.values(), reverse=True)
        prob_selected = prob_value[: len(positive_s)]

        for s in range(len(prob_selected)):
            one_key = random.choice([key for key, v in prob_dic.iteritems() if v == prob_selected[s]])
            if 'hsa:' in k:
                RN_samples.append([k, one_key])
            elif 'D' in k:
                RN_samples.append([one_key, k])
            else:
                print k + ' error'

    return RN_samples
non_interact_Matrix = []
comb = RNselection(drugTargets_dic, protein_Matrix, protein_names) + RNselection(proteinDrugs_dic, drug_Matrix, drug_names)
comb_set = [(a,b) for a,b in comb]
RN = [item for item, count in collections.Counter(comb_set).items() if count > 1]  # combinations occurring more than once
non_interact_Matrix = RN
RN_comb = list(set(comb_set) - set(RN))
np.random.shuffle(RN_comb)
non_interact_Matrix += list(RN_comb)[:(len(interact_Matrix) - len(RN))] #left protein, right drug   (len(interact_Matrix) - len(RN))

#Combine interaction features
interaction_features = []
for i in range(len(interact_Matrix)):
    protein_name = interact_Matrix[i][0]
    drug_name = interact_Matrix[i][1]
    interaction_features.append(protein_Matrix[protein_name] + drug_Matrix[drug_name])
non_interaction_features = []
for j in range(len(non_interact_Matrix)):
    protein_name = non_interact_Matrix[j][0]

    drug_name = non_interact_Matrix[j][1]

    non_interaction_features.append(protein_Matrix[protein_name] + drug_Matrix[drug_name])
validate_matrix = interact_Matrix + non_interact_Matrix
vm_set = [(a,b) for a,b in validate_matrix]
test = [item for item, count in collections.Counter(vm_set).items() if count > 1]
print len(protein_Matrix.values()[0])
print len(drug_Matrix.values()[0])
print len(test)
print interaction_features
print non_interaction_features
print len(non_interact_Matrix)

# write original arff
with open('Nuclear_receptor_is.arff', 'w') as out:
    out.write('@RELATION    Nuclear_receptor_is.arff' + '\n')
    out.write('@ATTRIBUTE class {1,2}' + '\n')
    for i in range(100):  # need to change
        j = i + 1
        out.write('@ATTRIBUTE   motif' + str(j) + '      NUMERIC' + '\n')
    for m in range(4860):  #need to change
        out.write('@ATTRIBUTE   csf' + str(m) + '    NUMERIC' + '\n')
    out.write('@DATA' + '\n')
    for m in range(len(interaction_features)):
        out.write('1' + '\t')
        for n in range(len(interaction_features[0])):
            out.write(interaction_features[m][n] + '\t')
        out.write('\n')
    for m in range(len(non_interaction_features)):
        out.write('2' + '\t')
        for n in range(len(non_interaction_features[0])):
            out.write(non_interaction_features[m][n] + '\t')
        out.write('\n')


# write as .csv format
'''
with open('mRMR_Enzyme_is.csv','w') as out:
    out.write('label,')
    for i in range(len(protein_Matrix.values()[0])):
        out.write('motif' + str(i+1) + ',')
    for j in range (len(drug_Matrix.values()[0])-1):
        out.write('csf' + str(j) + ',')
    out.write('csf' + str(len(drug_Matrix.values()[0])-1) +'\n')
    for m in range(len(interaction_features)):
        out.write('+1' + ',')
        for n in range(len(interaction_features[0])-1):
            out.write(interaction_features[m][n] + ',')
        out.write(interaction_features[m][len(interaction_features[0])-1] + '\n')
    for m in range(len(non_interaction_features)):
        out.write('-1' + ',')
        for n in range(len(non_interaction_features[0])-1):
            out.write(non_interaction_features[m][n] + ',')
        out.write(non_interaction_features[m][len(non_interaction_features[0])-1] + '\n')
'''

# all non-interacted data
'''
all_comb = []
for protein in protein_Matrix.keys():
    for drug in drug_Matrix .keys():
        all_comb.append([protein, drug])
for comb in interact_Matrix:
    if comb in all_comb:
        all_comb.remove(comb)
unknown_Matrix = all_comb
# Ramdonly select non-intraction data
np.random.shuffle(all_comb)
non_interact_Matrix = all_comb[:len(interact_Matrix)]

unknown_features = []
for i in range(len(unknown_Matrix)):
    protein_name2 = unknown_Matrix[i][0]
    drug_name2 = unknown_Matrix[i][1]
    unknown_features.append(protein_Matrix[protein_name2] + drug_Matrix[drug_name2])


#weighted combination
weight = len(non_interact_Matrix) / len(interact_Matrix)
with open ('1000000000wRN_GPCR_is.arff','w') as out:
    out.write('@RELATION    1000000000wRN_GPCR_is' + '\n') #need to change??
    out.write('@ATTRIBUTE class {1,2}' + '\n')
    for i in range(95): #need to change
        j = i + 1
        out.write('@ATTRIBUTE   motif' + str(j) + '      NUMERIC' + '\n')
    for m in range(881): #need to change
        out.write('@ATTRIBUTE   csf' + str(m) + '    NUMERIC' + '\n')
    out.write('@DATA' + '\n')
    for m in range(len(interaction_features)):
        out.write('1' + '\t')
        for n in range(len(interaction_features[0])):
            out.write(interaction_features[m][n] + '\t')
        out.write('{' + str(weight * 1000000000) + '}' + '\n')
    for m in range(len(non_interaction_features)):
        out.write('2' + '\t')
        for n in range(len(non_interaction_features[0])):
            out.write(non_interaction_features[m][n] + '\t')
        out.write('{1}' + '\n')
'''
'''
with open('0.6GPCR_sws.arff','w') as out:
    out.write('@RELATION    Nuclear_receptor_is' + '\n')
    out.write('@ATTRIBUTE class {1,2}' + '\n')
    for i in range(95):  # need to change
        j = i + 1
        out.write('@ATTRIBUTE   motif' + str(j) + '      NUMERIC' + '\n')
    for m in range(881):
        out.write('@ATTRIBUTE   csf' + str(m) + '    NUMERIC' + '\n')
    out.write('@DATA' + '\n')
    for m in range(len(interaction_features)):
        out.write('1' + '\t')
        for n in range(len(interaction_features[0])):
            out.write(interaction_features[m][n] + '\t')
        out.write('\n')
    for m in range(len(non_interaction_features)):
        out.write('2' + '\t')
        for n in range(len(non_interaction_features[0])):
            out.write(non_interaction_features[m][n] + '\t')
        out.write('\n')
'''
