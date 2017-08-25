#coding=utf-8
from xml.etree import ElementTree as ET
import linecache
import math


def value_I(P):

    if P ==1.000000:
        logPai = math.log(0.990000,2)
    else:
        logPai = math.log(float(P),2)
    value_I = - (float(P) * logPai)
    return value_I


tree = ET.parse('Nuclear_receptor200_meme.xml')
root = tree.getroot()


protein_keys = {}
i = 0
for sequence in root.iter('sequence'):
    name = sequence.get('name')
    protein_keys[i] = name
    i += 1


letters_keys = {}
k = 0
for alphabet in root.iter('alphabet'):
    for letter in alphabet.iter('letter'):
        letterId = letter.get('id')
        letters_keys[letterId] = k
        k += 1


# initialize score_matrix, all 0
score_matrix_all = []
for proteinId in range (0,len(protein_keys)):
    score_matrix_oneline = []
    for motif_n in range (200):
        score_matrix_oneline.append(0)
    score_matrix_all.append(score_matrix_oneline)

# initialize sitesn_matrix, all 0
sitesn_matrix = []
for proteinId in range (0, len(protein_keys)):
    score_matrix_oneline = []
    for motif_n in range (200):
        score_matrix_oneline.append(0)
    sitesn_matrix.append(score_matrix_oneline)


motif_number = 0
for motif in root.iter('motif'):
    motif_number += 1
    motif_id = int(motif.get('id').split('_')[1]) #motif_id starts from 1
    for probability in motif.iter('probabilities'):
        alphabet_array_score_all = []
        for matrix in probability.iter('alphabet_matrix'):
            for array in matrix.iter('alphabet_array'):
                alphabet_array_score_oneline = []
                for value in array.iter('value'):
                    value_text = value.text
                    alphabet_array_score_oneline.append(value_text)
                alphabet_array_score_all.append(alphabet_array_score_oneline)

    for sites in motif.iter('contributing_sites'):
        for site in sites.iter('contributing_site'):
            sequenceId = site.get('sequence_id')
            sequenceIdKey = int(sequenceId.split('_')[1])
            letter_ref_counter = 0
            site_I_value = 0
            for letter_ref in site.iter('letter_ref'):
                letter_id = letter_ref.get('letter_id')
                site_I_value += value_I(float(alphabet_array_score_all[letter_ref_counter][letters_keys[letter_id]]))
                letter_ref_counter += 1
            score_matrix_all[sequenceIdKey][motif_id - 1] += site_I_value / letter_ref_counter
            sitesn_matrix[sequenceIdKey][motif_id -1] += 1
print sitesn_matrix
for i in range(len(protein_keys)):
    for j in range(motif_number):
        if sitesn_matrix[i][j] == 0:
            score_matrix_all[i][j] = 0
        else:
            score_matrix_all[i][j] = score_matrix_all[i][j] / sitesn_matrix[i][j]

with open('Nuclear_receptor200_is.txt','a') as out:
    for n in range (len(protein_keys)):
        out.write(protein_keys[n] +'\t')
        for m in range (motif_number): # the number of motifs
            out.write(str(score_matrix_all[n][m]) + '\t')
        out.write('\n')

out.close()















