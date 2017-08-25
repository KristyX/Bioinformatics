#coding=utf-8
from xml.etree import ElementTree as ET

protein_keys = {}
keys_sequence = {}
with open('Nuclear_receptor.txt','r') as fp:

    data = fp.readlines()
fp.close()
m = 0
n = 0
for i in range(len(data)):
    if i % 2 == 0:
        protein_keys[m] = data[i].replace('>', '').replace('\n', '')
        m += 1
    else:
        keys_sequence[n] = str(data[i].replace('\n', ''))
        n += 1
tree = ET.parse('Nuclear_receptor200_meme.xml')#change the MEME output file name to combinedmeme.xml
root = tree.getroot()
letters_keys = {}
k = 0
for alphabet in root.iter('alphabet'):
    for letter in alphabet.iter('letter'):
        letterId = letter.get('id')
        letters_keys[letterId] = k
        k += 1

'''

def add_char(str):
    output_sequence = ['']
    i=0
    while i < len(str):
        if str[i] == '[':
            i=i+1
            temp_sequence = []
            while str[i] != ']':
                for output in output_sequence:
                    temp_sequence.append(output + str[i])
                i=i+1
            output_sequence = temp_sequence
        elif str[i] == ']':
            i=i+1

        else:
            for j in range(len(output_sequence)):
                output_sequence[j]=output_sequence[j]+str[i]
            i=i+1
    return output_sequence
'''
score_matrix_all = []
for proteinId in range(0,len(protein_keys)):
    score_matrix_oneline = []
    for motif_n in range(200):
        score_matrix_oneline.append(0)
    score_matrix_all.append(score_matrix_oneline)
motif_number = 0
for motif in root.iter('motif'):
    motif_number += 1
    motif_id = int(motif.get('id').split('_')[1]) #motif_id starts from 1
    '''
    sites = []
    for regular_expression in motif.iter('regular_expression'):
        expression = str(regular_expression.text)
        site = list(expression)
        for item in site:
            if item == '\n':
                site.remove(item)
            else:
                continue

        sites = add_char(site)
        #print sites
    '''
    for probability in motif.iter('probabilities'):
        alphabet_array_score_all = []
        for matrix in probability.iter('alphabet_matrix'):
            for array in matrix.iter('alphabet_array'):
                alphabet_array_score_oneline = []
                for value in array.iter('value'):
                    value_text = value.text
                    #print value_text
                    alphabet_array_score_oneline.append(value_text)
                alphabet_array_score_all.append(alphabet_array_score_oneline)


    for x in range (0, len(keys_sequence)):
        sequence = keys_sequence[x]
        counter = 0
        score = 0.0
        for slide_i in range (0, len(sequence) - len(alphabet_array_score_all) + 1):
            window_slide_sub_sequence = sequence[slide_i : slide_i + len(alphabet_array_score_all)]
            temp_score = 0.0
            for ii in range(0, len(alphabet_array_score_all)):     # current site's score in current sequence
                letterName = str(window_slide_sub_sequence)[ii]
                temp_score += float(alphabet_array_score_all[ii][letters_keys[letterName]])
            temp_score = temp_score / len(alphabet_array_score_all)
            if temp_score >= 0.4:  # lumda need to be changed here!!!
                counter += 1
                score = score + temp_score
            else:
                continue
        if counter == 0:
            score = 0
        else:
            score = score / counter
        score_matrix_all[x][motif_id - 1] += score

with open('0.4Nuclear_receptor200_sws.txt','a') as out:
    for n in range (len(protein_keys)):
        out.write(protein_keys[n] +'\t')
        for m in range (motif_number): # the number of motifs
            out.write(str(score_matrix_all[n][m]) + '\t')
        out.write('\n')
out.close()







