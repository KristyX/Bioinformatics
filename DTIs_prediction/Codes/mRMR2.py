import numpy as np
with open ('.../output_files/mRMR/features/GPCR_features.txt', 'r') as fp_features:
    data_mrmr = fp_features.readlines()
with open ('.../output_files/mRMR//csv/mRMR_GPCR_is.csv', 'r') as fp_matrix:
    data_matrix = fp_matrix.readlines()
mrmr_index = []
mrmr_feature = []
for i in range(len(data_mrmr)):
    if data_mrmr[i].split()[3] != '-0.000' and data_mrmr[i].split()[3] != '-0.001':
        mrmr_index.append(data_mrmr[i].split()[1])
    else:
        continue
print len(mrmr_index)
print mrmr_index
mrmr_matrix = []
for j in range(len(data_matrix)):
    mrmr_feature.append([])
for p in range(len(data_matrix)):
    mrmr_feature[p].append(data_matrix[p].split(',')[0])
    for m in range(95):  #need to change
        mrmr_feature[p].append(data_matrix[p].split(',')[m+1])
    for q in range(len(mrmr_index)):
        mrmr_feature[p].append(data_matrix[p].split(',')[int(mrmr_index[q])])
print mrmr_feature
print len(mrmr_feature[0])

'''
with open('mRMR2_GPCR_is.csv', 'w') as out:
    for i in range(len(mrmr_feature)):
        for j in range(len(mrmr_feature[0])):
            out.write(mrmr_feature[i][j] + ',')
        out.write('\n')
'''

#transfer into .arff format
with open('mRMR2_GPCR_is.arff', 'w') as out:
    out.write('@RELATION    mRMR2_GPCR_is' + '\n')
    out.write('@ATTRIBUTE class {+1,-1}' + '\n')
    for i in range(len(mrmr_feature[0])-1):
        j = i + 1
        out.write('@ATTRIBUTE   ' + mrmr_feature[0][j] +  '      NUMERIC' + '\n')
    out.write('@DATA' + '\n')
    for a in range(len(mrmr_feature)-1):
        for b in range(len(mrmr_feature[0])):
            out.write(mrmr_feature[a+1][b] + '\t')
        out.write('\n')

