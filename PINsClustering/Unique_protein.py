with open('protein_name.txt', 'w+') as out:
    with open('S.cerevisiae.txt') as fp:
        all_data = fp.readlines()
        # unique protein
        proteins = []
        for i in range(len(all_data)):
            proteins.append(all_data[i].split()[0])
            proteins.append(all_data[i].split()[1])
        protein = list(set(proteins))
        '''
        for i in range(len(protein)):
            out.write (protein[i] + "\n")
        '''
        adjc_matrix = []
        for i in range(len(protein)):
            adjc_matrix[0][i+1] = protein[i]
            adjc_matrix[i+1][0] = protein[i]


