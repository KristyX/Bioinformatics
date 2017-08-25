import networkx as nx
import operator
import numpy as np
import matplotlib.pyplot as plt



with open('S.cerevisiae.txt') as fp:
    data_all = fp.readlines()
fp.close()
miu = 0.125
module_size = 10
g = nx.Graph()
for i in range(len(data_all)):
    g.add_edge(data_all[i].split()[0], data_all[i].split()[1], weight=data_all[i].split()[2])
pk = []
for item in g.nodes():
    pk.append([item])
length = len(g.nodes()) + 1
adjc_matrix = [[0 for i in range(length)] for j in range(length)]
degree_kv = {}
# store weight between each two nodes
#w = nx.get_edge_attributes(g, 'weight')
for i in range(len(g.nodes())):
    adjc_matrix[0][i+1] = g.nodes()[i]
    adjc_matrix[i+1][0] = g.nodes()[i]
    # loop neighbors of g.node()[i], and mark their weight in adjc_matrix
    degree = 0
    for nb in g.neighbors(g.nodes()[i]):
        adjc_matrix[i+1][g.nodes().index(nb)+1] = float(g.get_edge_data(g.nodes()[i], nb)['weight'])
        adjc_matrix[g.nodes().index(nb)+1][i+1] = float(g.get_edge_data(g.nodes()[i], nb)['weight'])
        degree = degree + adjc_matrix[i+1][g.nodes().index(nb)+1]
    degree_kv[adjc_matrix[0][i+1]] = degree
#print adjc_matrix
#pk is non-overlapping subnetwork

#print pk
find_weight = []

kv_decreased = sorted(degree_kv.items(), key=operator.itemgetter(1), reverse=True)

def R_values(u, v, g):
    weight_up = float(g.get_edge_data(u, v)['weight'])
    #cns is a set of common neighbors
    cns = []
    for i in range(len(g.neighbors(u))):
        for j in range(len(g.neighbors(v))):
            if g.neighbors(u)[i] == g.neighbors(v)[j]:
                cns.append(g.neighbors(v)[j])

    if cns != []:
        for m in range(len(cns)):
            weight_up = weight_up + float(g.get_edge_data(u, cns[m])['weight']) + float(g.get_edge_data(v, cns[m])['weight'])
    for n in range(len(kv_decreased)):
        if u == kv_decreased[n][0]:
            weight_down = kv_decreased[n][1]
    r = weight_up / weight_down
    return r

clusters_matrix = []
for i in range(len(kv_decreased)):
    clusters_matrix.append([kv_decreased[i][0]])

assist_array = []
for j in range(len(kv_decreased)):
    assist_array.append(kv_decreased[j][0])



for p in range(len(kv_decreased)):
    vertex_v = kv_decreased[p][0]
    for vertex_u in g.neighbors(vertex_v):
        Ruv = R_values(vertex_u, vertex_v, g)
        Rvu = R_values(vertex_v, vertex_u, g)
        if Ruv >= Rvu and Ruv > 0.5*miu:
            if vertex_u in assist_array:
                clusters_matrix[p].append(vertex_u)
                assist_array.remove(vertex_u)
            else:
                continue
#clusters_matrix = [['a','d','e'],['b','c','f'],['c','h'],['d','g','k']]
def numpy_unique_ordered(seq):
    """Remove duplicate from a list while keeping order with Numpy.unique
    Required:
        seq: A list containing all items
    Returns:
        A list with only unique values. Only the first occurence is kept.
    """
    array_unique = np.unique(seq, return_index=True)
    dstack = np.dstack(array_unique)
    dstack.dtype = np.dtype([('v', dstack.dtype), ('i', dstack.dtype)])
    dstack.sort(order='i', axis=1)
    return dstack.flatten()['v'].tolist()
clusters_matrix.reverse()
length_cm = len(clusters_matrix)
j = 0
for i in range(length_cm):
    length_cc = len(clusters_matrix)
    for clusters in clusters_matrix[j+1: ]:
        if clusters_matrix[j][0] in clusters:
            clusters_matrix[clusters_matrix.index(clusters)] = numpy_unique_ordered((clusters_matrix[clusters_matrix.index(clusters)] + clusters_matrix[j]))
            clusters_matrix.remove(clusters_matrix[j])
    if length_cc != len(clusters_matrix):
        j = j
    else:
        j = j+1

modules = []
size_array = []
sum_size = 0
for module in clusters_matrix:
    if len(module) >= module_size:
        modules.append(module)
        print module
        size_array.append(len(module))
        sum_size = sum_size + len(module)

min_size = min(size_array)
max_size = max(size_array)
mean_size = sum_size / len(modules)
print "miu = %f, size = %d" % (miu, module_size)
print 'Amount of clusters: %d' % len(clusters_matrix)
print 'Amount of functional modules: %d' % len(modules)
print 'Minimum size is %d' % min_size
print 'Maximum size is %d' % max_size
print 'Average size is %d' % mean_size

'''
sys.setrecursionlimit(10000)
def getList(arr, flag):
    flag = 1
    for k in xrange(len(arr)):
        header = arr[k][0]
        if k == len(arr)-1:
            break;
        for i in xrange(k+1,len(arr),1):
            if header in arr[i]:
                flag = 0
                arr[i] = arr[i] + arr[k][1:]
                arr.remove(arr[k])
                return getList(arr,flag)
            else:
                continue
    if flag == 1:
        return arr
    k = getList(arr,0)
    print k
print getList(clusters_matrix, 0)
'''

#plot
#subgrf are clusters of g.
pos = nx.spring_layout(g)
subgrf = []
for m in range(len(modules)):
    subgrf.append([])
colors = ['b','g','r','c', 'm', 'y']
for i in range(len(modules)):
    subgrf[i] = g.subgraph(modules[i])
    nx.draw(subgrf[i], pos=pos, node_size=10, node_color=colors[i % 6], edge_color=colors[i % 6])
#plt.savefig('/Users/Kristy/Documents/2016Fall/539Database/project/graphs/SSS.png')
plt.show()
