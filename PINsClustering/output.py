import Matrix as mx
import requests
from bs4 import BeautifulSoup
import numpy as np
p_value = []
for i in range(len(mx.modules)):
    p_value.append([])
p1 = 0
p2 = 0
p3 = 0
p4 = 0
size1 = size2 = size3 = size4 = 0
count_error = 0

url = "http://www.yeastgenome.org/cgi-bin/GO/goTermFinder.pl/html"
for n in range(len(mx.modules)):
    str = ' '.join(mx.modules[n])
    payload = {'loci': str,
               'ontology': 'Function',
               'feature_type': ['ORF', 'blocked_reading_frame', 'ncRNA_gene', 'not in systematic sequence of S288C', 'pseudogene', 'rRNA_gene', 'snRNA_gene', 'snoRNA_gene', 'tRNA_gene', 'telomerase_RNA_gene'],
               'qualifier': ['transposable_element_gene', 'Dubious', 'Uncharacterized', 'Verified'],
               'manually curated': 'yes',
               'high-throughput': 'yes',
               'annotation_source': ['SGD', 'UniProt', 'InterPro', 'GO_Central', 'GOC', 'MGI', 'HGNC','CACAO'],
               'evidence_code': ['IC', 'IDA', 'IEP', 'IGC', 'IGI', 'IMP', 'IPI', 'ISA', 'ISM', 'ISO', 'ISS', 'NSA', 'ND', 'RCA', 'TAS'],
               'cutoff': '0.01',
               'fdr': 'fdr'}
    r = requests.post(url, payload)


    #with open("requests_results.html", "w") as f:
    #    f.write(r.text)

    soup = BeautifulSoup(r.text, 'lxml')

    try:

        p_value[n] = min([float(tag.text) for tag in soup.find_all(attrs={'class':'nowrap'})])
        if p_value[n] < 1e-15:
            p1= p1 + 1
            size1 = size1 + len(mx.modules[n])
        if (p_value[n] > 1e-15) and (p_value[n] < 1e-10):
            p2 = p2 + 1
            size2 = size2 + len(mx.modules[n])
        if (p_value[n] > 1e-10) and (p_value[n] < 1e-5):
            p3 = p3 +1
            size3 = size3 + len(mx.modules[n])
        if p_value[n] > 1e-5:
            p4 = p4 + 1
            size4 = size4 + len(mx.modules[n])

    except ValueError:
        print "Module %s is unknown." % str
        count_error = count_error +1

while ([] in p_value):
    p_value.remove([])
mean_p = np.mean(list(p_value))
if p1 != 0:
    ave_size1 = size1 / p1
else:
    ave_size1 = 0
if p2 != 0:
    ave_size2 = size2 / p2
else:
    ave_size2 = 0
if p3 != 0:
    ave_size3 = size3 / p3
else:
    ave_size3 = 0
if p4 != 0:
    ave_size4 = size4 / p4
else:
    ave_size4 = 0

print p_value
print "<E-15: N.modules: %d   Avg.size: %d" % (p1, ave_size1)
print "[E-15, E-10: N.modules: %d   Avg.size: %d" % (p2, ave_size2)
print "[E-10, E-5]: N.modules: %d   Avg.size: %d" % (p3, ave_size3)
print "[E-5, 1]: N.modules: %d   Avg.size: %d" % (p4, ave_size4)
print "Average of P values is %f" % mean_p
print "Amount of unknown modules: %d  (%.2f)" % (count_error, float(p1+p2+p3+p4)/len(mx.modules))




