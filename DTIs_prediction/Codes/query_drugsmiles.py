from bs4 import BeautifulSoup
import urllib,xml.dom.minidom,re,string
import requests
import lxml.html as LH

with open('Enzyme.txt','r') as fp:
    data=fp.readlines()
    #protein_name is hsa:column
    drug_name=[]
    for line in data:
        if line.split()[1] in drug_name:
            continue
        else:
            drug_name.append(line.split()[1])

fp.close()

print len(drug_name)

with open('Enzyme_smiles.txt','a') as out:
    smiles = []
    for dname in drug_name:
        web1='http://www.kegg.jp/dbget-bin/www_bget?dr:' + dname
        urllib.urlretrieve(web1,'sequencepart1.html')
        try:
            #read the elements in sequencepart.html;
            soup1=BeautifulSoup(open('sequencepart1.html'))
        except Exception, e:
            print 'cannot find the html for the drug according to drug_name: ', dname
            continue

        several_p = soup1.findAll(attrs={'style': 'margin-left:5em'})
        for i in range(len(several_p)):
            href_text=several_p[i].find('a').get('href')
            if 'pubchem' in href_text:
                pubchem_id= re.findall(r'\b\d+\b', href_text)[0]
                web2 = 'https://pubchem.ncbi.nlm.nih.gov/rest/rdf/descriptor/CID' + pubchem_id +'_Canonical_SMILES.html'
                #https://pubchem.ncbi.nlm.nih.gov/rest/rdf/descriptor/CID5831_Canonical_SMILES.html
                urllib.urlretrieve(web2,'sequencepart2.html')
                try:
                    soup2 = BeautifulSoup(open('sequencepart2.html'))
                except Exception, e:
                    print 'cannot find the drug according to pubchem_id: ', id
                    continue
                smiles.append(str(soup2.find('span',{'class':'value'}).get_text()))
            else:
                continue





        #smiles.append(smiles_content.find('div',{'class': 'section-content-item'}).text)
    out.write(str(smiles))

