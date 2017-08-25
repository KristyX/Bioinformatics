from bs4 import BeautifulSoup
import urllib,xml.dom.minidom,re,string

with open('.../Data/Enzyme.txt','r') as fp:
    data=fp.readlines()
    #protein_name is hsa:column
    protein_name=[]
    str='has:10'
    for line in data:
        if str != line.split()[0]:
            protein_name.append(line.split()[0])
            str=line.split()[0]
fp.close()
print len(protein_name)

with open('lon_channel.fasta','a') as out:

    for m in range(0, len(protein_name)):
        #set the web address
        web='http://www.genome.jp/dbget-bin/www_bget?-f+-n+a+'+protein_name[m]

        #change the website to 'sequencepart.html',rewrite
        urllib.urlretrieve(web,'sequencepart.html')

        try:
            #read the elements in sequencepart.html;
            soup=BeautifulSoup(open('sequencepart.html'))

        except Exception, e:

            print 'cannot find the html for the protein: ' , protein_name[m]
            continue

        #print 'Approved: ',protein_name[m]

        seq=soup.pre.text
        pn=seq.split('\n')[1].split(' ', 1)[0]
        ps=''.join(seq.split('\n')[2:])
        protein_seq=pn+'\n'+ps
        print pn
        out.write(protein_seq+'\n')









