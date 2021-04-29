import os
from Bio import SeqIO
from Bio import Entrez

#%%
##TRY 1
from Bio import Entrez
Entrez.email = "cyrusmallon@gmail.com" # Always tell NCBI who you are
handle = Entrez.esearch(db="nuccore", term="\"Escherichia coli O157\"[Organism] AND complete[properties] NOT plasmid", retmax="40", idtype="acc")
record1 = Entrez.read(handle)

#%%
##TRY 2
from Bio import Entrez
Entrez.email = "cyrusmallon@gmail.com" # Always tell NCBI who you are
handle2 = Entrez.esearch(db="nucleotide", term="Escherichia coli O157", retmax="40", idtype="acc")
record2 = Entrez.read(handle2)

#%%
print(record1["IdList"])
len(record1["IdList"])
print(record1["Count"])
#%%
##print length of dictionary
print(len(record1))
##print type of dictionary
print(type(record1))
##get keys for the dictionary
keys = record1.keys()
print(keys)
##get values for dictionary
values = record1.values()
print(values)
##get the list of key:value pairs
items = record1.items()
print(items)

#%%
##Are the search results of TRY 1 and TRY 2 the same?
print(record1["IdList"],record2["IdList"])

#%%
##Download the first 10 sequences
##Create Filenames with accession numbers
##Create  accession numbers
acc_files = []
acc_numbers = []
for index,id in enumerate(record1["IdList"]):
    print(index,id)
    acc_files.append(str(id + ".gbk"))
    acc_numbers.append(id)
    if index == 9:
        break
print(acc_files)
print(acc_numbers)

#%%
Entrez.email = "cyrusmallon@gmail.com"
for i in range(0,10):
    if not os.path.isfile(acc_files[i]): #check whether specified path is an existing regular file or not.
        net_handle = Entrez.efetch(
            db="nuccore",
            id=acc_numbers[i],
            rettype="gbwithparts",
            retmode="text"
        )
        out_handle = open(acc_files[i], "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print("Saved")

#%%
#download seq AND add them to a dictionary
seq_dict = {}
#Entrez.email = "cyrusmallon@gmail.com"
#for i in range(0,10):
#    if not os.path.isfile(acc_files[i]):
#        net_handle = Entrez.efetch(
#             db="nuccore",
#             id=acc_numbers[i],
#             rettype="gbwithparts",
#             retmode="text"
#         )
#         out_handle = open(acc_files[i], "w")
#         out_handle.write(net_handle.read())
#         out_handle.close()
#         net_handle.close()
#         print("Saved")
#        # seq_dict.update(SeqIO.parse(str(acc_files[i], "genbank"))

#%%
##add SeqObjects to a dictionary
seq_dict = {}
for i in range(len(acc_files)):
    print(i)
    seq_dict.update(SeqIO.to_dict(SeqIO.parse(acc_files[i], "genbank")))

#%%
chbOP = SeqIO.read("chbOP.fasta", "fasta")

#%%
l1 = []
#seq_dict = {}
for i in acc_numbers:
    print(i)
    l1.append(seq_dict[i])

#%%
l1.append(chbOP)

#%%
#make everything in a fasta file
fasta_file = SeqIO.write(l1, "seqs.fasta", "fasta")

#%%













#%%Add chbOP seq to dictionary
chbOP_dict = SeqIO.to_dict(SeqIO.parse("chbOP.fasta", "fasta"))
seq_dict.update(chbOP_dict)
#%%
if not os.path.isfile("NZ_CP012804.1.gbk"):
    net_handle = Entrez.efetch(
        db="nucleotide",
        id="NZ_CP012804.1",
        rettype="gbwithparts",
        retmode="text"
    )
    out_handle = open("NZ_CP012804.1.gbk", "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print("Saved")


#%%
##What type of instance is this part of the dictionary?
isinstance(record1["IdList"],list)
##it is indeed a list

#%%




#%%
from Bio import SeqIO

#%%
for seq_record in SeqIO.parse("B-Ancestral.gbk", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

#%%
##This doesn't work and returns an iterator.
##That's what it's supposed to return...
seqs = SeqIO.parse("B-Ancestral.gbk", "genbank")
#%%Add ancestral strain to dictionary
#But the below adds all scafolds...need to get them into a single genbank file...
seq_dict.update(SeqIO.to_dict(SeqIO.parse("B-Ancestral.gbk", "genbank")))

#check length of dictionary
len(seq_dict)

#check dictionary keys
list(seq_dict.keys())

#if you really want to look at the entire output
list(seq_dict.values())

#access a single record
seq_dict["chpOP"]