#%%write the sequence file
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#%%
#all ancestral E. coli and chpOP
with open('output_fasta.fa', 'w') as handle:
    SeqIO.write(seq_dict.values(), handle, 'fasta')
