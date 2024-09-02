from Bio import SeqIO
import os

SNP_dir = f'{os.getcwd()}/../1-genesel/intermediate'

with open(f'summary.txt','w') as outsum:
    print(f'Sample,Pct_passed,N_passed,N_failed',file=outsum)

for file in os.listdir(SNP_dir):
    if file.endswith('.processed.fasta'):
        isol = file.split('.')[0]
        with open(os.path.join(SNP_dir, file), 'r') as f, open(f'{isol}_gene_ambs.txt', 'w') as out, open(f'summary.txt','a') as outsum:
            records = list(SeqIO.parse(f, 'fasta'))
            print(f'Gene,PercAmb,NAmbiguous,NTotal,Notes', file=out)
            passed = 0
            allgenes = len(records)
            print(len(records))
            for record in records:
                amb = sum('?' in base for base in record.seq)
                if amb/len(record.seq) < 0.9 and amb!=0:
                    print(f'{record.id},{amb/len(record.seq)},{amb},{len(record.seq)}', file=out)
                    passed +=1
                else:
                    print(f'{record.id},{amb/len(record.seq)},{amb},{len(record.seq)},NULL gene', file=out)
            print(f'{isol},{passed/allgenes},{passed},{passed-allgenes}', file=outsum)






