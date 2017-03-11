from subprocess import Popen, PIPE
from warnings import warn
from sys import argv
import glob
from Bio import SeqIO

job = int(argv[1])-1
rp_id = argv[2]

inputs = list(glob.glob("/home/457031/rootknot_phylogenomics/raw/*.fasta"))
print(len(inputs))
OGs = [i.rpartition('/')[-1].replace(rp_id+'_','').replace('.fasta','') for i in inputs]
sorted_OGs = sorted(OGs)[4000:]
OG = sorted_OGs[job]

# Run mafft
input_f = "/home/457031/rootknot_phylogenomics/raw/{0}_{1}.fasta".format(rp_id, OG)
output_f = "/home/457031/rootknot_phylogenomics/non_cds_alignment/{0}_{1}.aln.fasta".format(rp_id, OG)

mafft = "/home/457031/mafft-7.299-with-extensions/core/mafft"
mafftcline = mafft + " --thread 28 --localpair --maxiterate 1000 {0} > {1}".format(input_f, output_f)

input_f = "/home/457031/rootknot_phylogenomics/non_cds_alignment/{0}_{1}.aln.fasta".format(rp_id, OG)
output_f = "/home/457031/rootknot_phylogenomics/non_cds_alignment-trimal0.7/{0}_{1}.aln.trm.fasta".format(rp_id, OG)

#p = Popen(mafftcline, shell=True, stdout=PIPE, stderr=PIPE)
#out, err = p.communicate()

#print(out)
#print(err)

# Run trimal

trimal = "/home/457031/trimAl/source/trimal"
trimalcline = trimal + " -in {0} -out {1} -gt 0.7 -st 0.001".format(input_f,output_f)

p = Popen(trimalcline, shell=True, stdout=PIPE, stderr=PIPE)
out, err = p.communicate()

print(out)
print(err)

# Check filter out alignments in which a single gene copy is split between contigs

def is_bad_alignment(aln):
    badaln = False
    alnlength = aln.get_alignment_length()
    for i in range(len(aln)):
        if badaln:
            break
        iseqrecord = aln[i]
        isp = iseqrecord.id.split('_')[0]
        for j in range(i+1, len(aln)):
            jseqrecord = aln[j]
            jsp = jseqrecord.id.split('_')[0]
            if isp == jsp:
                overlap = 0
                for x in range(alnlength):
                    if aln[i,x] != '-' and aln[j,x] != '-':
                        overlap += 1
                if overlap < 20:  
                    badaln = True
                    break
    return badaln


def check_alignment_file(f):

    from Bio import AlignIO
    aln = AlignIO.read(f,'fasta')
    aln_name = f.split('/')[-1].split('_')[0]
    if is_bad_alignment(aln):
        return aln_name
    else:
        return None

input_f = "/home/457031/rootknot_phylogenomics/non_cds_alignment/{0}_{1}.aln.fasta".format(rp_id,OG)
output_f = input_f.replace("non_cds_alignment", "non_cds_alignment-trimal0.7-badaln")

aln_name = check_alignment_file(input_f)

if aln_name:
    with open(output_f, 'wt') as hndl:
        hndl.write(aln_name)

baddies_path = "/home/457031/rootknot_phylogenomics/non_cds_alignment-trimal0.7-badaln/*_*.aln.fasta"
bad_alns = [pt.rpartition('/')[-1].split('_')[1].split('.')[0] for pt in list(glob.glob(baddies_path))]

# Run RAxML

raxmlpthreads = "/home/457031/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3"

input_f = "/home/457031/rootknot_phylogenomics/non_cds_alignment-trimal0.7/{0}_{1}.aln.trm.fasta".format(rp_id,OG)
raxmlcline = raxmlpthreads + " -w /home/457031/rootknot_phylogenomics/non_cds_alignment-trimal0.7-tree  -p 123 -s {0} -T 28 -n {1}_{2} -m GTRGAMMA".format(input_f, rp_id, OG)

if not OG in bad_alns:
    try:
        p = Popen(raxmlcline, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
    except:
        warn(OG + " tree failed")
