README.txt

The following files accompany the paper:

Lunt DH, Kumar S, Koustovoulos GD, Blaxter ML, 2014. The complex hybrid origins of the root knot nematodes revealed through comparative genomics. PeerJ

- mh.cds.fna nucleotide fasta file for transcripts from Meloidogyne hapla
- mi.cds.fna nucleotide fasta file for transcripts from Meloidogyne incognita
- mf.cds.fna nucleotide fasta file for transcripts from Meloidogyne floridensis
- mh.protein.faa aa fasta file for proteins from Meloidogyne hapla
- mi.protein.faa aa fasta file for proteins from Meloidogyne incognita
- mf.protein.faa aa fasta file for proteins from Meloidogyne floridensis
- InParanoid-mh-mi-mf.tgz results of running InParanoid on each pair of species
- QuickParanoid-mh-mi-mf.tgz results of running QuickParanoid on InParanoid output to get three-species clusters.
- raxml-mh-mi-mf.tgz - tar gz file of raxml results. Contains 7 subfolders with different numbers of proteins from each species in eah protein cluster, corresponding to the cells in Figure 4. Each subfolder has the following files:
    - faa.tgz protein fasta files, one file for each cluster 
    - fna.tgz nucleotide fasta files, one file for each cluster 
    - clustalo.tgz clustalo protein alignments for each cluster
    - tranalign.tgz tranalign nucleotide alignments for each cluster using clustalo alignments as guides
    - phy.tgz above alignments in phy format (for input to raxml)
    - RAxML.bipartitions output of raxml on .phy alignments above, concatenated into one file
