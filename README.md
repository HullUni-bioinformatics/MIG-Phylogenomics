# MIG-Phylogenomics

### 0. Dependencies
Notebook file name:  [`Dependencies.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/Dependencies.ipynb)
#### Related files:
None
### 1. CDSs and proteins from genome assemblies
Notebook file name: [`CDSs_and_proteins_from_genome_assemblies.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/CDSs_and_proteins_from_genome_assemblies.ipynb)
#### Related files:
`meloidogyne_assemblies`: contains fasta genome assemblies  
`annotation`: contain gff files for the assemblies in `assemblies`   
  
`<None | stopped | all>_<gene | cds | protein>_ref_<files | centroids | reviewed> `   
with `None` indicating that nothing is written.  
  
+  dirs that *start* with `None`: genes, cdss or proteins without premature stop codon
+  dirs that *start* with `stopped`: genes, cdss or proteins with a premature stop codon
+  dirs that *start* with `all`: a merge of `None` and `stopped`  
+  dirs that *end* with `files`: raw, as indicated in the gff
+  dirs that *end* with `centroinds`: cds files that were reduced with a vsearch step
+  dirs that *end* with `reviewed`: final treated datasets (see notebook)
+  `ref` in all the dir names indicate that these files are derived from a genome assembly annotation.

### 2. Map-assemble genes from read data for samples without assemblies
Notebook file name: [`Map_assemble_gene.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/Map_assemble_genes.ipynb)
#### Related files
`<sample name>_bwa/<sample name>.nt.fasta`: map-assembled gene files  
  
`<None | stopped | all >_<cdss | proteins | gffs>`   
with `None` indicating that nothing is written.  
  
+  dirs that *start* with `None`: gffs, cdss or proteins without premature stop codon
+  dirs that *start* with `stopped`: gffs, cdss or proteins with a premature stop codon
+  dirs that *start* with `all`: a merge of `None` and `stopped`  
  
### 3. Orthology clustering
Notebook file name: [`Orthology_clustering.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/Orthology_clustering.ipynb)
#### Related files:
`orthofinder/all_inputs/<sample name><None|_ref>.aa.fasta`: links to protein sequences of all the samples. They will need to be regenerated locally (step included in the notebook).   
  
`orthofinder/all_inputs/Results_Jul02/<inflation value>_OrthologousGroup.csv`: Orthology clusters, with `<inflation value>` representign the mcl inflation parameter, except for 0, representing an inflation of 1.5, and 1, representing inflation of 1.1.   
  
`orthofinder/all_inputs/Results_Jul02/WorkingDirectory`: OrthoFinder inputs and putputs of the Blast step.  
### 4. Nuclear phylogenomics
Notebook file name: [`Nuclear_phylogenomics.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/Nuclear_phylogenomics.ipynb)
#### Related files:
`orthofinder/all_inputs/Results_Jul02/OGs_I2_1-4.gb.gz`: A genbank file with coding and protein sequences of orthology clusters with 1 to 4 gene copies for each reference sample`.  
  
`orthofinder/all_inputs/Results_Jul02/OGs_I2_1-4.gb.loci.<csv|txt>`: ReproPhylo formated list of the loci that are in the genbank file.
  
`orthofinder/all_inputs/Results_Jul02/I2_rootknot_phylogenomics`: Input and output files of the gene tree phylogenetic pipeline, with trimal settings of gt=0.7 and st=0.01`  
  
`orthofinder/all_inputs/Results_Jul02/I2_3X2_gt0.7_st_0.01_alns_<1-4 | all2 | flo2>`: Sequence alignments of orthology clusters in which inparalogues are collapsed into a single sequence.  
`1-4`: all the orthology clusters  in which there are at least 3 reference samples with 2 gene copies.
`all2`: a subset of `1-4` in which all the reference samples have two gene copies.  
`flo2`: a subset of `1-4` in which all MfloSJF1 has two gene copies.  
    
`orthofinder/all_inputs/Results_Jul02/I2_3X2_gt0.7_st_0.01_alns_1-4/<astralshuffeled | raxmlshuffled>`: randomization analyses in which homeologue 1 and homeologue 2 are randomly determined for each gene.  
`astralshuffeled`: 100 astral runs, in which hom 1 and 2 were randomly determined for each gene.   
`raxmlshuffled`: 100 raxml supermatrix trees, in which hom 1 and 2 were randomly determined for each gene, prior to the concatenation of the supermatrix.  
  
`orthofinder/all_inputs/Results_Jul02/I2_3X2_gt0.7_st_0.01_alns_1-4/trees.txt`: a list of gene trees that were used for astral (non randomized)  
  
`orthofinder/all_inputs/Results_Jul02/I2_3X2_gt0.7_st_0.01_alns_1-4/raxmlshuffled/trees.txt`: a list of randomized supermatrix trees.     
  
`orthofinder/all_inputs/Results_Jul02/I2_3X2_gt0.7_st_0.01_alns_1-4/ RAxML_StrictConsensusTree<AstStrict | RaxStrict>`: strict consensus trees that resulted from the two randomization analyses with astral and raxml.  

`orthofinder/all_inputs/Results_Jul02/I2_3X2_gt0.7_st_0.01_alns_1-4/ RAxML_bipartitions.f2f51893cb434759411b04a1ca40a4b9ad85c46b_2`:  
A through raxml tree reconstrction of a supermatrix as produced by treeCl.  
 
### 5. Mitochondrial genome assembly
Notebook file name:[ `Mitochondrial_genome_assembly.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/Mitochondrial_genome_assembly.ipynb)
### 6. Mitochondrial genome annotation
Notebook file name: [`Mitochondrial_genomes_annotation.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/Mitochondrial_genomes_annotation.ipynb)
### 7. Mitochondrial genome phylogenomics
Notebook file name: [`Mitochondrial_genomes_tree.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/Mitochondrial_genomes_tree.ipynb)
### 8. Intra-genome identity among homeologue gene pairs
Notebook file name:[ `Intra_genome_sequence_divergence.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/Intra_genome_sequence_divergence.ipynb)
### 9. Coverage ratio between homeologue contigs within a genome
Notebook file name: [`Median_ratio.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/Median_ratio.ipynb)
### 10. Gene conversion in *M. floridensis* MfloSJF1
Notebook file name: [`Non-conversion_tracts_in_MfloSJF1s.ipynb`](https://github.com/HullUni-bioinformatics/MIG-Phylogenomics/blob/master/Non-conversion_tracts_in_MfloSJF1.ipynb)