# trio-scripts
A repository of scripts associated with the <i>Neotoma</i> trio binning project.

## GFF3 commands
These scripts were used to reformat the BRAKER GFF3 and assign gene names.

#### Add a header and correct formatting issues with the BRAKER GFF3:
```
GFF3/braker_gff3_reformatter.py Neotoma_bryanti.Contigs.fasta braker.gff3 | gt gff3 -sort yes -tidy yes -checkids yes -retainids yes > Neotoma_bryanti.BRAKER.gff3
```

#### Generate protein sequences and run InterProScan
```
GFF3/seqextractor.py -f Neotoma_bryanti.Contigs.fasta -g Neotoma_bryanti.BRAKER.gff3 -t protein -i all --no_stops -o Neotoma_bryanti.BRAKER.Protein.No_Stops.fasta
interproscan -i Neotoma_bryanti.BRAKER.Protein.No_Stops.fasta -dp -t p --goterms --pa -f TSV -cpu 24
```

#### Remove BRAKER gene predictions that lack IPR assignments
```
GFF3/filter_gff3_by_interproscan.py Neotoma_bryanti.BRAKER.Protein.No_Stops.fasta.tsv Neotoma_bryanti.BRAKER.gff3 > Neotoma_bryanti.BRAKER.IPR.gff3
```

#### Retrieve protein sequences and align them against other species with BLAST
```
GFF3/seqextractor.py -f Neotoma_bryanti.Contigs.fasta -g Neotoma_bryanti.BRAKER.IPR.gff3 -t protein -i all -o Neotoma_bryanti.BRAKER.IPR.Protein.fasta
blastp -db GFF3/Combined.Protein.fasta -query Neotoma_bryanti.BRAKER.IPR.Protein.fasta -out Neotoma_bryanti.BRAKER.IPR.Protein.outfmt6 -evalue 1e-5 -outfmt 6 -num_threads 24
```

#### Annotate the GFF3 with InterPro group IDs and BLAST descriptions
```
GFF3/annotate_gff3.py Neotoma_bryanti.BRAKER.IPR.gff3 Neotoma_bryanti.BRAKER.IPR.Protein.fasta Neotoma_bryanti.BRAKER.Protein.No_Stops.fasta.tsv GFF3/Combined.Descriptions.txt Neotoma_bryanti.BRAKER.IPR.Protein.outfmt6 NBRY > Neotoma_bryanti.Contigs.gff3
```

## Scaffolding commands
These scripts were used to generate the chromosome scaffold FASTAs and transform the GFF3 coordinates to match.

#### Generate a FASTA containing all contigs for Nbry_chromosome_1 and merge that into a single sequence with 500 N characters separating each contig ("\*" denotes a contig in reverse orientation).
```
Scaffolding/extract_and_orient_sequences.py Neotoma_bryanti.Contigs.fasta "Nbry_tig00000890" "*Nbry_tig00136238" "Nbry_tig00000726" "Nbry_tig00001382" "Nbry_tig00001371" "*Nbry_tig00000015" "*Nbry_tig00136208" "*Nbry_tig00136204" > Neotoma_bryanti.Contigs.chromosome_1.fasta
Scaffolding/multi-fasta_joiner.py Neotoma_bryanti.Contigs.chromosome_1.fasta Nbry_chromosome_1 500 > Neotoma_bryanti.chromosome_1.fasta
```

#### Transform the <i>Neotoma bryanti</i> contig GFF3 to match the chromosome coordinates.
```
Scaffolding/sequences_to_chromosomes.py 500 Scaffolding/Maps/Neotoma_bryanti.txt Neotoma_bryanti.Contigs.fasta Neotoma_bryanti.Contigs.gff3 > Neotoma_bryanti.Chromosomes.gff3
```

## MUMmer commands
These scripts were used to generate SNP and indel metrics based on MUMmer alignents.

#### Align the <i>Neotoma</i> mitochondrial sequences against each other
```
nucmer --prefix N_lepida_vs_N_bryanti.Mitochondrion Neotoma_bryanti.Mitochondrion.fasta Neotoma_lepida.Mitochondrion.fasta
delta-filter -r -q N_lepida_vs_N_bryanti.delta > N_lepida_vs_N_bryanti.Mitochondrion.filter
show-snps -Clr N_lepida_vs_N_bryanti.Mitochondrion.filter > N_lepida_vs_N_bryanti.Mitochondrion.Repeats.snps
```

#### Determine the amount of sequence aligned
```
MUMmer/mummer_alignment_length.py N_lepida_vs_N_bryanti.Mitochondrion.filter
```

#### Tally the number of SNPs and indels
```
MUMmer/mummer_snp_indel_counter.py N_lepida_vs_N_bryanti.Mitochondrion.Repeats.snps
```
