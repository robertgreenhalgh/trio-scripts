# trio-scripts
A repository of scripts associated with the <i>Neotoma</i> trio binning project.

## GFF3 generation commands
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

