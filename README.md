
 OXPHOS Gene Extraction and Probe Design Pipeline

This project identifies and extracts ~80 OXPHOS genes from the North American Mallard genome using a reference from Pekin duck, maps them via BLAST, and designs 120mer hybridization probes for targeted sequencing.

##  Step-by-Step Instructions

### 🔹 Step 1: Prepare OXPHOS Gene List
Create a file `oxphos_gene_list.txt` from KEGG pathway (map00190) with gene symbols like:
```
NDUFA6
COX4I1
ATP5F1A
...
```
🗂 File also saved as: `OXPHOS_gene_List.csv`

---

###  Step 2: Download Duck Genomes
```bash
# Create directory
mkdir -p pekin_duck_annotation && cd pekin_duck_annotation

# Download and unzip Pekin duck genome and annotations
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/663/525/GCF_047663525.1_IASCAAS_PekinDuck_T2T/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/663/525/GCF_047663525.1_IASCAAS_PekinDuck_T2T/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/663/525/GCF_047663525.1_IASCAAS_PekinDuck_T2T/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gff.gz
gunzip *.gz
```

---

###  Step 3: Extract GFF Hits for OXPHOS Genes
```bash
mkdir -p gff_matches

while read gene; do
  echo "Searching for $gene..."
  grep -i -w "$gene" pekin_duck_annotation/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gff >> gff_matches/${gene}.gff
done < oxphos_gene_list.txt

cd gff_matches
cat *.gff > ../oxphos_combined.gff
```

Check the combined file:
```bash
wc -l ../oxphos_combined.gff
head ../oxphos_combined.gff
```

---

###  Step 4: Create BED File and Extract Sequences
```bash
awk '$3 == "gene" {print $1"	"($4-1)"	"$5"	"$9"	.	"$7}' ../oxphos_combined.gff > oxphos_genes.bed

bedtools getfasta   -fi pekin_duck_annotation/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna   -bed oxphos_genes.bed   -s -name   -fo oxphos_gene_sequences.fasta

grep ">" oxphos_gene_sequences.fasta | wc -l
```

---

###  Step 5: Map to NAwild Genome with BLAST
```bash
# Prepare BLAST database
makeblastdb -in GCA_030704485.1_NAwild_v1.0_genomic.fna -dbtype nucl -out nawild_db

# Run BLAST
blastn -query oxphos_gene_sequences.fasta   -db nawild_db   -out oxphos_vs_nawild.blastn.out   -evalue 1e-10   -outfmt 6   -num_threads 4

# Extract useful columns
cut -f1,2,9,10 oxphos_vs_nawild.blastn.out | sort | uniq > oxphos_gene_hits_coords.txt

# Convert to BED format
awk '{start=($3<$4)?$3:$4; end=($3>$4)?$3:$4; print $2"\t"start-1"\t"end"\t"$1}' oxphos_gene_hits_coords.txt > oxphos_blast_hits.bed
```

---

###  Step 6: Design 120mer Probes
```bash
# Tile 120mer windows with 60 bp overlap
bedtools makewindows -b oxphos_blast_hits.bed -w 120 -s 60 > oxphos_120mer_windows.bed

# Extract probe sequences
bedtools getfasta   -fi GCA_030704485.1_NAwild_v1.0_genomic.fna   -bed oxphos_120mer_windows.bed   -s -name   -fo oxphos_120mer_probes.fasta

# Count total number of probes
grep ">" oxphos_120mer_probes.fasta | wc -l
```

---

###  Step 7: Rename Probes for Synthesis
```bash
awk 'NR % 2 == 0' oxphos_120mer_probes.fasta > seqs.txt
seq 1 $(wc -l < seqs.txt) | awk '{printf ">probe_%05d\n", $1}' > headers.txt
paste -d "\n" headers.txt seqs.txt > oxphos_120mer_probes_synthesis.fasta
```

---

###  Step 8: Link Probes to Coordinates (Metadata File)
```bash
awk '{print $1"\t"$2"\t"$3"\t"$6}' oxphos_120mer_windows.bed > probe_coords.txt
paste headers.txt probe_coords.txt > probe_metadata.tsv
```

---

You now have:
- `oxphos_120mer_probes_synthesis.fasta`: probes for hybridization
- `probe_metadata.tsv`: probe ID to genome coordinate table
- `oxphos_gene_sequences.fasta`: extracted gene sequences

