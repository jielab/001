# VCF-to-protein and AlphaFold-related notes

This page records a compact workflow for connecting genetic variation, transcript/protein sequence changes, and protein structure prediction.

## 1. Input resources

Useful resources:

- [1000 Genomes high-coverage VCF files](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)
- [Ensembl reference FASTA](https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna)
- [GENCODE GTF/GFF3 files](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/)
- [Ensembl GTF files](https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/)
- [vcf2prot](https://github.com/ikmb/vcf2prot)

Typical preparation:

```bash
samtools faidx reference.fa
```

## 2. Conceptual workflow

```text
VCF variants
    ↓
Reference genome + transcript annotation
    ↓
Altered coding sequence
    ↓
Protein sequence
    ↓
Protein structure prediction or comparison
```

## 3. cDNA, mRNA and protein sequence

For AlphaFold-style protein structure prediction, the input is usually an amino-acid sequence, not genomic DNA. If starting from DNA, the relevant sequence is the coding transcript/cDNA sequence derived from mRNA, which does not contain introns.

Tools such as SnapGene or EditSeq can convert coding DNA sequence into protein sequence, but transcript selection and variant consequence annotation should be checked carefully.

## 4. Protein structure comparison

For comparing protein 3D structures, including structures from AlphaFold or cryo-EM, a common tool is TM-align.

## 5. Remove duplicate protein sequences

Many `_1` and `_2` protein sequences may be identical. Before running AlphaFold or other expensive prediction tools, remove duplicates.

```bash
seqkit rmdup input.fa > input.rmdup.fa
```
