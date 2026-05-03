# GWAS signal annotation resources

After identifying GWAS loci, annotation helps answer practical questions:

```text
Which gene is nearby?
Is the variant coding or non-coding?
Is it common or rare in public reference datasets?
Does it affect splicing, regulation, or expression?
Has it been associated with other traits before?
```

## Useful resources

| Resource | Main use |
|---|---|
| [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) | Variant Effect Predictor; coding consequence, transcript consequence, gene annotation |
| [PhenoScanner](http://www.phenoscanner.medschl.cam.ac.uk/) | Search previous genotype–phenotype associations |
| [SpliceAI Lookup](https://spliceailookup.broadinstitute.org/) | Predict potential splice-altering variants |
| [gnomAD](https://gnomad.broadinstitute.org/) | Population allele frequency and constraint metrics |
| [ENCODE SCREEN](https://screen.wenglab.org/) | Regulatory elements and enhancer/promoter annotations |
| [BRAVO](https://bravo.sph.umich.edu/) | Variant frequency browser based on TOPMed data |
| [UCSC Genome Browser](https://www.genome.ucsc.edu/) | Genome browser tracks, conservation, regulatory annotation |
| [All of Us Data Browser](https://databrowser.researchallofus.org/) | Public browsing of All of Us survey, EHR, and genomic data summaries |
| [dbSNP](https://www.ncbi.nlm.nih.gov/snp) | rsID lookup and variant metadata |

## Suggested annotation order

1. Confirm the genome build and variant identity.
2. Check whether the variant is coding, splice-region, or non-coding.
3. Check allele frequency in multiple reference datasets.
4. Search previous phenotype associations.
5. Review regional genes and regulatory elements.
6. Use locus-level evidence such as fine-mapping, eQTL/pQTL colocalization, or MR when available.
