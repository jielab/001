# pip install alphagenome biopython
# 启明 module load python/anaconda3/2020.7; source activate; conda activate alphagenome

# %% -------------------- imports & setup --------------------
import os
import re
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from alphagenome.models import dna_client
from alphagenome.data import genome

# %% 
dir0 = Path("D:/")
API_KEY = os.getenv("ALPHA_GENOME_API_KEY") or "XXXXXXXX"
model = dna_client.create(API_KEY)


# %% -------------------- parameters & paths --------------------
prot = "ABO" 
CHR = "chr9" 
prot_START = 133_233_278
prot_END   = 133_276_024
fa_path  = dir0 / "analysis/vcf2prot" / prot / "dna/HG00096.H1.fa"
gtf_path = dir0 / "data/annot/ensembl/Homo_sapiens.GRCh38.115.gtf.gz"
out_fa  = dir0 / "analysis/vcf2prot" / prot / "alphaG/HG00096.H1.aa.fa"
out_fa.parent.mkdir(parents=True, exist_ok=True)
assert fa_path.exists(), f"FASTA not found: {fa_path}"
assert gtf_path.exists(), f"GTF not found:   {gtf_path}"


# %% -------------------- load haplotype & build window --------------------
rec = next(SeqIO.parse(str(fa_path), "fasta"))
hap = str(rec.seq).upper().replace("U", "T")
hap_len = len(hap)
expected_len = prot_END - prot_START + 1
print(f"Haplotype length: {hap_len:,} bp")
if hap_len != expected_len:
    print(f"[warn] H1 length ({hap_len}) != {prot} expected span ({expected_len}). Proceeding anyway.")
L = dna_client.SEQUENCE_LENGTH_1MB
center = (prot_START + prot_END) // 2
iv_start = center - (L // 2)
iv_end   = iv_start + L
interval = genome.Interval(CHR, iv_start, iv_end)
left_pad_len = prot_START - iv_start
assert 0 <= left_pad_len < L, "Left pad out of range; check coords/length L."
right_pad_len = L - left_pad_len - hap_len
assert right_pad_len >= 0, "Haplotype longer than window; increase L."
seq_for_model = "N" * left_pad_len + hap + "N" * right_pad_len
assert len(seq_for_model) == L
print(f"Interval: {interval.chromosome}:{interval.start:,}-{interval.end:,} (length={L:,})")
print(f"Placed haplotype at hg38 {CHR}:{prot_START:,}-{prot_START + hap_len - 1:,}")


# %% -------------------- run AlphaGenome (RNA-seq) --------------------
requested_outputs = [dna_client.OutputType.RNA_SEQ]
out = model.predict_sequences(
    sequences=[seq_for_model],
    organism=dna_client.Organism.HOMO_SAPIENS,
    requested_outputs=requested_outputs,
    ontology_terms=None,
    intervals=[interval],
    progress_bar=True,
)[0]
print("RNA-seq array shape:", out.rna_seq.values.shape)
print(out.rna_seq.metadata[['name','ontology_curie','biosample_name']].head())


# %% -------------------- GTF loader + attribute parser --------------------
COLS = ["seqname","source","feature","start","end","score","strand","frame","attributes"]
def load_gencode_gtf(gtf_path: str | Path) -> pd.DataFrame:
    return pd.read_csv(
        gtf_path, sep="\t", comment="#", header=None,
        names=COLS, compression="infer", low_memory=False
    )
_attr_re = re.compile(r'(\S+)\s+"([^"]+)"')
def parse_attr(attr: str) -> dict:
    return {k: v for k, v in _attr_re.findall(attr)}
gtf = load_gencode_gtf(gtf_path)
gtf = gtf[(gtf["feature"] == "CDS") & (gtf["seqname"].isin([CHR, CHR.replace("chr", "")]))].copy()
attrs = gtf["attributes"].apply(parse_attr)
gtf["gene_name"]     = attrs.apply(lambda d: d.get("gene_name"))
gtf["transcript_id"] = attrs.apply(lambda d: d.get("transcript_id"))
gtf["frame"]         = pd.to_numeric(gtf["frame"], errors="coerce")
cds = gtf[gtf["gene_name"] == prot].copy()
assert len(cds), f"No CDS entries found for {prot}."
cds = cds[(cds["end"] >= iv_start) & (cds["start"] <= iv_end)].copy()


# %% -------------------- build AA sequence(s) per transcript --------------------
def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())
def extract_hap(seg_start: int, seg_end: int) -> str:
    # GTF end inclusive -> slice end-exclusive +1
    i0 = max(seg_start, iv_start) - iv_start
    i1 = min(seg_end,   iv_end)   - iv_start + 1
    if i0 < 0 or i1 > len(seq_for_model) or i0 >= i1:
        return ""
    return seq_for_model[i0:i1]
def clip_by_frame(piece: str, frame, strand: str) -> str:
    if pd.isna(frame) or frame == 0:
        return piece
    frame = int(frame)
    return piece[frame:] if strand == "+" else (piece[:-frame] if frame <= len(piece) else "")
aa_records = []
for tid, sub in cds.groupby("transcript_id", sort=False):
    strand = sub["strand"].iloc[0]
    sub = sub.sort_values("start", ascending=(strand == "+")).copy()
    pieces = []
    for _, row in sub.iterrows():
        piece = extract_hap(int(row["start"]), int(row["end"]))
        piece = clip_by_frame(piece, row["frame"], strand)
        pieces.append(piece)
    cds_seq = "".join(pieces)
    if strand == "-":
        cds_seq = revcomp(cds_seq)
    if len(cds_seq) % 3 != 0:
        drop = len(cds_seq) % 3
        print(f"[warn] {tid}: CDS length not multiple of 3; trimming {drop} base(s).")
        cds_seq = cds_seq[:len(cds_seq) - drop]
    aa = str(Seq(cds_seq).translate(to_stop=False))
    aa_records.append((tid, strand, len(aa), aa))
aa_df = pd.DataFrame(aa_records, columns=["transcript_id","strand","aa_len","aa_seq"]).sort_values("aa_len", ascending=False)
print(aa_df.head(10))


# %% -------------------- write protein FASTA --------------------
with out_fa.open("w") as fh:
    for _, row in aa_df.iterrows():
        header = f">{row.transcript_id}|{prot}|len={row.aa_len}"
        fh.write(header + "\n")
        for i in range(0, len(row.aa_seq), 60):
            fh.write(row.aa_seq[i:i+60] + "\n")
print("Wrote:", out_fa)
