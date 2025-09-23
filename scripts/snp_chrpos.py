#!/bin/python3
from itertools import compress
from typing import Optional, Iterable, Tuple, IO, Union, List, Dict, Set
from multiprocessing import cpu_count, Pool
from gzip import open as gzip_open
from time import time
from tempfile import mkdtemp
import os
import argparse
from collections import OrderedDict

# -------------------- Globals (mostly kept from original) --------------------
sep_str = '\t'
NA_str = 'NA'
rs_id_str = 'SNP'
chrom_order_list = [str(i) for i in range(1, 22 + 1)] + ['X', 'Y', 'MT']
chrom_order = {chrom: index + 1 for index, chrom in enumerate(chrom_order_list)}
chrom_aliases = {'23': 'X', '24': 'Y', '25': 'MT', 'M': 'MT'}
# Names the script can recognize in input GWAS
columns_name = OrderedDict({
    'CHR': ['#CHROM', 'chr', 'CHR', 'chrom'],
    'POS': ['BEG', 'BEGIN', 'pos', 'POS'],
    'REF': ['reference', 'ref', 'REF'],
    'ALT': ['alt', 'alternate', 'ALT']
})

# -------------------- Exceptions & basic utils -------------------------------
class PAGEANTERROR(Exception):
    def __init__(self, *args):
        super(PAGEANTERROR, self).__init__(*args)

def open_(file: str, mode: str = 'rb', **kwargs) -> IO:
    if file.split('.')[-1].lower() == 'gz':
        return gzip_open(file, mode=mode, **kwargs)
    else:
        return open(file, mode=mode, **kwargs)

def rm_dir(rmdir: str) -> None:
    if os.path.isdir(rmdir):
        for file in os.listdir(rmdir):
            os.remove(os.path.join(rmdir, file))
        os.rmdir(rmdir)

def select_list(ob_list: list, index) -> list:
    return [ob_list[i] if i < len(ob_list) else None for i in index]

def equal_seq(seq1: str, seq2: str) -> bool:
    # treat 'N' as wildcard
    if seq1 == seq2:
        return True
    elif seq1 in seq2 or seq2 in seq1:
        return True
    elif len(seq1) == len(seq2) and ('N' in seq1 or 'N' in seq2):
        return all(b1 == b2 or b1 == 'N' or b2 == 'N' for b1, b2 in zip(seq1, seq2))
    else:
        return False

def get_ind(pos: List[int], pos_length: List[int]) -> List[Optional[int]]:
    idx_rev = []
    for p in pos:
        if p is None:
            idx_rev.append(None)
        else:
            i = 0
            while p >= pos_length[i]:
                p -= pos_length[i]
                i += 1
            idx_rev.append(i)
    return [idx_rev.index(i) if i in idx_rev else None for i in range(len(pos_length))]

def get_columns_index(col_names: List[str], need_cols: List[List[str]]) -> List[Optional[int]]:
    needs = [need for need_col in need_cols for need in need_col]
    need_pos = [len(need_col) for need_col in need_cols]
    need_idx = [needs.index(name) if name in needs else None for name in col_names]
    return get_ind(need_idx, need_pos)

def file_div(file: str, n_div: Optional[int] = None, start: int = 0, min_step: int = 1) -> Iterable[Tuple[int, int]]:
    if not n_div:
        n_div = cpu_count()
    with open_(file) as f:
        end = 0
        for _ in f:
            end += 1
    print(f'Total line of {file}: {end}')
    step = int(end / n_div + 1)
    step = step if step > min_step else min_step
    while start < end:
        if start + step >= end:
            yield start, end
        else:
            yield start, start + step
        start += step

def read_file_part(file: str, start: int, end: Optional[int] = None) -> Iterable[Tuple[int, bytes]]:
    if not end:
        end = float('inf')
    with open_(file) as f:
        i = 0
        for line in f:
            if i >= start:
                if i < end:
                    yield i, line
                else:
                    break
            i += 1

# -------------------- Genomic position classes (kept) ------------------------
class POS:
    def __init__(self, chrom: Union[str, int], pos: Union[str, int]):
        if isinstance(chrom, str):
            if chrom in chrom_order:
                self.chrom = chrom_order[chrom]
            else:
                try:
                    self.chrom = chrom_order[chrom_aliases[chrom]]
                except KeyError as e:
                    raise PAGEANTERROR(f'Cannot recognize chromosome "{e.args[0]}"!') from e
        else:
            self.chrom = chrom
        if isinstance(chrom, str):
            try:
                self.pos = int(pos)
            except ValueError as e:
                raise PAGEANTERROR(f'Cannot recognize position {e.args[0].split(":")[-1]}') from e
        else:
            self.pos = pos
        assert self.pos > 0, PAGEANTERROR(f'Bad position ({self.pos})')

    def __lt__(self, other):
        return (self.chrom, self.pos) < (other.chrom, other.pos)
    def __eq__(self, other):
        return (self.chrom, self.pos) == (other.chrom, other.pos)
    def __ne__(self, other):
        return not (self == other)
    def __repr__(self):
        return f"{self.chrom}:{self.pos}"
    def __str__(self):
        return repr(self)
    def tuple(self):
        return self.chrom, self.pos
    def export(self, *args: str, sep: str = '\t'):
        return f"{sep.join([str(self.chrom), str(self.pos)] + list(args))}"

class POS_Idx(POS):
    def __init__(self, pos: Union[POS, Tuple[Union[str,int], Union[str,int]]]):
        if isinstance(pos, POS):
            pos = pos.tuple()
        super(POS_Idx, self).__init__(*pos)
        self.pos //= 100000  # 100kb bin

# -------------------- Original indexer for pos->SNP (kept) -------------------
def get_index_part(file: str, start: int, end: int, index: int, output: str) -> None:
    with open(os.path.join(output, f'{index}.idx'), 'w') as f_idx:
        pos_t = POS_Idx((0, 1))
        for idx, line in read_file_part(file, start, end):
            line = line.decode().strip().split('\t')
            if pos_t < POS_Idx(line[0:2]):
                pos_t = POS_Idx(line[0:2])
                f_idx.write(pos_t.export(str(idx)) + '\n')

def get_index(file: str):
    temp_dir = mkdtemp(suffix='idx')
    with Pool(processes=cpu_count()) as pool:
        i = 0
        for start, end in file_div(file):
            i += 1
            pool.apply_async(get_index_part, args=(file, start, end, i, temp_dir))
        pool.close()
        pool.join()
    init = (0, 0)
    with open(file + '.idx', 'w') as f:
        for idx in range(1, i + 1):
            with open(os.path.join(temp_dir, f'{idx}.idx')) as f_div:
                start = f_div.readline()
                if tuple(int(num) for num in start.split('\t')[:2]) > init:
                    f.write(start)
                for line in f_div:
                    start = line
                    f.write(start)
    rm_dir(temp_dir)

# -------------------- Variant / line wrappers (extended) ---------------------
class Variant:
    def __init__(self, information: List[str] or Tuple[str, str, Optional[str], Optional[str]]):
        # information: [CHR, POS, REF, ALT] ; REF/ALT can be None
        self.pos = POS(*information[:2])
        self.ref = information[2] if len(information) > 2 else None
        self.alt = information[3] if len(information) > 3 else None

    def __eq__(self, other):
        if self.pos == other.pos:
            # if either side lacks alleles, treat as positional match
            if (self.ref is None or other.ref is None) and (self.alt is None or other.alt is None):
                return True
            # otherwise need allele-consistent match
            if equal_seq(self.ref, other.ref) and equal_seq(self.alt, other.alt):
                return True
            # allow swapped REF/ALT (strand flips are not handled here)
            if equal_seq(self.alt, other.ref) and equal_seq(self.ref, other.alt):
                return True
        return False

    def __lt__(self, other):
        return self.pos < other.pos

class LINE:
    def __init__(self, line: List[str], idx: Optional[List[int]] = None):
        if idx:
            self.variant = Variant(select_list(line, idx))
        else:
            self.variant = Variant(line[:4])
        self.line = line
        self.id: str = NA_str

    def __eq__(self, other):
        return self.variant == other.variant
    def __lt__(self, other):
        return self.variant < other.variant
    def export(self, sep: str = '\t'):
        return sep.join(self.line + [self.id])

class RSID_LINE(LINE):
    """
    dbSNP / mapping line wrapper.
    Accepts either:
      [CHR, POS, RSID]  or  [CHR, POS, RSID, REF, ALT]
    """
    def __init__(self, fields: List[str]):
        # Detect layout
        if len(fields) < 3:
            raise PAGEANTERROR("DB line with <3 columns encountered.")
        chrom = fields[0]
        pos   = fields[1]
        rsid  = fields[2]
        ref   = fields[3] if len(fields) >= 4 else None
        alt   = fields[4] if len(fields) >= 5 else None

        # Build variant as [CHR, POS, REF, ALT] (REF/ALT can be None)
        super(RSID_LINE, self).__init__([chrom, pos, ref, alt])
        self.id = rsid

# -------------------- pos->SNP (kept; now tolerant to 3/5-col DB) ------------
def get_line_start(pos: POS_Idx, db_file: str):
    with open(db_file + '.idx') as f:
        for line in f:
            if line.startswith(pos.export()):
                return int(line.strip().split(sep_str)[-1])

def get_rsid_part(raw_file: str, start: int, end: int, idx: int, db_file: str, temp_dir: str,
                  columns_idx: List[int], sep_str: str) -> None:
    with open(os.path.join(temp_dir, str(idx)), 'w', encoding='UTF-8') as f_write_part:
        rsid_group_reader = None
        rsid_group = None
        for _, line in read_file_part(raw_file, start, end):
            line = line.decode().strip().split(sep_str)
            process = LINE(line, columns_idx)
            # lazy init of db cursor near this bin
            if rsid_group_reader is None:
                rsid_group_reader = read_file_part(db_file, get_line_start(POS_Idx(process.variant.pos), db_file))
                try:
                    _, rsid_line = next(rsid_group_reader)
                    rsid_group = RSID_LINE(rsid_line.decode().strip().split('\t'))
                except StopIteration:
                    rsid_group = None
            # advance db cursor until >= process
            while rsid_group is not None and rsid_group < process:
                try:
                    _, rsid_line = next(rsid_group_reader)
                    rsid_group = RSID_LINE(rsid_line.decode().strip().split('\t'))
                except StopIteration:
                    rsid_group = None
                    break
            # scan equals (same position; optionally alleles)
            while rsid_group is not None and rsid_group.variant.pos == process.variant.pos:
                if rsid_group == process:
                    process.id = rsid_group.id
                    break
                else:
                    try:
                        _, rsid_line = next(rsid_group_reader)
                        rsid_group = RSID_LINE(rsid_line.decode().strip().split('\t'))
                    except StopIteration:
                        rsid_group = None
                        break
            f_write_part.write(process.export(sep=sep_str) + '\n')

def get_rsid(raw_file: str, db_file: str, output_file: str, processes: Optional[int] = None):
    time0 = time()
    if not processes:
        processes = cpu_count()
    print("Getting/Checking index of db file...")
    if not os.path.isfile(db_file + '.idx'):
        get_index(db_file)
    print('Index ready.')

    temp_dir = mkdtemp(suffix='rsid')
    with open_(raw_file) as f_raw, open_(output_file, 'wb') as f_write:
        header = f_raw.readline()
        f_write.write(header.strip() + sep_str.encode() + rs_id_str.encode() + b'\n')
        header = header.decode().strip().split(sep_str)
        columns_idx = get_columns_index(header, list(columns_name.values()))
        assert None not in columns_idx, \
            f"Cannot find columns " \
            f"{list(compress(['CHR','POS','REF','ALT'], map(lambda a: a is None, columns_idx)))} in input!"
    print("Start annotation (pos -> rsID)...")
    with Pool(processes=processes) as pool:
        idx = 0
        for start, end in file_div(raw_file, processes, 1):
            idx += 1
            pool.apply_async(get_rsid_part, args=(raw_file, start, end, idx, db_file, temp_dir, columns_idx, sep_str))
        pool.close(); pool.join()
    with open_(output_file, 'ab') as f_write:
        for i in range(1, idx + 1):
            with open(os.path.join(temp_dir, str(i)), 'rb') as f_write_part:
                for line in f_write_part:
                    f_write.write(line)
    rm_dir(temp_dir)
    print(f'Finish! Used time: {time() - time0:.0f}s')

# -------------------- NEW: SNP->pos (stream DB once) -------------------------
def parse_db_line_to_tuple(fields: List[str]) -> Tuple[str, str, str, Optional[str], Optional[str]]:
    """
    Returns (CHR, POS, RSID, REF, ALT); REF/ALT may be None.
    Accept 3-col or 5-col layouts.
    """
    if len(fields) < 3:
        raise PAGEANTERROR("DB line with <3 columns encountered.")
    chrom = fields[0]
    pos   = fields[1]
    rsid  = fields[2]
    ref   = fields[3] if len(fields) >= 4 else None
    alt   = fields[4] if len(fields) >= 5 else None
    return chrom, pos, rsid, ref, alt

def add_pos_from_rsid(raw_file: str, db_file: str, output_file: str,
                      sep: str, snp_colname: str,
                      na: str = 'NA') -> None:
    """
    Read input GWAS (must have a SNP column), stream db file to build a small mapping
    only for SNPs present in GWAS, then append CHR & POS (and REF/ALT if available).
    """
    t0 = time()
    # 1) Read header, collect needed SNPs
    with open_(raw_file) as f:
        header = f.readline().decode().rstrip('\n')
        cols   = header.split(sep)
        # locate SNP column (case-insensitive)
        snp_idx = None
        for i, c in enumerate(cols):
            if c.upper() == snp_colname.upper():
                snp_idx = i; break
        if snp_idx is None:
            raise PAGEANTERROR(f"SNP column '{snp_colname}' not found in input file.")
        need: Set[str] = set()
        for _, line in read_file_part(raw_file, 1, None):
            rs = line.decode().rstrip('\n').split(sep)[snp_idx]
            if rs and rs != na:
                need.add(rs)
    print(f"Need positions for {len(need):,} rsIDs")

    # 2) Stream db file and cache only needed rsIDs
    mapping: Dict[str, Tuple[str,str,Optional[str],Optional[str]]] = {}
    with open_(db_file) as fdb:
        # Detect header: if third token doesn't start with 'rs', assume header present -> skip
        first = fdb.readline()
        if not first:
            raise PAGEANTERROR("Empty DB file.")
        first_fields = first.decode().rstrip('\n').split('\t')
        has_header = (len(first_fields) >= 3 and not first_fields[2].startswith('rs'))
        if not has_header:
            # treat first line as data
            chrom, pos, rsid, ref, alt = parse_db_line_to_tuple(first_fields)
            if rsid in need and rsid not in mapping:
                mapping[rsid] = (chrom, pos, ref, alt)
        # continue streaming
        for _, b in enumerate(fdb):
            fields = b.decode().rstrip('\n').split('\t')
            if len(fields) < 3: continue
            chrom, pos, rsid, ref, alt = parse_db_line_to_tuple(fields)
            if rsid in need and rsid not in mapping:
                mapping[rsid] = (chrom, pos, ref, alt)
            # small early-exit if fully covered
            if len(mapping) == len(need):
                break
    print(f"Mapped {len(mapping):,} / {len(need):,} rsIDs")

    # 3) Write output: append CHR POS (and REF ALT if available)
    with open_(raw_file) as fin, open_(output_file, 'wb') as fout:
        header = fin.readline().decode().rstrip('\n')
        cols   = header.split(sep)
        has_refalt = any(v[2] is not None for v in mapping.values())
        extra_cols = ['CHR','POS'] + (['REF','ALT'] if has_refalt else [])
        fout.write((header + sep + sep.join(extra_cols) + '\n').encode())
        snp_idx = None
        for i, c in enumerate(cols):
            if c.upper() == snp_colname.upper():
                snp_idx = i; break
        for _, b in read_file_part(raw_file, 1, None):
            row = b.decode().rstrip('\n').split(sep)
            rs  = row[snp_idx]
            if rs in mapping:
                chrom, pos, ref, alt = mapping[rs]
                extra = [chrom, pos] + ([ref or na, alt or na] if has_refalt else [])
            else:
                extra = [na, na] + (['NA','NA'] if has_refalt else [])
            fout.write((sep.join(row + extra) + '\n').encode())
    print(f"Done SNP->POS in {time() - t0:.0f}s; output -> {output_file}")

# -------------------- CLI glue (two modes) -----------------------------------
def addrsid_run(raw_file: str, db_file: str, output_file: str, processes: Optional[int] = None,
                sep: str = '\t', chr_col: Optional[str] = None, pos_col: Optional[str] = None,
                ref: str = 'REF', alt: str = 'ALT',
                rs_id: str = 'SNP', na: str = 'NA',
                snp_col: Optional[str] = None) -> None:
    """
    If chr_col and pos_col are provided (optionally ref/alt): fill SNP
    If snp_col is provided: fill CHR & POS (and REF/ALT if available)
    """
    global columns_name, sep_str, NA_str, rs_id_str
    sep_str = eval(repr(sep).replace('\\\\', '\\'))
    NA_str  = na
    rs_id_str = rs_id

    # Mode detection & validation
    mode_pos2snp = (chr_col is not None and pos_col is not None)
    mode_snp2pos = (snp_col is not None)
    if mode_pos2snp and mode_snp2pos:
        raise PAGEANTERROR("Please choose ONE mode: either --snp OR (--chr and --pos).")
    if not mode_pos2snp and not mode_snp2pos:
        raise PAGEANTERROR("You must specify either --snp <col> OR --chr <col> --pos <col>.")

    if mode_pos2snp:
        # augment recognizable column names so get_columns_index() can find them
        columns_name['CHR'] = list(set(columns_name['CHR'] + [chr_col]))
        columns_name['POS'] = list(set(columns_name['POS'] + [pos_col]))
        columns_name['REF'] = list(set(columns_name['REF'] + [ref]))
        columns_name['ALT'] = list(set(columns_name['ALT'] + [alt]))
        get_rsid(raw_file, db_file, output_file, processes)
    else:
        # SNP -> POS path
        add_pos_from_rsid(raw_file, db_file, output_file, sep=sep, snp_colname=snp_col, na=na)

# -------------------- Main ---------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='snp_chrpos',
        description='Fill rsID from CHR/POS (and optional REF/ALT), or fill CHR/POS from rsID.'
    )
    parser.add_argument('-i', '--input', required=True, help="Input GWAS file (.txt/.tsv, optionally .gz)")
    parser.add_argument('-d', '--database', required=True, help="dbSNP/rsids mapping file (.tsv/.tsv.gz)")
    parser.add_argument('-o', '--output', required=True, default='./annotated.tsv.gz',
                        help="Output file (.tsv or .tsv.gz)")
    parser.add_argument('-s', '--sep', default='\t', help="Delimiter of input file (default: tab)")
    parser.add_argument('-c', '--core', default=None, help="Number of processes for pos->SNP mode")
    parser.add_argument('--na', default='NA', help="NA string in output")
    parser.add_argument('--rsid', default='SNP', help="Output column name for rsID when filling SNP")

    # Mode A: CHR/POS -> SNP
    parser.add_argument('--chr', dest='chr_col', help="Column name of chromosome in input GWAS")
    parser.add_argument('--pos', dest='pos_col', help="Column name of position in input GWAS")
    parser.add_argument('--ref', default='REF', help="Column name of reference allele in input (optional)")
    parser.add_argument('--alt', default='ALT', help="Column name of alternate allele in input (optional)")

    # Mode B: SNP -> CHR/POS
    parser.add_argument('--snp', dest='snp_col', help="Column name of SNP(rsID) in input GWAS")

    args = parser.parse_args()
    addrsid_run(
        raw_file=args.input,
        db_file=args.database,
        output_file=args.output,
        processes=args.core,
        sep=args.sep,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        ref=args.ref,
        alt=args.alt,
        rs_id=args.rsid,
        na=args.na,
        snp_col=args.snp_col
    )
