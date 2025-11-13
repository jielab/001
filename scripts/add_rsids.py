#!/usr/bin/env python3
import argparse, gzip
from typing import IO, Dict, Tuple, Set, Optional, List, DefaultDict
from collections import defaultdict

# ------------ basic utils ------------
CHROM_ALIASES = {'23': 'X', '24': 'Y', '25': 'MT', 'M': 'MT'}
NA_STR = 'NA'

def open_any(path: str, mode: str = 'rt', encoding: Optional[str] = 'utf-8') -> IO:
    """Open plain or .gz file in text or binary mode."""
    if path.endswith('.gz'):
        if 'b' in mode:
            return gzip.open(path, mode)
        return gzip.open(path, mode, encoding=encoding)
    else:
        return open(path, mode, encoding=None if 'b' in mode else encoding)

def norm_sep(s: str) -> str:
    """Allow bash-like $'\\t' escaping in --sep."""
    return eval(repr(s).replace('\\\\', '\\'))

def detect_sep(first_line: str) -> str:
    """Auto-detect dbSNP delimiter (tab, comma, else whitespace)."""
    if '\t' in first_line: return '\t'
    if ',' in first_line: return ','
    return ' '  # whitespace

def split_with_sep(line: str, sep: str) -> List[str]:
    if sep == ' ':
        return line.rstrip('\n').split()
    return line.rstrip('\n').split(sep)

def norm_chr(ch: str) -> str:
    """Normalize chromosome notation: strip 'chr', map 23/24/25/M."""
    if ch.lower().startswith('chr'):
        ch = ch[3:]
    ch = ch.strip()
    if ch in CHROM_ALIASES:
        ch = CHROM_ALIASES[ch]
    return ch

# ------------ input & db key collection ------------
def find_col_index(header_cols: List[str], target_name: Optional[str]) -> Optional[int]:
    """Return exact index of target_name in header (case sensitive as given)."""
    if not target_name:
        return None
    try:
        return header_cols.index(target_name)
    except ValueError:
        return None

def collect_rsids(in_path: str, sep: str, idx_snp: Optional[int]) -> Set[str]:
    rsids: Set[str] = set()
    if idx_snp is None:
        return rsids
    with open_any(in_path, 'rt') as f:
        _ = f.readline()
        for line in f:
            if not line.strip(): continue
            arr = split_with_sep(line, sep)
            if idx_snp < len(arr):
                v = arr[idx_snp]
                if v and v.upper() != NA_STR:
                    rsids.add(v)
    return rsids

def collect_poskeys(in_path: str, sep: str, idx_chr: Optional[int], idx_pos: Optional[int]) -> Set[Tuple[str,int]]:
    pos_keys: Set[Tuple[str,int]] = set()
    if idx_chr is None or idx_pos is None:
        return pos_keys
    with open_any(in_path, 'rt') as f:
        _ = f.readline()
        for line in f:
            if not line.strip(): continue
            arr = split_with_sep(line, sep)
            if idx_chr < len(arr) and idx_pos < len(arr):
                ch = arr[idx_chr]; po = arr[idx_pos]
                if ch and po and po.upper() != NA_STR:
                    chn = norm_chr(ch)
                    try:
                        p = int(po)
                    except:
                        continue
                    pos_keys.add((chn, p))
    return pos_keys

# ------------ dbSNP parsing (headerless or header) ------------
def guess_db_indices_from_first_row(first_row: str, db_sep: str) -> Dict[str,int]:
    """
    Headerless db: assume CHR POS SNP [REF] [ALT]
    """
    cols = split_with_sep(first_row, db_sep)
    idx: Dict[str,int] = {}
    if len(cols) >= 3:
        idx['CHR'] = 0
        idx['POS'] = 1
        idx['SNP'] = 2
        if len(cols) >= 4:
            idx['REF'] = 3
        if len(cols) >= 5:
            idx['ALT'] = 4
    return idx

def parse_db_header(header_line: str, db_sep: str) -> Dict[str,int]:
    """
    Very permissive: accept any casing; match by equality on upper-cased tokens.
    Expected keys: CHR, POS, SNP; optional REF, ALT.
    """
    cols = [c.strip() for c in split_with_sep(header_line, db_sep)]
    up = [c.upper() for c in cols]
    idx: Dict[str,int] = {}
    for key in ['CHR','POS','SNP','REF','ALT']:
        if key in up:
            idx[key] = up.index(key)
    return idx

def scan_db_for_matches(db_path: str, rsids: Set[str], pos_keys: Set[Tuple[str,int]]
) -> Tuple[
    DefaultDict[str, List[Tuple[str,int,str,str]]],
    DefaultDict[Tuple[str,int], List[Tuple[str,str,str]]]
]:
    """
    Build:
      rsid2pos_multi: rsid -> [(chr,pos,ref,alt), ...]
      pos2rsid_multi: (chr,pos) -> [(snp,ref,alt), ...]
    Only store entries for keys present in input to bound memory.
    """
    rsid2pos_multi: DefaultDict[str, List[Tuple[str,int,str,str]]] = defaultdict(list)
    pos2rsid_multi: DefaultDict[Tuple[str,int], List[Tuple[str,str,str]]] = defaultdict(list)

    with open_any(db_path, 'rt') as f:
        first = f.readline().rstrip('\n')
        db_sep = detect_sep(first)

        # Try headered parse first
        idx = parse_db_header(first, db_sep)
        feed_first_data = False
        if not (('CHR' in idx) and ('POS' in idx) and ('SNP' in idx)):
            # treat first line as data; headerless order CHR POS SNP [REF] [ALT]
            idx = guess_db_indices_from_first_row(first, db_sep)
            if ('CHR' in idx) and ('POS' in idx) and ('SNP' in idx):
                feed_first_data = True
            else:
                tokens = split_with_sep(first, db_sep)
                raise ValueError(
                    "dbSNP columns cannot be recognized (neither header nor headerless CHR POS SNP ...). "
                    f"Detected sep: {repr(db_sep)}; first line tokens: {tokens}"
                )

        iCHR = idx['CHR']; iPOS = idx['POS']; iSNP = idx['SNP']
        iREF = idx.get('REF', None); iALT = idx.get('ALT', None)

        def handle_row(arr: List[str]):
            if max(iCHR, iPOS, iSNP) >= len(arr): return
            ch = norm_chr(arr[iCHR])
            try:
                po = int(arr[iPOS])
            except:
                return
            snp = arr[iSNP]
            ref = (arr[iREF] if iREF is not None and iREF < len(arr) else NA_STR)
            alt = (arr[iALT] if iALT is not None and iALT < len(arr) else NA_STR)

            if rsids and snp in rsids:
                rsid2pos_multi[snp].append((ch, po, ref, alt))
            if pos_keys:
                key = (ch, po)
                if key in pos_keys:
                    pos2rsid_multi[key].append((snp, ref, alt))

        if feed_first_data:
            handle_row(split_with_sep(first, db_sep))
        for line in f:
            if not line.strip(): continue
            arr = split_with_sep(line, db_sep)
            handle_row(arr)

    return rsid2pos_multi, pos2rsid_multi

# ------------ allele logic ------------
def parse_db_alt_to_set(db_ref: str, db_alt_multi: str) -> Set[str]:
    s = {db_ref.upper()} if db_ref else set()
    if db_alt_multi:
        s |= {x.strip().upper() for x in db_alt_multi.split(',') if x.strip() != ''}
    return s

def alleles_match_subset(input_ref: str, input_alt: str, db_ref: str, db_alt_multi: str) -> bool:
    """
    True iff BOTH input alleles are present among db allele choices.
    Unordered, case-insensitive. Supports multi-allelic ALT like 'A,G'.
    """
    if not input_ref or not input_alt:
        return False
    inp = {input_ref.upper(), input_alt.upper()}
    dbset = parse_db_alt_to_set(db_ref, db_alt_multi)
    return inp.issubset(dbset)

# ------------ writer ------------
def write_output(
    in_path: str, out_path: str, sep: str,
    # exact input column indices for user's names
    idx_snp: Optional[int], idx_chr: Optional[int], idx_pos: Optional[int],
    idx_ref: Optional[int], idx_alt: Optional[int],
    # db maps
    rsid2pos_multi: DefaultDict[str, List[Tuple[str,int,str,str]]],
    pos2rsid_multi: DefaultDict[Tuple[str,int], List[Tuple[str,str,str]]],
    # desired output names (exactly as user specified)
    out_snp_name: str, out_chr_name: str, out_pos_name: str,
    # whether we need to add these columns
    need_add_snp: bool, need_add_chr: bool, need_add_pos: bool,
    # enforce allele if both provided and present
    enforce_allele: bool
) -> None:

    with open_any(in_path, 'rt') as fin, open_any(out_path, 'wt') as fout:
        header = split_with_sep(next(fin), sep)
        new_header = header[:]

        # Ensure arg-named columns exist in OUTPUT header; append if missing
        if need_add_snp and out_snp_name not in header:
            new_header.append(out_snp_name)
        if need_add_chr and out_chr_name not in header:
            new_header.append(out_chr_name)
        if need_add_pos and out_pos_name not in header:
            new_header.append(out_pos_name)

        fout.write(sep.join(new_header) + '\n')

        for line in fin:
            if not line.strip(): continue
            arr = split_with_sep(line, sep)
            out_arr = arr[:]

            # cache lookups per row to avoid double work
            derived_chr = None
            derived_pos = None
            derived_snp = None

            # If we need CHR/POS and they don't exist as user-named columns, try SNP→POS
            if (need_add_chr or need_add_pos) and (out_chr_name not in header or out_pos_name not in header):
                if idx_snp is not None and idx_snp < len(arr):
                    snpval = arr[idx_snp]
                    if snpval and snpval.upper() != NA_STR:
                        candidates = rsid2pos_multi.get(snpval, [])
                        if candidates:
                            if enforce_allele and idx_ref is not None and idx_alt is not None and idx_ref < len(arr) and idx_alt < len(arr):
                                a_ref = arr[idx_ref]; a_alt = arr[idx_alt]
                                for (c_ch, c_po, c_ref, c_alt) in candidates:
                                    if alleles_match_subset(a_ref, a_alt, c_ref, c_alt):
                                        derived_chr, derived_pos = c_ch, c_po
                                        break
                            if derived_chr is None or derived_pos is None:
                                # fall back to first candidate if no allele enforcement or no match
                                c_ch, c_po, _, _ = candidates[0]
                                derived_chr, derived_pos = c_ch, c_po

            # If we need SNP and it doesn't exist as user-named column, try POS→SNP
            if need_add_snp and out_snp_name not in header:
                if idx_chr is not None and idx_pos is not None and idx_chr < len(arr) and idx_pos < len(arr):
                    ch_val = arr[idx_chr]; po_val = arr[idx_pos]
                else:
                    ch_val = None; po_val = None

                # If we computed derived CHR/POS above, use them
                if ch_val is None or po_val is None or not ch_val or not po_val or po_val.upper() == NA_STR:
                    if derived_chr is not None and derived_pos is not None:
                        ch_val, po_val = derived_chr, str(derived_pos)

                snp_out = NA_STR
                if ch_val and po_val and po_val.upper() != NA_STR:
                    chn = norm_chr(ch_val)
                    try:
                        p = int(po_val)
                        candidates = pos2rsid_multi.get((chn, p), [])
                        if candidates:
                            if enforce_allele and idx_ref is not None and idx_alt is not None and idx_ref < len(arr) and idx_alt < len(arr):
                                a_ref = arr[idx_ref]; a_alt = arr[idx_alt]
                                for (c_snp, c_ref, c_alt) in candidates:
                                    if alleles_match_subset(a_ref, a_alt, c_ref, c_alt):
                                        snp_out = c_snp
                                        break
                            if snp_out == NA_STR:
                                snp_out = candidates[0][0]
                    except:
                        pass
                derived_snp = snp_out

            # Append values for the columns we added (in the same order)
            if need_add_snp and out_snp_name not in header:
                out_arr.append(derived_snp if derived_snp is not None else NA_STR)
            if need_add_chr and out_chr_name not in header:
                vchr = NA_STR
                if idx_chr is not None and idx_chr < len(arr) and arr[idx_chr]:
                    vchr = arr[idx_chr]
                elif derived_chr is not None:
                    vchr = derived_chr
                out_arr.append(str(vchr))
            if need_add_pos and out_pos_name not in header:
                vpos = NA_STR
                if idx_pos is not None and idx_pos < len(arr) and arr[idx_pos]:
                    vpos = arr[idx_pos]
                elif derived_pos is not None:
                    vpos = derived_pos
                out_arr.append(str(vpos))

            fout.write(sep.join(out_arr) + '\n')

# ------------ main ------------
def main():
    parser = argparse.ArgumentParser(
        prog='snp_chrpos',
        description=('Ensure OUTPUT has columns named exactly as --snp/--chr/--pos. '
                     'If any are missing, add them by mapping via dbSNP (SNP↔CHR,POS). '
                     'Alleles are enforced only when both --ref and --alt are provided and present in input.')
    )
    parser.add_argument('-i', '--input', required=True, help='input table (gzip ok); input must HAVE a header')
    parser.add_argument('-d', '--database', required=True, help='dbSNP (gzip ok). Headerless allowed: CHR POS SNP [REF] [ALT]')
    parser.add_argument('-o', '--output', required=True, help='output file (gzip ok)')
    parser.add_argument('-s', '--sep', default='\t', help="delimiter for input/output (supports $'\\t')")

    # Exact OUTPUT column names to ensure
    parser.add_argument('--snp', default='SNP', help='exact SNP column name to ensure in OUTPUT')
    parser.add_argument('--chr', dest='chr_col', default='CHR', help='exact CHR column name to ensure in OUTPUT')
    parser.add_argument('--pos', dest='pos_col', default='POS', help='exact POS column name to ensure in OUTPUT')

    # Optional allele columns for enforcement (must be present in input under these exact names)
    parser.add_argument('--ref', dest='ref_col', default=None, help='exact REF column name (optional; if provided, enforce allele match)')
    parser.add_argument('--alt', dest='alt_col', default=None, help='exact ALT column name (optional; if provided, enforce allele match)')

    parser.add_argument('--na', default='NA', help="NA string (default 'NA')")

    args = parser.parse_args()

    global NA_STR
    NA_STR = args.na
    sep = norm_sep(args.sep)

    # Read input header
    with open_any(args.input, 'rt') as f:
        header_line = f.readline().rstrip('\n')
    header_cols = [c.strip() for c in split_with_sep(header_line, sep)]

    # Exact indices for user-specified names (only these matter)
    idx_snp = find_col_index(header_cols, args.snp)
    idx_chr = find_col_index(header_cols, args.chr_col)
    idx_pos = find_col_index(header_cols, args.pos_col)
    idx_ref = find_col_index(header_cols, args.ref_col) if args.ref_col else None
    idx_alt = find_col_index(header_cols, args.alt_col) if args.alt_col else None

    # Which arg-named columns are already present?
    have_arg_snp = (idx_snp is not None)
    have_arg_chr = (idx_chr is not None)
    have_arg_pos = (idx_pos is not None)

    # We must ensure OUTPUT has exactly these names; add any that are missing
    need_add_snp = not have_arg_snp
    need_add_chr = not have_arg_chr
    need_add_pos = not have_arg_pos

    # If all three are missing, we still can add them only if at least one side exists to derive:
    # - If SNP exists (under your exact name), can derive CHR/POS.
    # - If CHR & POS exist (under your exact names), can derive SNP.
    # Since we only look at exact names (no aliases), we compute keys accordingly.
    rsids: Set[str] = set()
    pos_keys: Set[Tuple[str,int]] = set()

    # If we might derive CHR/POS from SNP for new columns, collect rsids
    if (need_add_chr or need_add_pos) and have_arg_snp:
        rsids = collect_rsids(args.input, sep, idx_snp)

    # If we might derive SNP from CHR+POS for new column, collect (chr,pos)
    if need_add_snp and have_arg_chr and have_arg_pos:
        pos_keys = collect_poskeys(args.input, sep, idx_chr, idx_pos)

    # If nothing to add, just copy input to output as-is
    if not (need_add_snp or need_add_chr or need_add_pos):
        with open_any(args.input, 'rt') as fin, open_any(args.output, 'wt') as fout:
            for line in fin:
                fout.write(line if sep == '\t' else line.replace('\t', sep))
        return

    # Build db maps only if needed
    if rsids or pos_keys:
        rsid2pos_multi, pos2rsid_multi = scan_db_for_matches(args.database, rsids, pos_keys)
    else:
        rsid2pos_multi = defaultdict(list)
        pos2rsid_multi = defaultdict(list)

    # Enforce allele only if both --ref and --alt are provided AND those columns exist
    enforce_allele = (idx_ref is not None) and (idx_alt is not None) and (args.ref_col is not None) and (args.alt_col is not None)

    # Write output
    write_output(
        args.input, args.output, sep,
        idx_snp, idx_chr, idx_pos, idx_ref, idx_alt,
        rsid2pos_multi, pos2rsid_multi,
        out_snp_name=args.snp, out_chr_name=args.chr_col, out_pos_name=args.pos_col,
        need_add_snp=need_add_snp, need_add_chr=need_add_chr, need_add_pos=need_add_pos,
        enforce_allele=enforce_allele
    )

if __name__ == '__main__':
    main()
