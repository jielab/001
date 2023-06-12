from itertools import compress
from typing import Optional, Iterable, Tuple, IO, Union, List
from multiprocessing import cpu_count, Pool
from gzip import open as gzip_open
from time import time
from tempfile import mkdtemp
import os
from collections import OrderedDict


sep_str = '\t'
NA_str = 'NA'
rs_id_str = 'RS_ID'
chrom_order_list = [str(i) for i in range(1, 22 + 1)] + ['X', 'Y', 'MT']
chrom_order = {chrom: index + 1 for index, chrom in enumerate(chrom_order_list)}
chrom_aliases = {'23': 'X', '24': 'Y', '25': 'MT', 'M': 'MT'}
columns_name = OrderedDict({'CHR': ['#CHROM', 'chr', 'CHR', 'chrom'],
                            'BP': ['BEG', 'BEGIN', 'pos', 'BP'],
                            'REF': ['reference', 'ref', 'REF'],
                            'ALT': ['alt', 'alternate', 'ALT']})


class PAGEANTERROR(Exception):
    def __init__(self, *args):
        super(PAGEANTERROR, self).__init__(*args)


def open_(file: str, mode: str = 'rb', **kwargs) -> IO:
    if file.split('.')[-1] == 'gz':
        return gzip_open(file, mode=mode, **kwargs)
    else:
        return open(file, mode=mode, **kwargs)


def rm_dir(rmdir: str) -> None:
    if os.listdir(rmdir):
        for file in os.listdir(rmdir):
            os.remove(os.path.join(rmdir, file))
    os.rmdir(rmdir)


def select_list(ob_list: list, index) -> list:
    return [ob_list[i] if i < len(ob_list) else None for i in index]


def equal_seq(seq1: str, seq2: str) -> bool:
    if seq1 == seq2:
        return True
    elif seq1 in seq2 or seq2 in seq1:
        return True
    elif len(seq1) == len(seq2) and 'N' in seq1 or 'N' in seq2:
        return all(b1 == b2 or b1 == 'N' or b2 == 'N' for b1, b2 in zip(seq1, seq2))
    else:
        return False


def get_ind(pos: List[int], pos_length: List[int]) -> List[int or None]:
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


def get_columns_index(col_names: List[str], need_cols: List[List[str]]) -> List[int or None]:
    # col_names = [name.upper().replace(' ', '_') for name in col_names]
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


class POS:
    def __init__(self, chrom: str or int, pos: str or int):
        if type(chrom) == str:
            if chrom in chrom_order:
                self.chrom = chrom_order[chrom]
            else:
                try:
                    self.chrom = chrom_order[chrom_aliases[chrom]]
                except KeyError as e:
                    raise PAGEANTERROR(f'Cannot recognize chromosome "{e.args[0]}"!') from e
        else:
            self.chrom = chrom
        if type(chrom) == str:
            try:
                self.pos = int(pos)
            except ValueError as e:
                raise PAGEANTERROR(f'Cannot recognize position{e.args[0].split(":")[-1]}') from e
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
    def __init__(self, pos: Union[POS, Tuple[str or int, str or int]]):
        if type(pos) == POS:
            pos = pos.tuple()
        super(POS_Idx, self).__init__(*pos)
        self.pos //= 100000


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


class Variant:
    def __init__(self, information: List[str] or Tuple[str, str, str, str]):
        self.pos = POS(*information[:2])
        self.ref = information[2]
        self.alt = information[3]

    def __eq__(self, other):
        if self.pos == other.pos:
            if equal_seq(self.ref, other.ref):
                return True
            elif equal_seq(self.alt, other.ref):
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
    def __init__(self, line: List[str]):
        rs_id = line.pop(2)
        super(RSID_LINE, self).__init__(line)
        self.id = rs_id


def get_line_start(pos: POS_Idx, db_file: str):
    with open(db_file + '.idx') as f:
        for line in f:
            if line.startswith(pos.export()):
                return int(line.strip().split(sep_str)[-1])


def get_rsid_part(raw_file: str, start: int, end: int, idx: int, db_file: str, temp_dir: str,
                  columns_idx: List[int], sep_str: str) -> None:
    with open(os.path.join(temp_dir, str(idx)), 'w', encoding='UTF-8') as f_write_part:
        rsid_group_reader = None
        for _, line in read_file_part(raw_file, start, end):
            line = line.decode().strip().split(sep_str)
            process = LINE(line, columns_idx)
            if not rsid_group_reader:
                rsid_group_reader = read_file_part(db_file, get_line_start(POS_Idx(process.variant.pos),
                                                                           db_file))
                _, rsid_line = next(rsid_group_reader)
                rsid_group = RSID_LINE(rsid_line.decode().strip().split('\t'))
            while rsid_group < process:
                try:
                    _, rsid_line = next(rsid_group_reader)
                    rsid_group = RSID_LINE(rsid_line.decode().strip().split('\t'))
                except StopIteration:
                    break
            while rsid_group.variant.pos == process.variant.pos:
                if rsid_group == process:
                    process.id = rsid_group.id
                    break
                else:
                    try:
                        _, rsid_line = next(rsid_group_reader)
                        rsid_group = RSID_LINE(rsid_line.decode().strip().split('\t'))
                    except StopIteration:
                        break
            f_write_part.write(process.export(sep=sep_str) + '\n')


def get_rsid(raw_file: str, db_file: str, output_file: str, processes: Optional[int] = None):
    time0 = time()
    if not processes:
        processes = cpu_count()
    print("Getting index of rs_id file...")
    if not os.path.isfile(db_file + '.idx'):
        get_index(db_file)
    print('Finish.')
    temp_dir = mkdtemp(suffix='rsid')
    with open_(raw_file) as f_raw,\
            open_(output_file, 'wb') as f_write:
        header = f_raw.readline()
        f_write.write(header.strip() + sep_str.encode() + rs_id_str.encode() + b'\n')
        header = header.decode().strip().split(sep_str)
        columns_idx = get_columns_index(header, list(columns_name.values()))
        assert None not in columns_idx,\
            f"Cannot find column " \
            f"{list(compress(['CHR', 'POS', 'REF', 'ALT'], map(lambda a: a is None, columns_idx)))} in file!"
    print("Start annotation...")
    with Pool(processes=processes) as pool:
        idx = 0
        for start, end in file_div(raw_file, processes, 1):
            idx += 1
            pool.apply_async(get_rsid_part, args=(raw_file, start, end, idx, db_file, temp_dir, columns_idx, sep_str))
        pool.close()
        pool.join()
    with open_(output_file, 'ab') as f_write:
        for i in range(1, idx + 1):
            with open(os.path.join(temp_dir, str(i)), 'rb') as f_write_part:
                for line in f_write_part:
                    f_write.write(line)
    rm_dir(temp_dir)
    print(f'Finish! Used time: {time() - time0:.0f}s')


def addrsid_run(raw_file: str, db_file: str, output_file: str, processes: Optional[int] = None,
                sep: str = '\t', chr: str = 'CHR', pos: str = 'BP', ref: str = 'REF', alt: str = 'ALT',
                rs_id: str = 'RS_ID', na: str = 'NA') -> None:
    global columns_name, sep_str, NA_str, rs_id_str
    columns_name['CHR'].append(chr)
    columns_name['BP'].append(pos)
    columns_name['REF'].append(ref)
    columns_name['ALT'].append(alt)
    sep_str = eval(repr(sep).replace('\\\\', '\\'))
    NA_str = na
    rs_id_str = rs_id
    get_rsid(raw_file, db_file, output_file, processes)
