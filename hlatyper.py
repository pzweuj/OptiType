import pandas as pd
import numpy as np
import re
import pylab
import warnings
from collections import OrderedDict
from datetime import datetime
import sys

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

VERBOSE = False

def now(start=datetime.now()):
    return str(datetime.now() - start)[:-4]

def memoize(f):
    class MemoDict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret
    return MemoDict().__getitem__

CIGAR_SLICER = re.compile(r'[0-9]+[MD]')

@memoize
def length_on_reference(cigar_string):
    return sum(int(p[:-1]) for p in re.findall(CIGAR_SLICER, cigar_string))

def feature_order(feature):
    feat_type, feat_number = feature
    assert feat_type in ('intron', 'exon', 'UTR'), 'Unknown feature in list (function accepts intron/exon/UTR)'
    assert isinstance(feat_number, int), 'Feature number has to be integer'
    return 0 if feature == ('UTR', 1) else 999 if feature == ('UTR', 2) else (feat_number * 2 + (feat_type == 'intron'))

def store_dataframes(out_hdf, **kwargs):
    complevel = kwargs.pop('complevel', 9)
    complib = kwargs.pop('complib', 'zlib')

    if VERBOSE:
        print(f"{now()} Storing {len(kwargs)} DataFrames in file {out_hdf} with compression settings {complevel} {complib}...")

    store = pd.HDFStore(out_hdf, complevel=complevel, complib=complib)
    for table_name, dataframe in kwargs.items():
        store[table_name] = dataframe
    store.close()

    if VERBOSE:
        print(f"{now()} DataFrames stored in file.")

def load_hdf(in_hdf, as_dict=False, *args):
    store = pd.HDFStore(in_hdf, 'r')
    if len(args):
        if as_dict:
            to_return = {table: store[table] for table in args}
        else:
            to_return = tuple(store[table] for table in args)
        store.close()
        return to_return
    else:
        return store

def sam_to_hdf(samfile):
    if VERBOSE:
        print(f"{now()} Loading alleles and read IDs from {samfile}...")

    read_ids, allele_ids = [], []
    first_hit_row = True
    total_hits = 0

    with open(samfile, 'rb') as f:
        last_read_id = None
        for line in f:
            line = line.decode('utf-8')
            if line.startswith('@'):
                if line.startswith('@SQ'):
                    allele_ids.append(line.split('\t')[1][3:])  # SN:HLA:HLA00001
                continue

            total_hits += 1
            read_id = line.split('\t')[0]
            if last_read_id != read_id:
                read_ids.append(read_id)
                last_read_id = read_id

            if first_hit_row:
                first_hit_row = False
                columns = line.split()
                try:
                    nm_index = list(map(lambda x: x.startswith('NM:'), columns)).index(True)
                except ValueError:
                    print('\tNo NM-tag found in SAM file!')
                    nm_index = None

    if VERBOSE:
        print(f"{now()} {len(allele_ids)} alleles and {len(read_ids)} reads found.")
        print(f"{now()} Initializing mapping matrix...")

    matrix_pos = pd.DataFrame(np.zeros((len(read_ids), len(allele_ids)), dtype=np.uint16), columns=allele_ids, index=read_ids)
    read_details = OrderedDict()

    if VERBOSE:
        print(f"{now()} {len(read_ids)}x{len(allele_ids)} mapping matrix initialized. Populating {total_hits} hits from SAM file...")

    milestones = [x * total_hits / 10 for x in range(1, 11)]

    with open(samfile, 'rb') as f:
        counter = 0
        percent = 0
        for line in f:
            line = line.decode('utf-8')
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            read_id, allele_id, position, cigar, nm = (fields[i] for i in (0, 2, 3, 5, nm_index))

            if read_id not in read_details:
                read_details[read_id] = (int(nm[5:]), length_on_reference(cigar))

            matrix_pos[allele_id][read_id] = int(position)

            counter += 1
            if counter in milestones:
                percent += 10
                if VERBOSE:
                    print(f'\t{percent}% completed')

    if VERBOSE:
        print(f"{now()} {counter} elements filled. Matrix sparsity: 1 in {matrix_pos.shape[0] * matrix_pos.shape[1] / float(counter):.2f}")

    matrix_pos.rename(columns=lambda x: x.replace('HLA:', ''), inplace=True)

    details_df = pd.DataFrame.from_dict(read_details, orient='index')
    details_df.columns = ['mismatches', 'read_length']

    return matrix_pos, details_df

def pysam_to_hdf(samfile):
    if not PYSAM_AVAILABLE:
        print("Warning: PySam not available on the system. Falling back to primitive SAM parsing.")
        return sam_to_hdf(samfile)

    sam_or_bam = 'rb' if samfile.endswith('.bam') else 'r'
    sam = pysam.AlignmentFile(samfile, sam_or_bam)
    is_yara = (sam.header['PG'][0]['ID'] in ('Yara', 'yara'))
    xa_tag = is_yara and (' -os ' not in sam.header['PG'][0]['CL'])

    nref = sam.nreferences
    hits = OrderedDict()
    allele_id_to_index = {aa: ii for ii, aa in enumerate(sam.references)}
    read_details = OrderedDict()

    if VERBOSE:
        print(f"{now()} Loading {samfile} started. Number of HLA reads loaded (updated every thousand):")

    read_counter = 0
    hit_counter = 0

    for aln in sam:
        if aln.qname not in hits:
            hits[aln.qname] = np.zeros(nref, dtype=np.uint16)
            read_details[aln.qname] = (aln.get_tag('NM'), aln.query_length)
            read_counter += 1
            if VERBOSE and not (read_counter % 1000):
                sys.stdout.write(f'{len(hits) // 1000}K...')
                sys.stdout.flush()
            if xa_tag and aln.has_tag('XA'):
                current_row = hits[aln.qname]
                subtags = aln.get_tag('XA').split(';')[:-1]
                hit_counter += len(subtags)
                for subtag in subtags:
                    allele, pos = subtag.split(',')[:2]
                    current_row[allele_id_to_index[allele]] = int(pos)

        hits[aln.qname][aln.reference_id] = aln.reference_start + 1
        hit_counter += 1

    if VERBOSE:
        print(f"\n{now()} {len(hits)} reads loaded. Creating dataframe...")
    pos_df = pd.DataFrame.from_dict(hits, orient='index')
    pos_df.columns = sam.references[:]
    details_df = pd.DataFrame.from_dict(read_details, orient='index')
    details_df.columns = ['mismatches', 'read_length']
    if VERBOSE:
        print(f"{now()} Dataframes created. Shape: {pos_df.shape[0]} x {pos_df.shape[1]}, hits: {np.sign(pos_df).sum().sum()} ({hit_counter}), sparsity: 1 in {pos_df.shape[0] * pos_df.shape[1] / float(hit_counter):.2f}")
    return pos_df, details_df

def get_compact_model(hit_df, weak_hit_df=None, weight=None):
    hit_df = hit_df.loc[hit_df.any(axis=1)]
    occurence = {r[0]: len(r) for r in hit_df.groupby(hit_df.columns.tolist()).groups.values()}

    if weak_hit_df is not None:
        weak_hit_df = weak_hit_df.loc[weak_hit_df.any(axis=1)]
        assert 0 < weight <= 1, 'weak hit weight must be in (0, 1]'
        weak_occ = {r[0]: len(r) * weight for r in weak_hit_df.groupby(weak_hit_df.columns.tolist()).groups.values()}
        occurence.update(weak_occ)
        unique_mtx = pd.concat([hit_df.drop_duplicates(), weak_hit_df.drop_duplicates()])
    else:
        unique_mtx = hit_df.drop_duplicates()
    
    return unique_mtx, occurence

def mtx_to_sparse_dict(hit_df):
    all_hits = {}
    for read_id, alleles in hit_df.iterrows():
        hit_alleles = alleles[alleles != 0].index
        for hit_allele in hit_alleles:
            all_hits[(read_id, hit_allele)] = 1
    return all_hits

def create_allele_dataframes(imgt_dat, fasta_gen, fasta_nuc):
    from Bio import SeqIO
    if VERBOSE:
        print(f"{now()} Loading IMGT allele dat file...")

    alleles = OrderedDict()

    with open(imgt_dat, 'r') as handle:
        for i, record in enumerate(SeqIO.parse(handle, "imgt")):
            record.id = record.id.split('.')[0]
            alleles[record.id] = record

    if VERBOSE:
        print(f"{now()} Initializing allele DataFrame...")

    allele_info = 'id type 4digit locus flags len_dat len_gen len_nuc full_gen full_nuc'
    table = pd.DataFrame(index=list(alleles.keys()), columns=allele_info.split())
    sequences = []

    if VERBOSE:
        print(f"{now()} Filling DataFrame with allele data...")

    all_features = []

    for allele in alleles.values():
        allele_type = allele.description.replace('HLA-', '').split(',')[0]
        table.loc[allele.id]['id'] = allele.id
        table.loc[allele.id]['type'] = allele_type
        table.loc[allele.id]['4digit'] = ':'.join(allele_type.split(':')[:2])
        table.loc[allele.id]['locus'] = allele_type.split('*')[0]
        table.loc[allele.id]['flags'] = 0 if allele_type[-1].isdigit() else 1
        table.loc[allele.id]['len_dat'] = len(str(allele.seq))
        table.loc[allele.id]['len_gen'] = 0
        table.loc[allele.id]['len_nuc'] = 0
        table.loc[allele.id]['full_gen'] = 0
        table.loc[allele.id]['full_nuc'] = 0
        sequences.append((allele.id, 'dat', str(allele.seq)))

        features = [f for f in allele.features if f.type in ('exon', 'intron', 'UTR')]
        for feature in features:
            if feature.type in ('exon', 'intron'):
                feature_num = int(feature.qualifiers['number'][0])
            else:
                assert feature.location.start == 0 or feature.location.end == len(allele.seq)
                feature_num = 1 if feature.location.start == 0 else 2

            all_features.append(
                (allele.id, feature.type, feature_num,
                int(feature.location.start), int(feature.location.end), len(feature), feature_order((feature.type, feature_num)))
            )

        cds = [f for f in allele.features if f.type == 'CDS']
        if cds:
            if sum(map(len, [f for f in features if f.type == 'exon'])) != len(cds[0]):
                if VERBOSE:
                    print("\tCDS length doesn't match sum of exons for", allele.id, allele_type)
                table.loc[allele.id]['flags'] += 2
        else:
            if VERBOSE:
                print("\tNo CDS found for", allele.id, allele_type)
            table.loc[allele.id]['flags'] += 2

    if VERBOSE:
        print(f"{now()} Loading gen and nuc files...")

    with open(fasta_gen, 'r') as fasta_gen:
        for record in SeqIO.parse(fasta_gen, 'fasta'):
            allele_id = record.id.replace('HLA:', '')
            table.loc[allele_id]['len_gen'] = len(record.seq)
            sequences.append((allele_id, 'gen', str(record.seq)))

    with open(fasta_nuc, 'r') as fasta_nuc:
        for record in SeqIO.parse(fasta_nuc, 'fasta'):
            allele_id = record.id.replace('HLA:', '')
            table.loc[allele_id]['len_nuc'] = len(record.seq)
            sequences.append((allele_id, 'nuc', str(record.seq)))

    all_features = pd.DataFrame(all_features, columns=['id', 'feature', 'number', 'start', 'end', 'length', 'order'])
    sequences = pd.DataFrame(sequences, columns=['id', 'source', 'sequence'])

    joined = pd.merge(table, all_features, how='inner', on='id')

    exons_for_locus = {}
    for i_locus, i_group in joined.groupby('locus'):
        exons_for_locus[i_locus] = i_group[i_group['feature'] == 'exon']['number'].max()

    if VERBOSE:
        print(f"{now()} Checking dat features vs gen/nuc sequences...")

    for allele, features in joined.groupby('id'):
        row = features.iloc[0]
        sum_features_length = features['length'].sum()
        sum_exons_length = features.loc[features['feature'] == 'exon']['length'].sum()
        if row['len_gen'] > 0 and row['len_gen'] != sum_features_length:
            if VERBOSE:
                print("\tFeature lengths don't add up to gen sequence length", allele, row['len_gen'], sum_features_length, row['type'])
            table.loc[allele]['flags'] += 4
        if row['len_nuc'] > 0 and row['len_nuc'] != sum_exons_length:
            if VERBOSE:
                print("\tExon lengths don't add up to nuc sequence length", allele, row['len_nuc'], sum_exons_length, row['type'])
            table.loc[allele]['flags'] += 8

    if VERBOSE:
        print(f"{now()} Sanity check finished. Computing feature sequences...")

    ft_seq_lookup = OrderedDict()
    ft_seq_lookup['---DUMMY---'] = 0
    ft_counter = 1
    all_ft_counter = 0
    all_features['seq_id'] = 0
    for i_id, i_features in all_features.groupby('id'):
        seq = sequences.loc[(sequences['id'] == i_id) & (sequences['source'] == 'dat')].iloc[0]['sequence']
        for ft_idx, feature in i_features.iterrows():
            ft_seq = seq[feature['start']:feature['end']]
            all_ft_counter += 1
            if ft_seq not in ft_seq_lookup:
                ft_seq_lookup[ft_seq] = ft_counter
                all_features.loc[ft_idx, 'seq_id'] = ft_counter
                ft_counter += 1
            else:
                all_features.loc[ft_idx, 'seq_id'] = ft_seq_lookup[ft_seq]

    feature_sequences = pd.DataFrame(list(ft_seq_lookup.keys()), columns=['sequence'])

    return table, all_features, sequences, feature_sequences

def prune_identical_alleles(binary_mtx, report_groups=False):
    hash_columns = binary_mtx.transpose().dot(np.random.rand(binary_mtx.shape[0]))
    if report_groups:
        grouper = hash_columns.groupby(hash_columns)
        groups = {g[1].index[0]: g[1].index.tolist() for g in grouper}
    alleles_to_keep = hash_columns.drop_duplicates().index
    return binary_mtx[alleles_to_keep] if not report_groups else (binary_mtx[alleles_to_keep], groups)

def prune_identical_reads(binary_mtx):
    reads_to_keep = binary_mtx.dot(np.random.rand(binary_mtx.shape[1])).drop_duplicates().index
    return binary_mtx.loc[reads_to_keep]

def prune_overshadowed_alleles(binary_mtx):
    if (binary_mtx.shape[0] < np.iinfo(np.uint16).max) or (binary_mtx.sum(axis=0).max() < np.iinfo(np.uint16).max):
        bb = binary_mtx if all(binary_mtx.dtypes == np.uint16) else binary_mtx.astype(np.uint16)
    else:
        bb = binary_mtx.astype(np.uint32)

    covariance = bb.transpose().dot(bb)
    diagonal = pd.Series([covariance[ii][ii] for ii in covariance.columns], index=covariance.columns)
    new_covariance = covariance[covariance.columns]
    for ii in new_covariance.columns:
        new_covariance[ii][ii] = 0
    overshadowed = []
    for ii in new_covariance.columns:
        potential_superiors = new_covariance[ii][new_covariance[ii] == diagonal[ii]].index
        if any(diagonal[potential_superiors] > diagonal[ii]):
            overshadowed.append(ii)
    non_overshadowed = covariance.columns.difference(overshadowed)
    return non_overshadowed

def create_paired_matrix(binary_1, binary_2, id_cleaning=None):
    if id_cleaning is not None:
        binary_1.index = list(map(id_cleaning, binary_1.index))
        binary_2.index = list(map(id_cleaning, binary_2.index))

    common_read_ids = binary_1.index.intersection(binary_2.index)
    only_1 = binary_1.index.difference(binary_2.index)
    only_2 = binary_2.index.difference(binary_1.index)

    b_1 = binary_1.loc[common_read_ids]
    b_2 = binary_2.loc[common_read_ids]
    b_12 = b_1 * b_2
    b_ispaired = b_12.any(axis=1)
    b_paired = b_12.loc[b_ispaired]
    b_mispaired = b_1.loc[~b_ispaired] + b_2.loc[~b_ispaired]
    b_unpaired = pd.concat([binary_1.loc[only_1], binary_2.loc[only_2]])

    if VERBOSE:
        print(f"{now()} Alignment pairing completed. {b_paired.shape[0]} paired, {b_unpaired.shape[0]} unpaired, {b_mispaired.shape[0]} discordant")

    return b_paired, b_mispaired, b_unpaired

def get_features(allele_id, features, feature_list):
    if '_' in allele_id:
        partial_allele, complete_allele = allele_id.split('_')
    else:
        complete_allele = allele_id

    feats_complete = {(of['feature'], of['number']): of for _, of in features.loc[features['id'] == complete_allele].iterrows()}
    feats_partial = {(of['feature'], of['number']): of for _, of in features.loc[features['id'] == partial_allele].iterrows()} if '_' in allele_id else feats_complete

    feats_to_include = []

    for feat in sorted(feature_list, key=feature_order):
        if feat in feats_partial:
            feats_to_include.append(feats_partial[feat])
        elif feat in feats_complete:
            feats_to_include.append(feats_complete[feat])
        else:
            warnings.warn(f'Feature {feat} not found for allele {allele_id}')

    return pd.DataFrame(feats_to_include)

def calculate_coverage(alignment, features, alleles_to_plot, features_used):
    assert len(alignment) in (2, 4, 5), ("Alignment tuple either has to have 2, 4 or 5 elements. First four: pos, read_details "
        "once or twice depending on single or paired end, and an optional binary DF at the end for PROPER paired-end plotting")
    has_pairing_info = (len(alignment) == 5)

    if len(alignment) == 2:
        matrix_pos, read_details = alignment[:2]
        pos = matrix_pos[alleles_to_plot]
        pairing = np.sign(pos) * 2
        hit_counts = np.sign(pos).sum(axis=1)
        max_ambiguity = max(1, hit_counts.max())
        to_process = [(pos, read_details, hit_counts, pairing)]
    elif len(alignment) == 4:
        matrix_pos1, read_details1, matrix_pos2, read_details2 = alignment
        pos1 = matrix_pos1[alleles_to_plot]
        pos2 = matrix_pos2[alleles_to_plot]
        pairing1 = np.sign(pos1) * 2
        pairing2 = np.sign(pos2)
        hit_counts1 = np.sign(pos1).sum(axis=1)
        hit_counts2 = np.sign(pos2).sum(axis=1)
        max_ambiguity = max(1, hit_counts1.max(), hit_counts2.max())
        to_process = [(pos1, read_details1, hit_counts1, pairing1), (pos2, read_details2, hit_counts2, pairing2)]
    else:
        matrix_pos1, read_details1, matrix_pos2, read_details2, pairing_binaries = alignment
        bin_p, bin_u, bin_m = pairing_binaries
        pos1 = matrix_pos1[alleles_to_plot]
        pos2 = matrix_pos2[alleles_to_plot]
        pairing = pd.concat([bin_p[alleles_to_plot], bin_u[alleles_to_plot] * 2, bin_m[alleles_to_plot] * 3])
        pairing1 = pairing.loc[pos1.index]
        pairing2 = pairing.loc[pos2.index]
        hit_counts1 = np.sign(pairing1).sum(axis=1)
        hit_counts2 = np.sign(pairing2).sum(axis=1)
        max_ambiguity = max(1, hit_counts1.max(), hit_counts2.max())
        to_process = [(pos1, read_details1, hit_counts1, pairing1), (pos2, read_details2, hit_counts2, pairing2)]

    coverage_matrices = []

    for allele in alleles_to_plot:
        allele_features = get_features(allele, features, features_used)
        allele_length = allele_features['length'].sum()

        coverage = np.zeros((2, 3, max_ambiguity, allele_length), dtype=int)

        for pos, read_details, hit_counts, pairing_info in to_process:
            reads = pos[pos[allele] != 0].index
            for i_pos, i_read_length, i_mismatches, i_hitcount, i_pairing in zip(
                    pos.loc[reads][allele],
                    read_details.loc[reads]['read_length'],
                    read_details.loc[reads]['mismatches'],
                    hit_counts[reads],
                    pairing_info.loc[reads][allele]):
                if not i_pairing:
                    continue
                coverage[int(bool(i_mismatches))][i_pairing - 1][i_hitcount - 1][i_pos - 1:i_pos - 1 + i_read_length] += 1

        coverage_matrices.append((allele, coverage))
    return coverage_matrices

def plot_coverage(outfile, coverage_matrices, allele_data, features, features_used, columns=2):
    def start_end_zeros(cov_array):
        return np.append(np.append([0], cov_array), [0])

    def allele_sorter(allele_cov_mtx_tuple):
        allele, _ = allele_cov_mtx_tuple
        return allele_data.loc[allele.split('_')[0]]['type']

    def get_allele_locus(allele):
        return allele_data.loc[allele.split('_')[0]]['locus']

    number_of_loci = len(set((get_allele_locus(allele) for allele, _ in coverage_matrices)))

    dpi = 50
    box_size = (7, 1)
    subplot_rows = 3 * number_of_loci + 1

    area_colors = [
        (0.26, 0.76, 0.26),
        (0.40, 0.84, 0.40),
        (0.99, 0.75, 0.20),
        (0.99, 0.75, 0.20),
        (0.99, 0.85, 0.35),
        (0.99, 0.85, 0.35),
        (0.99, 0.23, 0.23),
        (0.99, 0.49, 0.49),
        (0.14, 0.55, 0.72),
        (0.14, 0.55, 0.72),
        (0.33, 0.70, 0.88),
        (0.33, 0.70, 0.88)
    ]

    figure = pylab.figure(figsize=(box_size[0] * columns, box_size[1] * subplot_rows), dpi=dpi)

    coverage_matrices = sorted(coverage_matrices, key=allele_sorter)
    prev_locus = ''
    i_locus = -1

    for allele, coverage in coverage_matrices:
        if '_' in allele:
            partial, complete = allele.split('_')
            plot_title = f'{allele_data.loc[partial]["type"]} (introns from {allele_data.loc[complete]["type"]})'
        else:
            plot_title = allele_data.loc[allele]['type']

        if prev_locus != get_allele_locus(allele):
            i_locus += 1
            i_allele_in_locus = 0
        else:
            i_allele_in_locus = 1

        prev_locus = get_allele_locus(allele)

        plot = pylab.subplot2grid((subplot_rows, columns), (3 * i_locus, i_allele_in_locus), rowspan=3, adjustable='box')

        _, _, max_ambig, seq_length = coverage.shape

        shared_weighting = np.reciprocal(np.arange(max_ambig) + 1.0)
        shared_weighting[0] = 0

        perfect_paired_unique = start_end_zeros(coverage[0][0][0])
        mismatch_paired_unique = start_end_zeros(coverage[1][0][0])
        perfect_unpaired_unique = start_end_zeros(coverage[0][1][0])
        mismatch_unpaired_unique = start_end_zeros(coverage[1][1][0])
        perfect_mispaired_unique = start_end_zeros(coverage[0][2][0])
        mismatch_mispaired_unique = start_end_zeros(coverage[1][2][0])

        perfect_paired_shared = start_end_zeros(shared_weighting.dot(coverage[0][0]))
        mismatch_paired_shared = start_end_zeros(shared_weighting.dot(coverage[1][0]))
        perfect_unpaired_shared = start_end_zeros(shared_weighting.dot(coverage[0][1]))
        mismatch_unpaired_shared = start_end_zeros(shared_weighting.dot(coverage[1][1]))
        perfect_mispaired_shared = start_end_zeros(shared_weighting.dot(coverage[0][2]))
        mismatch_mispaired_shared = start_end_zeros(shared_weighting.dot(coverage[1][2]))

        i_start = 1
        for _, ft in get_features(allele, features, features_used).iterrows():
            if ft['feature'] == 'exon':
                plot.axvspan(i_start, i_start + ft['length'], facecolor='black', alpha=0.1, linewidth=0, zorder=1)
            i_start += ft['length']

        areas = plot.stackplot(np.arange(seq_length + 2),
            perfect_paired_unique + 0.001,
            perfect_paired_shared,
            perfect_unpaired_unique,
            perfect_mispaired_unique,
            perfect_unpaired_shared,
            perfect_mispaired_shared,
            mismatch_paired_unique,
            mismatch_paired_shared,
            mismatch_unpaired_unique,
            mismatch_mispaired_unique,
            mismatch_unpaired_shared,
            mismatch_mispaired_shared,
            linewidth=0, colors=area_colors, zorder=5)

        for aa in areas:
            aa.set_edgecolor(aa.get_facecolor())

        plot.tick_params(axis='both', labelsize=10, direction='out', which='both', top=False)

        plot.text(.015, 0.97, plot_title, horizontalalignment='left', verticalalignment='top', transform=plot.transAxes, fontsize=10, zorder=6)
        _, _, _, y2 = plot.axis()
        plot.axis((0, seq_length, 1, y2))
        plot.set_yscale('log')
        plot.set_ylim(bottom=0.5)

    legend = pylab.subplot2grid((subplot_rows, columns), (subplot_rows - 1, 0), colspan=2, adjustable='box')
    ppp = pylab.matplotlib.patches
    legend.add_patch(ppp.Rectangle((0, 2), 2, 2, color=area_colors[0]))
    legend.add_patch(ppp.Rectangle((0, 0), 2, 2, color=area_colors[1]))
    legend.add_patch(ppp.Rectangle((25, 2), 2, 2, color=area_colors[2]))
    legend.add_patch(ppp.Rectangle((25, 0), 2, 2, color=area_colors[4]))
    legend.add_patch(ppp.Rectangle((50, 2), 2, 2, color=area_colors[6]))
    legend.add_patch(ppp.Rectangle((50, 0), 2, 2, color=area_colors[7]))
    legend.add_patch(ppp.Rectangle((75, 2), 2, 2, color=area_colors[8]))
    legend.add_patch(ppp.Rectangle((75, 0), 2, 2, color=area_colors[10]))
    legend.text(2.5, 3, 'paired, no mismatches, unique', va='center', size='smaller')
    legend.text(2.5, 1, 'paired, no mismatches, ambiguous', va='center', size='smaller')
    legend.text(27.5, 3, 'unpaired, no mismatches, unique', va='center', size='smaller')
    legend.text(27.5, 1, 'unpaired, no mismatches, ambiguous', va='center', size='smaller')
    legend.text(52.5, 3, 'paired, mismatched, unique', va='center', size='smaller')
    legend.text(52.5, 1, 'paired, mismatched, ambiguous', va='center', size='smaller')
    legend.text(77.5, 3, 'unpaired, mismatched, unique', va='center', size='smaller')
    legend.text(77.5, 1, 'unpaired, mismatched, ambiguous', va='center', size='smaller')
    legend.set_xlim(0, 100)
    legend.set_ylim(0, 4)
    legend.axison = False

    figure.tight_layout()
    figure.savefig(outfile)
