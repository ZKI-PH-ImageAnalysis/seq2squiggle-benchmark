import pysam
import re
import sys
import numpy as np 
import os
import rich_click as click
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
from matplotlib.lines import Line2D
from pysam import FastxFile
from tqdm import tqdm
import seaborn as sns
import pandas as pd
from matplotlib.collections import PathCollection
import logging
import warnings
logger = logging.getLogger("basecalling-stats")


base2index = {"A":0, "C":1, "G":2, "T":3, "_":4}


def update_err_dist(err_dist, symbol, homo_flag):
    symbol2name = {"X": "mis", "D": "del", "I": "ins"}
    if symbol in symbol2name:
        if homo_flag:
            err_dist["ho_" + symbol2name[symbol]] += 1
        else:
            err_dist["he_" + symbol2name[symbol]] += 1

def read_samfile(samfile, subsample = 10000000, mapq_thres = 20, seqlen_thres = 100, return_ids = False):
    seqs, supp_count, ids = [], 0, set()
    failed_mapping = 0
    count_s = 0
    lowq_count = 0
    with open(samfile, "r") as f:
        for line in f:
            if len(seqs) == subsample:
                break
            if line[0] != "@":
                records = line.strip().split("\t")
                name, flag, start, mapq, cigar, seq, chromo, qscore = records[0], records[1], records[3], records[4], records[5], records[9], records[2], records[10]
                count_s += 1
                if int(flag) in set([4, 73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137, 77 , 141]):
                    failed_mapping += 1
                elif int(mapq) < mapq_thres and len(seq) < seqlen_thres:
                    lowq_count += 1
                else:
                    if return_ids == True:
                        ids.add(name)
                    else:
                        seqs.append((name, flag, start, cigar, seq, mapq, chromo, qscore))
    if return_ids == True:
        return ids
    return seqs, failed_mapping, count_s, lowq_count


def get_read_stats(cigar_splitted, strand):
    match, mismatch, insertion, deletion, n_clips = [], [], [], [], [0, 0]
    for index in range(0, len(cigar_splitted)-1, 2): # items[-1] is actually ""
        count, symbol = int(cigar_splitted[index]), cigar_splitted[index+1]
        if symbol == "=":
            match.append(int(count))
        elif symbol == "X":
            mismatch.append(int(count))
        elif symbol == "I":
            insertion.append(int(count))
        elif symbol == "D":
            deletion.append(int(count))
        elif symbol == "S":
            if strand == "+":
                if index == 0:
                    n_clips[0] += count
                else:
                    n_clips[1] += count
            if strand == "-":
                if index == 0:
                    n_clips[1] += count
                else:
                    n_clips[0] += count
    return match, mismatch, insertion, deletion, n_clips

def symbol2qscore(symbols):
    quals = np.fromiter(
            (ord(x) - 33 for x in symbols),
            dtype=int, count=len(symbols))
    mean_p = np.mean(np.power(10, quals/-10))
    return -10*np.log10(mean_p)


def get_clipped_bases(cigar_splitted):
    clip_start = int(cigar_splitted[0]) if cigar_splitted[1] == "S" else 0
    clip_end = int(cigar_splitted[-3]) if cigar_splitted[-2] == "S" else 0
    return clip_start, clip_end



def get_sample_stats(reference, sam, qscores=None, lengths=None, subsample = 10000000, return_counts = False):
    # print("Processing sample", sam)
    seqs, failed_map, total_c, lowq_c = read_samfile(sam, subsample, mapq_thres=20, seqlen_thres=100)

    read_stats, read_q, empirical_q = {}, [], []
    n_seqs = np.minimum(len(seqs), subsample)
    
    read_stats["metrics"] = np.empty((n_seqs, 3), dtype = np.float64)
    read_stats["acc"] = np.empty((n_seqs), dtype = np.float64)
    read_stats["qscores"] = np.empty((n_seqs, 3), dtype = np.float64)
    
    read_stats["align_ratio"] = np.empty(n_seqs, dtype = np.float64)
    read_stats["original_read_length"] = np.empty(n_seqs, dtype = np.int64)
    read_stats["aligned_read_length"] = np.empty(n_seqs, dtype = np.int64)
    read_stats["read_q"] = np.empty(n_seqs, dtype = np.float64)
    read_stats["empirical_q"] = np.empty(n_seqs, dtype = np.float64)
    read_stats["n_clips"] = np.empty((n_seqs, 2), dtype = np.float64)
    
    mismatches, insertions, deletions = [], [], []
    
    for seq_num in tqdm(range(n_seqs)):
        name, flag, start, cigar, sample_seq, mapq, chromo, qs_base = seqs[seq_num]
        # qscore, length = qscores[name], lengths[name]
        
        strand = '-' if int(flag) & 0x10 else '+'
        
        cigar_splitted = re.split('([^0-9])',cigar.strip())
        match, mismatch, insertion, deletion, n_clips = get_read_stats(cigar_splitted, strand = strand)
        
        if return_counts == True:
            mismatches += mismatch
            insertions += insertion
            deletions += deletion
        
        match, mismatch, insertion, deletion = np.sum(match), np.sum(mismatch), np.sum(insertion), np.sum(deletion)
        align_ratio = 1 - (n_clips[0] + n_clips[1]) / len(sample_seq)
        read_l = match + mismatch + insertion + deletion
        read_stats["original_read_length"][seq_num] = len(sample_seq)
        read_stats["metrics"][seq_num, 0] = mismatch/read_l
        read_stats["metrics"][seq_num, 1] = insertion/read_l
        read_stats["metrics"][seq_num, 2] = deletion/read_l
        read_stats["acc"][seq_num] = match / read_l
        
        read_stats["align_ratio"][seq_num] = align_ratio
        read_stats["n_clips"][seq_num, :] = n_clips


        #read_stats["read_q"][seq_num] = qscore
        read_stats["empirical_q"][seq_num] = symbol2qscore(qs_base)
        clip_start, clip_end = get_clipped_bases(cigar_splitted)
        qscore_clipped = symbol2qscore(qs_base[clip_start : -clip_end-1])
        #read_stats["qscores"][seq_num, 0] = qscore
        #read_stats["qscores"][seq_num, 1] = qscore_clipped

    
    if return_counts == True:
        return read_stats, np.array(mismatches), np.array(insertions), np.array(deletions), failed_map
    else:
        return read_stats, failed_map, total_c, lowq_c, len(seqs)


def init_kmer_dict(k = 5):
    kmer_dict, kmers = {}, ["A", "C", "G", "T"]
    alphabet = set("ACGT")
    
    for i in range(k-1):
        kmers_extended = []
        for kmer in kmers:
            for letter in alphabet:
                kmers_extended.append(kmer+letter)
        kmers = kmers_extended
        
    for kmer in kmers:
        kmer_dict[kmer] = dict()
        
    return kmer_dict

def init_homodict(min_size, max_size):
    homos = set()
    alphabet = ["A", "C", "G", "T"]
    for i in range(min_size, max_size + 1):
        for letter in alphabet:
            homos.add(letter * i)
    return homos

def get_homo_pos(seq, max_len = 9):
    homo2pos = []
    for homopolymer in re.finditer(r'([ACGT])\1{1,}', seq):
        homo_seq = homopolymer.group()
        homo_start = homopolymer.start()
        l = len(homo_seq)
        if l <= max_len:
            homo2pos.append((homo_seq, homo_start, l))
    return homo2pos

def is_homop(seq, letter):
    if "_" in seq:
        seq = "".join(seq.split("_"))
    if (len(set(seq)) == 1 and letter == seq[0]) or seq == '':
        return seq
    return False


def rev_compl(seq):
    base2comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(base2comp[base] for base in reversed(seq))



def get_homop_profiles(samfile, ref, subsample = 100000, min_l = 2, max_l = 7):
    seqs, _, _,_ = read_samfile(samfile, subsample)

    kmer_dict = {kmer:[0, 0, []] for kmer in init_homodict(min_l, max_l)}
    for seq_num in tqdm(range(np.minimum(len(seqs), subsample))):
        name, flag, start, cigar, sample_seq, mapq, chromo, qs_base = seqs[seq_num]
        cigar_splitted = re.split('([^0-9])',cigar.strip())
        cigar_expanded, n_clips = '', 0
        for index in range(0, len(cigar_splitted)-1, 2): # items[-1] is actually ""
            count, symbol = int(cigar_splitted[index]), cigar_splitted[index+1]
            if symbol == "S":
                n_clips += count
            else:
                cigar_expanded = cigar_expanded + ''.join([symbol * count])
        ref_start = int(start) - 1 # samfiles are 1-based
        ref_end = ref_start + len(sample_seq) - cigar_expanded.count("I") + cigar_expanded.count("D") + cigar_expanded.count("N") - n_clips
        
        sample_seq = sample_seq[int(cigar_splitted[0]):] if cigar_splitted[1] == "S" else sample_seq
        sample_seq = sample_seq[:-int(cigar_splitted[-3])] if cigar_splitted[-2] == "S" else sample_seq
        ref_seq = ref.fetch(reference = chromo, start = ref_start, end = ref_end)
        
        if any(base not in "ACTG" for base in ref_seq) or ref_seq.upper() != ref_seq:
            continue
            
        strand = '-' if int(flag) & 0x10 else '+'
        if strand == '-':
            sample_seq = rev_compl(sample_seq)
            ref_seq = rev_compl(ref_seq)
            cigar_expanded = cigar_expanded[::-1]
        
        # Expand sample sequence by adding deletions _ and removing insertions
        expanded_sample, p_sample, p_ref, insert_pos = '', 0, 0, {}
        p_cigar = 0 if cigar_splitted[1] != "S" else 2
        
        for i in range(len(cigar_expanded)):
            if cigar_expanded[i] == "=" or cigar_expanded[i] == "X":
                expanded_sample += sample_seq[p_sample]
                p_sample += 1
                p_ref += 1
            elif cigar_expanded[i] == "D":
                expanded_sample += "_"
                p_ref += 1
            elif cigar_expanded[i] == "N":
                expanded_sample += "N"
                p_ref += 1
            elif cigar_expanded[i] == "I":
                if cigar_expanded[i-1] != "I":
                    insert_pos[p_ref] = sample_seq[p_sample] if sample_seq[p_sample] != "T" else "T" 
                    p_sample += 1
                    step = 1
                    while cigar_expanded[i+step] == "I":
                        insert_pos[p_ref] += sample_seq[p_sample] if sample_seq[p_sample] != "T" else "T"
                        step += 1
                        p_sample += 1                        
        
        if len(expanded_sample) != len(ref_seq):
            print("expansion didn't work!!")
        
        for kmer, homo_start, l in get_homo_pos(ref_seq, max_l):
            if "N" not in kmer:
                sample_kmer = expanded_sample[homo_start:homo_start + l]
                ins_flag, base2insert = False, []
                for j in range(homo_start, homo_start+l):
                    if j in insert_pos:
                        ins_flag = True
                        base2insert.append((j-i, insert_pos[j]))
                if ins_flag == True:
                    for loc, insert in base2insert:
                        inserted_kmer = sample_kmer[:loc] + insert + sample_kmer[loc:]
                    sample_kmer = inserted_kmer
                    # print(i, kmer, sample_kmer, inserted_kmer, base2insert)

                kmer_dict[kmer][1] += 1
                if sample_kmer == kmer:
                    kmer_dict[kmer][0] += 1
                sample_kmer = is_homop(sample_kmer, kmer[0]) 
                if sample_kmer != False:
                    kmer_dict[kmer][2].append(len(sample_kmer))
    return kmer_dict
    

def get_kmer_profiles(samfile, ref, subsample = 100000, kmer_size = 5):
    # print("Processing sample", samfile)
    seqs, _, _, _ = read_samfile(samfile, subsample)
    
    kmer_dict, confusion = init_kmer_dict(kmer_size), np.zeros((5, 4), dtype = int)
    motif_acc = {kmer: {} for kmer in ['GA', 'AG', 'GT', 'TG', 'GC', 'CG', 'AC', 'CA', 'AT','TA', 'TC', 'CT',]}
    err_dist = {"ho_mis": 0, "ho_del": 0, "ho_ins": 0, "he_mis": 0, "he_del": 0, "he_ins": 0}
    
    for seq_num in tqdm(range(np.minimum(len(seqs), subsample))):
        name, flag, start, cigar, sample_seq, mapq, chromo, qs_base = seqs[seq_num]
        cigar_splitted, cigar_expanded, n_clips = re.split('([^0-9])', cigar.strip()), '', 0
        for index in range(0, len(cigar_splitted)-1, 2): # items[-1] is ""
            count, symbol = int(cigar_splitted[index]), cigar_splitted[index+1]
            if symbol == "S": # Is hard-clipping H dealt with?
                n_clips += count
            else:
                cigar_expanded = cigar_expanded + ''.join([symbol * count])
                
        ref_start = int(start) - 1 # samfiles are 1-based
        ref_end = ref_start + len(sample_seq) - cigar_expanded.count("I") + cigar_expanded.count("D") + cigar_expanded.count("N") - n_clips
        ref_seq = ref.fetch(reference = chromo, start = ref_start, end = ref_end)
        if "N" in ref_seq or ref_seq.upper() != ref_seq: # Filtering non-standard reference sequences?
            continue
        
        # to clip the soft-clipped bases
        sample_seq = sample_seq[int(cigar_splitted[0]):] if cigar_splitted[1] == "S" else sample_seq
        sample_seq = sample_seq[:-int(cigar_splitted[-3])] if cigar_splitted[-2] == "S" else sample_seq
            
        if int(flag) & 0x10 == '-':
            sample_seq, ref_seq, cigar_expanded = rev_compl(sample_seq), rev_compl(ref_seq), cigar_expanded[::-1]
        
        # to find the positions of homopolymers on the read, from left -1 to right + 1
        homo_pos = set()
        for kmer, homo_start, l in get_homo_pos(ref_seq, 1000000):
            for increment in range(-1, l+1, 1):
                homo_pos.add(homo_start + increment)
            
        # Expand sample sequence by adding deletions _ and removing insertions
        expanded_sample, p_sample, p_ref, insert_pos = '', 0, 0, {}
        for i in range(len(cigar_expanded)):
            
            if cigar_expanded[i] == "=" or cigar_expanded[i] == "X":
                expanded_sample += sample_seq[p_sample]
                update_err_dist(err_dist, cigar_expanded[i], p_ref in homo_pos)
                p_sample += 1
                p_ref += 1
                
            elif cigar_expanded[i] == "D":
                expanded_sample += "_"
                update_err_dist(err_dist, cigar_expanded[i], p_ref in homo_pos)
                p_ref += 1
                
            elif cigar_expanded[i] == "N":
                expanded_sample += "N"
                p_ref += 1
                
            elif cigar_expanded[i] == "I":
                insert_pos[p_ref] = sample_seq[p_sample] if sample_seq[p_sample] != "T" else "T" 
                update_err_dist(err_dist, cigar_expanded[i], p_ref in homo_pos)
                p_sample += 1
                
#         if len(expanded_sample) != len(ref_seq):
#             print("expansion didn't work!!")
        
        ref_seq = ref_seq.replace("T", "T")
        expanded_sample = expanded_sample.replace("T", "T")
                    
        for i in range(len(ref_seq)):
            if expanded_sample[i] != "N" and ref_seq[i] in base2index:
                    confusion[base2index[expanded_sample[i]], base2index[ref_seq[i]]] += 1
        
        step = 2
        for i in range(len(ref_seq)-step+1):
            kmer = ref_seq[i:i+step]
            sample_kmer = expanded_sample[i:i+step]
            
            ins_flag, inserts = False, []
            for j in range(i+1, i+step):
                if j in insert_pos:
                    ins_flag = True
                    inserts.append((j-i-1, insert_pos[j]))
            if ins_flag == True:
                for loc, insert in inserts:
                    inserted_kmer = sample_kmer[:loc+1] + insert + sample_kmer[loc+1:]
                sample_kmer = inserted_kmer

            if "N" not in sample_kmer and kmer in motif_acc:
                if sample_kmer not in motif_acc[kmer]:
                    motif_acc[kmer][sample_kmer] = 0
                motif_acc[kmer][sample_kmer] += 1
    
    return kmer_dict, motif_acc, confusion, err_dist



@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def main():
    """
    # base-calling stats

    Calculates the basecalling stats from a sam and reference
    """

@main.command()
@click.option(
    "--sam",
    type=click.Path(exists=True, dir_okay=False),
    help="Path to sam file"
)
@click.option(
    "--ref",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--out",
    required=True,
    type=click.Path(dir_okay=False),
    help="Path to reference exported npy file",
)
def evaluate(
    sam,
    ref,
    out
):   
    setup_logging("warning")
    logger.info("Parsing the entire dataset to compute general stats...")
    ref = pysam.Fastafile(ref)

    read_stats, failed_read_c, total_read_c, lowq_c, passed_read_c = get_sample_stats(ref, sam)
    # Name

    # Failed mapping ratio
    failed_map_r = failed_read_c / total_read_c
    # Low quality read ratio
    low_q_r = lowq_c / total_read_c
    # Pass rate
    pass_r = passed_read_c / total_read_c
    # Align ratio
    align_r = read_stats["align_ratio"]
    
    # match rate
    read_stats["acc"][:]
    # mismatch rate
    read_stats["metrics"][:, 0]
    # insetion rate
    read_stats["metrics"][:, 1]
    # deletion rate
    read_stats["metrics"][:, 2]


    # Base quality score distribution
    #read_stats["read_q"]
    # Aligned Read length distribution
    #read_stats["aligned_read_length"]

    # Homopolymer stats
    kmer_profiles, motif_acc, confusion, err_dist = get_kmer_profiles(sam, ref, subsample = 10000000000, kmer_size=3)
    homo_total = err_dist['ho_mis'] + err_dist['ho_ins'] + err_dist['ho_del']
    hetero_total = err_dist['he_mis'] + err_dist['he_ins'] + err_dist['he_del']
    total = homo_total + hetero_total


    homop_profiles = get_homop_profiles(sam, ref, subsample = 10000000000)

    # homopolymer error rate 
    if total > 0:
        ho_r = homo_total / total * 100
    else:
        ho_r = 0
    # mismatch rate
    err_dist['ho_mis']
    # insertion rate
    err_dist['ho_ins']
    # deletion rate
    err_dist['ho_del']


    # heteropolymer error
    if total > 0:
        he_r = hetero_total / total * 100
    else:
        he_r = 0
    # mismatch rate
    err_dist['he_mis']
    # insertion rate
    err_dist['he_ins']
    # deletion rate
    err_dist['he_del']
    

    # Homopolymer error rate for single base
    motif_stats = {}
    for motif in motif_acc:
        total, n_del, n_ins = np.sum(list(motif_acc[motif].values())), 0, 0
        for basecall in motif_acc[motif]:
            if "_" in basecall:
                n_del += motif_acc[motif][basecall]
            if len(basecall) > 2:
                n_ins += motif_acc[motif][basecall]

        deletion = n_del / total
        insertion = n_ins / total
        if not motif in motif_acc[motif]:
            mismatch = (total - n_del - n_ins) / total
        else:
            mismatch = (total - motif_acc[motif][motif] - n_del - n_ins) / total
        motif_stats[motif] = [mismatch, insertion, deletion]

    # Stats for specific motifs
    # motif_stats


    # 5-mer performance specific homopolymers
    confusion_matrix = np.flip(confusion / np.sum(confusion, axis = 0)[None, :], axis = 0)
    basecalled = ["del", "T", "G", "C", "A"] 


    # Read Length
    read_stats["original_read_length"]

    np.savez(out, 
            failed_map_r=failed_map_r, 
            low_q_r=low_q_r,
            pass_r=pass_r, 
            read_length=read_stats["original_read_length"],
            align_r=align_r,
            match_r=read_stats["acc"][:],
            mismatch_r=read_stats["metrics"][:, 0],
            ins_r=read_stats["metrics"][:, 1],
            del_r=read_stats["metrics"][:, 2],
            q_score=read_stats["empirical_q"],
            homop_profile=homop_profiles,
            ho_r=ho_r,
            ho_mis=err_dist['ho_mis'],
            ho_ins = err_dist['ho_ins'],
            ho_del = err_dist['ho_del'],
            he_r=he_r,
            he_mis = err_dist['he_mis'],
            he_ins = err_dist['he_ins'],
            he_del = err_dist['he_del'],
            motif_stats=motif_stats,
            confusion_matrix=confusion_matrix,      
    )


@main.command()
@click.argument(
    "npz",
    required=True,
    nargs=-1,
    type=click.Path(dir_okay=False),
)
@click.option(
    "--outdir",
    required=True,
    type=click.Path(dir_okay=True),
    help="Path to reference exported picture",
)
def plot(
    npz,
    outdir,
):  
    
    setup_logging("info")

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    

    names = [os.path.basename(file).split(".npz")[0] for file in npz]
    npzfiles = [np.load(file, allow_pickle=True) for file in npz]

    # READ LENGTH distribution violin plot
    read_len_l = [(npzfile["read_length"]) for npzfile in npzfiles]
    median_read_len = [np.median(read_len) for read_len in read_len_l]
    max_read_len = [np.max(read_len) for read_len in read_len_l]
    min_read_len = [np.min(read_len) for read_len in read_len_l]

    df_read_len = pd.DataFrame(read_len_l).T
    # Optionally, you can specify column names
    df_read_len.columns = names

    # READ LENGTH VIOLIN PLOT
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.violinplot(data=df_read_len, fill=True, zorder=9)
    for artist in ax.lines:
        artist.set_zorder(10)
    for artist in ax.findobj(PathCollection):
        artist.set_zorder(11)
        
    # ax = sns.stripplot(data=df_read_len, ax=ax ,jitter=0.4, size=1.5)
    ax.set_title("Read length distribution")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    outfile = os.path.join(outdir, "read_length_dist-violin.png")
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()


    # READ LENGTH distribution histogram plot
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.histplot(
        data=df_read_len,
        element="step", fill=False,
        cumulative=True, stat="density", common_norm=False, legend=True,
    )
    ax.set_title("Read length distribution")
    #ax.set_xticklabels(names)
    outfile = os.path.join(outdir, "read_length_dist-hist.png")
    # plt.legend(labels=names)
    # Rotate the labels
    ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()
    
    
    # READ LENGTH distribution kde plot
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.displot(data=df_read_len, kind="kde", label='small', legend=True)
    outfile = os.path.join(outdir, "read_length_dist-kde.png")
    plt.tight_layout()  
    plt.savefig(outfile, dpi=600)
    plt.close()
    
    # READ LENGTH BOXENPLOT
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.boxenplot(data=df_read_len, width_method="linear")
    #ax = sns.displot(data=read_len_l, multiple="dodge", kde=True, label='small', legend=False)
    ax.set_title("Read length distribution")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    outfile = os.path.join(outdir, "read_length_dist-boxenplt.png")
    #plt.legend(labels=names)
    plt.tight_layout()  
    plt.savefig(outfile, dpi=600)
    plt.close()

    # Failed mapping bar plot
    failed_map_l = [100*float(npzfile["failed_map_r"]) for npzfile in npzfiles]
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.barplot(x=names, y=failed_map_l, hue=names)
    ax.set_title("Failed mapping")
    for i in range(len(names)):
        ax.bar_label(ax.containers[i], padding=-32, color="white", fmt='%.2f%%') 
    outfile = os.path.join(outdir, "failed_mapping_reads.png")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()

    # Low quality reads bar plot
    low_q_l = [100*float(npzfile["low_q_r"]) for npzfile in npzfiles]
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.barplot(x=names, y=low_q_l, hue=names)
    ax.set_title("Low quality read (MAPQ < 20; read length < 50)")
    for i in range(len(names)):
        ax.bar_label(ax.containers[i], padding=10, color="black", fmt='%.2f%%') 
    outfile = os.path.join(outdir, "lowQC_read.png")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()
    
    # Passed reads bar plot
    pass_l = [float(100*npzfile["pass_r"]) for npzfile in npzfiles]
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.barplot(x=names, y=pass_l, hue=names)
    ax.set_title("Passed")
    for i in range(len(names)):
        ax.bar_label(ax.containers[i], padding=-32, color="white", fmt='%.2f%%')  
    outfile = os.path.join(outdir, "pass_reads.png")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()

    # Align ratio boxplot
    align_r_l = [100*npzfile["align_r"] for npzfile in npzfiles]
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.boxplot(data=align_r_l)
    ax.set_title("Align ratio")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    outfile = os.path.join(outdir, "alignment_ratio.png")
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()

    # Match ratio boxplot
    match_r_l = [100*npzfile["match_r"] for npzfile in npzfiles]
    median_match = [np.median(match_r) for match_r in match_r_l]

    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.boxplot(data=match_r_l)
    ax.set_title("Match rate")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    outfile = os.path.join(outdir, "match_rate.png")
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()

    # Mismatch ratio boxplot
    mismatch_r_l = [100*npzfile["mismatch_r"] for npzfile in npzfiles]
    median_mismatch = [np.median(mismatch_r) for mismatch_r in mismatch_r_l]

    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.boxplot(data=mismatch_r_l)
    ax.set_title("Mismatch rate")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    outfile = os.path.join(outdir, "mismatch_rate.png")
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()

    # Deletion ratio boxplot
    del_r_l = [100*npzfile["del_r"] for npzfile in npzfiles]
    median_del = [np.median(del_r) for del_r in del_r_l]
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.boxplot(data=del_r_l)
    ax.set_title("Deletion rate")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    outfile = os.path.join(outdir, "deletion_rate.png")
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()

    # Insertion ratio boxplot
    ins_r_l = [100*npzfile["ins_r"] for npzfile in npzfiles]
    median_ins = [np.median(ins_r) for ins_r in ins_r_l]
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.boxplot(data=ins_r_l)
    ax.set_title("Insertion rate")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    outfile = os.path.join(outdir, "insertion_rate.png")
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()

    # Phred Quality Violinplot
    q_score_l = [npzfile["q_score"] for npzfile in npzfiles]
    median_q_score = [np.median(q_score) for q_score in q_score_l]
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.violinplot(data=q_score_l, fill=True, zorder=9)
    for artist in ax.lines:
        artist.set_zorder(10)
    for artist in ax.findobj(PathCollection):
        artist.set_zorder(11)
    # ax = sns.stripplot(data=q_score_l, ax=ax ,jitter=0.4, size=0.5, zorder=8)
    ax.set_title("PHRED Quality score distribution")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    outfile = os.path.join(outdir, "PHRED_dist.png")
    plt.tight_layout()  
    plt.savefig(outfile, dpi=300)
    plt.close()

    # PHRED distribution kde plot
    df_qscore = pd.DataFrame(q_score_l).T
    # Optionally, you can specify column names
    df_qscore.columns = names

    # READ LENGTH BOXENPLOT
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.boxenplot(data=df_qscore, width_method="linear")
    #ax = sns.displot(data=read_len_l, multiple="dodge", kde=True, label='small', legend=False)
    ax.set_title("PHRED Quality score distribution")
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)  
    outfile = os.path.join(outdir, "PHRED-boxenplt.png")
    #plt.legend(labels=names)
    plt.tight_layout()  
    plt.savefig(outfile, dpi=600)
    plt.close()


    
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.0)
    ax = sns.displot(data=df_qscore, kind="kde", label='small', legend=True)
    outfile = os.path.join(outdir, "PHRED-kde.png")
    plt.savefig(outfile, dpi=600)
    plt.close()

    # Homopolymer error rate for A C T G
    # TODO Replace with boxplot at some point
    for homo_len in range(2, 7):
        x, y, hue = [], [], []
        for npzfile in npzfiles:
            prof = npzfile["homop_profile"].item()
            for kmer in sorted(prof.keys()):
                if len(kmer) == homo_len:
                    hue.append(kmer)
                    if prof[kmer][1] > 50:
                        y.append(prof[kmer][0] / prof[kmer][1])
                    else:
                        y.append(0)
        sns.set_context("paper", font_scale = 1.0)
        x = [4*[name] for name in names]
        x = [x for xs in x for x in xs]
        sns.set_theme(style="whitegrid")
        sns.set_context("paper", font_scale = 1.0)
        ax = sns.barplot(x=x, y=y, hue=hue)
        ax.set_title(f"Homopolymer accuracy (length = {homo_len})")
        ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6)

        outfile = os.path.join(outdir, f"homopolymer_acc_{homo_len}.png")
        plt.tight_layout()
        plt.savefig(outfile, dpi=300)
        plt.close()

    # AUC by sorting match rate by PHRED score
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1.15)
    fig, ax1 = plt.subplots(figsize=(6, 8))
    # colorpalettetab10 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    names = [os.path.basename(file).split(".npz")[0] for file in npz]
    aucs = []

    # Using a standard color palette from Seaborn
    palette = sns.color_palette("deep")

    for j, model in enumerate(names):
        dict = {'Match_rate': match_r_l[j], 'Q-score': q_score_l[j]} 
        def integrate(x, y):
            sm = 0
            for i in range(1, len(x)):
                h = x[i] - x[i-1]
                sm += h * (y[i-1] + y[i]) / 2

            return sm
        df = pd.DataFrame(dict)
        df = df.sort_values('Q-score', ascending=True)
        df = df.reset_index()
        coords = {
                'fraction': list(),
                'match_rate': list(),
                'phred_mean': list()
            }
        AUC_STEP = 10
        for i in range(0, len(df), AUC_STEP):
                sub_df = df.loc[i:, :]
                coords['fraction'].append(len(sub_df)/len(df))
                coords['match_rate'].append(np.mean(sub_df['Match_rate']))
                coords['phred_mean'].append(np.min(sub_df['Q-score']))
        
        auc = -integrate(coords['fraction'], coords['match_rate'])
        aucs.append(auc)
        coords = pd.DataFrame(coords)
        #ax2 = ax1.twinx()
        ax1.set_xlabel('Fraction of reads')
        ax1.set_ylabel('Average match rate')
        #ax2.set_ylabel('Minimum read-mean PhredQ score')
        ax1.plot(coords['fraction'], coords['match_rate'], color= palette[j])
        #ax2.plot(coords['fraction'], coords['phred_mean'], color= colorpalettetab10[j], linestyle='--')
    
    # Add legend below the plot
    # Add legend below the plot
    handles = []
    for k, (c, modelname) in enumerate(zip(palette, names)):
        handles.extend([Line2D([0], [0], label=modelname + ' (' + str(round(aucs[k], 2)) + ')', color=c, linewidth=2.5)])
    ax1.legend(handles=handles, edgecolor='black', loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=1)

    # Adjust layout
    #
    plt.subplots_adjust(bottom=0.4) 
    # plt.tight_layout() 
    outfile = os.path.join(outdir, "AUC.png")
    plt.savefig(outfile, dpi=600)
    plt.close()


    # Export median values as table!
    
    median_df = pd.DataFrame([median_match, median_mismatch, median_ins, median_del, median_q_score, median_read_len], 
                                   index=["median_match", "median_mismatch", "median_insertion", "median_deletion",
                                          "median_qscore", "median_read_length"])
    # Optionally, you can specify column names
    median_df.columns = names
    median_df.to_csv(os.path.join(outdir, "median_stats.tsv"), sep='\t', index=True)
    print(median_df)



@main.command()
@click.argument(
    "npz",
    required=True,
    nargs=-1,
    type=click.Path(dir_okay=False),
)
@click.option(
    "--outdir",
    required=True,
    type=click.Path(dir_okay=True),
    help="Path to reference exported picture",
)
def groupplot(
    npz,
    outdir,
):  
    
    setup_logging("info")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # TODO They need to be in correct order seq2squiggle - squigulator - experimental

    if len(npz) != 9:
        logger.error("Something is incorrect. You should select 3 npz files! Please check!")


    names = [os.path.basename(file).split(".npz")[0] for file in npz]
    npzfiles = [np.load(file, allow_pickle=True) for file in npz]
    group_labels = ["seq2squiggle", "squigulator", "experimental"]
    experiment_labels = ["Human Read Mode", "Human Genome Mode", "D. melanogaster Genome Mode"]
    groups = len(npz) // 3


    # GROUPED BOX PLOTS
    align_r_l = [100*npzfile["align_r"] for npzfile in npzfiles]
    match_r_l = [100*npzfile["match_r"] for npzfile in npzfiles]
    mismatch_r_l = [100*npzfile["mismatch_r"] for npzfile in npzfiles]
    del_r_l = [100*npzfile["del_r"] for npzfile in npzfiles]
    ins_r_l = [100*npzfile["ins_r"] for npzfile in npzfiles]
    q_score_l = [npzfile["q_score"] for npzfile in npzfiles]
    # Prepare DataFrame for grouped boxplot
    # Dictionary to store all metrics
    metrics = {
        "Align ratio": align_r_l,
        "Match Rate": match_r_l,
        "Mismatch Rate": mismatch_r_l,
        "Deletion Rate": del_r_l,
        "Insertion Rate": ins_r_l
    }

    for metric_name, metric_data in metrics.items():
        # Prepare DataFrame for grouped boxplot
        data = []
        for i in range(groups):
            for j in range(3):
                group_name = group_labels[j]
                rates = metric_data[i * 3 + j]
                for rate in rates:
                    data.append([group_name, rate, experiment_labels[i]])
        
        df = pd.DataFrame(data, columns=['Group', 'Rate', 'Experiment'])

        sns.set_theme(style="whitegrid")
        sns.set_context("paper", font_scale=1.0)
        
        plt.figure(figsize=(10, 6))
        ax = sns.boxplot(x='Experiment', y='Rate', hue='Group', data=df, showfliers=False)
        ax.set_title(f"{metric_name}")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        outfile = os.path.join(outdir, f"{metric_name.lower().replace(' ', '_')}_grouped.png")
        plt.savefig(outfile, dpi=600)
        plt.close()

    # Prepare DataFrame for grouped violin plot for q_score_l
    data = []
    for i in range(groups):
        for j in range(3):
            group_name = group_labels[j]
            q_scores = q_score_l[i * 3 + j]
            for score in q_scores:
                data.append([group_name, score, experiment_labels[i]])
    
    df = pd.DataFrame(data, columns=['Group', 'Q Score', 'Experiment'])

    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.0)
    
    plt.figure(figsize=(10, 6))
    ax = sns.violinplot(x='Experiment', y='Q Score', hue='Group', data=df, fill=True, zorder=9)
    for artist in ax.lines:
        artist.set_zorder(10)
    for artist in ax.findobj(PathCollection):
        artist.set_zorder(11)
    # Uncomment the following line to add strip plot on top of the violin plot
    # ax = sns.stripplot(x='Experiment', y='Q Score', hue='Group', data=df, jitter=0.4, size=0.5, zorder=8)
    
    ax.set_title("PHRED Quality Score Distribution")
    #ax.set_xticklabels(experiment_labels, ha='right')
    plt.legend(loc='upper right')
    plt.tight_layout()  
    outfile = os.path.join(outdir, "PHRED_dist.png")
    plt.savefig(outfile, dpi=600)
    plt.close()


    ### AUC FACET GRID
    def integrate(x, y):
        sm = 0
        for i in range(1, len(x)):
            h = x[i] - x[i-1]
            sm += h * (y[i-1] + y[i]) / 2
        return sm

    # Prepare data for facet grid
    data = []
    for i in range(groups):
        for j in range(3):
            index = i * 3 + j
            group_name = group_labels[j]
            dataset_name = experiment_labels[i]
            
            data_dict = {'Match_rate': match_r_l[index], 'Q-score': q_score_l[index]}
            df = pd.DataFrame(data_dict)
            df = df.sort_values('Q-score', ascending=True).reset_index(drop=True)
            
            coords = {
                'fraction': list(),
                'match_rate': list(),
                'phred_mean': list()
            }
            AUC_STEP = 10
            for k in range(0, len(df), AUC_STEP):
                sub_df = df.loc[k:, :]
                coords['fraction'].append(len(sub_df) / len(df))
                coords['match_rate'].append(np.mean(sub_df['Match_rate']))
                coords['phred_mean'].append(np.min(sub_df['Q-score']))
            
            auc = -integrate(coords['fraction'], coords['match_rate'])
            coords = pd.DataFrame(coords)
            coords['Group'] = group_name
            coords['Dataset'] = dataset_name
            data.append(coords)
    
    all_data = pd.concat(data)

    # Plotting
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.0)
    
    g = sns.FacetGrid(all_data, col="Dataset", hue="Group", col_wrap=3, height=6, aspect=1.5)
    g.map(sns.lineplot, 'fraction', 'match_rate')

    g.set_axis_labels('Fraction of reads', 'Average match rate')
    # g.add_legend(title='Tool')
    
    # Move legend to the bottom
    g.fig.subplots_adjust(bottom=0.2)
    g.add_legend(title='Tool', bbox_to_anchor=(0.5, -0.1), loc='lower center', ncol=3)


    # Save plot
    outfile = os.path.join(outdir, "AUC_facet.png")
    g.savefig(outfile, dpi=600)
    plt.close()

    ### AUC
    # AUC by sorting match rate by PHRED score
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.0)
    
    # Define a color palette
    color_palette = sns.color_palette("tab10", n_colors=3)
    color_map = {
        'seq2squiggle': color_palette[0],
        'squigulator': color_palette[1],
        'experimental': color_palette[2]
    }
    
    aucs = []
    f, ax1 = plt.subplots(figsize=(8, 5))
    
    # Create a combined legend
    handles = []
    for group_name in group_labels:
        color = color_map[group_name]
        handles.append(Line2D([0], [0], label=group_name, color=color))
    
    for i in range(groups):
        for j in range(3):
            index = i * 3 + j
            group_name = group_labels[j]
            dataset_name = experiment_labels[i]
            color = color_map[group_name]
            
            data_dict = {'Match_rate': match_r_l[index], 'Q-score': q_score_l[index]} 
            def integrate(x, y):
                sm = 0
                for k in range(1, len(x)):
                    h = x[k] - x[k-1]
                    sm += h * (y[k-1] + y[k]) / 2
                return sm
            
            df = pd.DataFrame(data_dict)
            df = df.sort_values('Q-score', ascending=True)
            df = df.reset_index(drop=True)
            coords = {
                    'fraction': list(),
                    'match_rate': list(),
                    'phred_mean': list()
                }
            AUC_STEP = 10
            for k in range(0, len(df), AUC_STEP):
                    sub_df = df.loc[k:, :]
                    coords['fraction'].append(len(sub_df) / len(df))
                    coords['match_rate'].append(np.mean(sub_df['Match_rate']))
                    coords['phred_mean'].append(np.min(sub_df['Q-score']))
            
            auc = -integrate(coords['fraction'], coords['match_rate'])
            aucs.append(auc)
            coords = pd.DataFrame(coords)
            
            ax1.plot(coords['fraction'], coords['match_rate'], label=f"{dataset_name} ({group_name})", color=color, linestyle='-' if i == 0 else '--' if i == 1 else ':')
    
    ax1.set_xlabel('Fraction of reads')
    ax1.set_ylabel('Average match rate')
    
    # Add combined legend
    ax1.legend(handles=handles, edgecolor='black', bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
    
    plt.tight_layout()
    outfile = os.path.join(outdir, "AUC.png")
    plt.savefig(outfile, dpi=300)
    plt.close()


    # Calculate AUC
    aucs = []
    for i in range(groups):
        for j in range(3):
            index = i * 3 + j
            group_name = group_labels[j]
            dataset_name = experiment_labels[i]
            
            data_dict = {'Match_rate': match_r_l[index], 'Q-score': q_score_l[index]} 
            def integrate(x, y):
                sm = 0
                for k in range(1, len(x)):
                    h = x[k] - x[k-1]
                    sm += h * (y[k-1] + y[k]) / 2
                return sm
            
            df = pd.DataFrame(data_dict)
            df = df.sort_values('Q-score', ascending=True)
            df = df.reset_index(drop=True)
            coords = {
                    'fraction': list(),
                    'match_rate': list(),
                    'phred_mean': list()
                }
            AUC_STEP = 10
            for k in range(0, len(df), AUC_STEP):
                    sub_df = df.loc[k:, :]
                    coords['fraction'].append(len(sub_df) / len(df))
                    coords['match_rate'].append(np.mean(sub_df['Match_rate']))
                    coords['phred_mean'].append(np.min(sub_df['Q-score']))
            
            auc = -integrate(coords['fraction'], coords['match_rate'])
            aucs.append([group_name, dataset_name, auc])

    # Create DataFrame for AUC values
    auc_df = pd.DataFrame(aucs, columns=['Tool', 'Dataset', 'AUC'])

    # Grouped barplot
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.0)
    
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x="Dataset", y="AUC", hue="Tool", data=auc_df)
    ax.set_title("AUC Values for Each Tool Across Datasets")
    plt.legend(title='Tool', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    outfile = os.path.join(outdir, "AUC_grouped_barplot.png")
    plt.savefig(outfile, dpi=300)
    plt.close()


    # Passed
    pass_l = [100*npzfile["pass_r"] for npzfile in npzfiles]
    failed_map_l = [100*float(npzfile["failed_map_r"]) for npzfile in npzfiles]
    low_q_l = [100*float(npzfile["low_q_r"]) for npzfile in npzfiles]

    metrics = {
        "Passed reads": pass_l,
        "Failed reads": failed_map_l,
        "Low quality reads": low_q_l,
        # Add more metrics here if needed
    }

    for metric_name, metric_data in metrics.items():
        # Prepare DataFrame for grouped barplot
        data = []
        for i in range(groups):
            for j in range(3):
                group_name = group_labels[j]
                rate = metric_data[i * 3 + j]  # Access the rate directly
                data.append([group_name, rate, experiment_labels[i]])
        
        df = pd.DataFrame(data, columns=['Group', 'Rate', 'Experiment'])

        # Plotting
        sns.set_theme(style="whitegrid")
        sns.set_context("paper", font_scale=1.0)
        
        plt.figure(figsize=(10, 6))
        ax = sns.barplot(x="Experiment", y="Rate", hue="Group", data=df)
        ax.set_title(f"{metric_name}")
        plt.legend(title='Tool', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        outfile = os.path.join(outdir, f"{metric_name.lower().replace(' ', '_')}_grouped.png")
        plt.savefig(outfile, dpi=300)
        plt.close()

@main.command()
@click.argument(
    "npz",
    required=True,
    nargs=-1,
    type=click.Path(dir_okay=False),
)
@click.option(
    "--outdir",
    required=True,
    type=click.Path(dir_okay=True),
    help="Path to reference exported picture",
)
def groupplotfigure(
    npz,
    outdir,
):  
    setup_logging("info")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if len(npz) != 9:
        logger.error("Something is incorrect. You should select 3 npz files! Please check!")

    names = [os.path.basename(file).split(".npz")[0] for file in npz]
    npzfiles = [np.load(file, allow_pickle=True) for file in npz]
    group_labels = ["seq2squiggle", "squigulator", "experimental"]
    experiment_labels = ["Human\nRead Mode", "Human\nGenome Mode", "D. melanogaster\nGenome Mode"]
    groups = len(npz) // 3

    # Extract metrics
    align_r_l = [100*npzfile["align_r"] for npzfile in npzfiles]
    match_r_l = [100*npzfile["match_r"] for npzfile in npzfiles]
    mismatch_r_l = [100*npzfile["mismatch_r"] for npzfile in npzfiles]
    del_r_l = [100*npzfile["del_r"] for npzfile in npzfiles]
    ins_r_l = [100*npzfile["ins_r"] for npzfile in npzfiles]
    q_score_l = [npzfile["q_score"] for npzfile in npzfiles]
    pass_l = [100*npzfile["pass_r"] for npzfile in npzfiles]
    read_len_l = [(npzfile["read_length"]) for npzfile in npzfiles]

    # Combine AUC computation in one step
    aucs = []
    data = []
    for i in range(groups):
        for j in range(3):
            index = i * 3 + j
            group_name = group_labels[j]
            dataset_name = experiment_labels[i]
            
            data_dict = {'Match_rate': match_r_l[index], 'Q-score': q_score_l[index]} 
            df = pd.DataFrame(data_dict)
            df = df.sort_values('Q-score', ascending=True).reset_index(drop=True)
            
            coords = {'fraction': [], 'match_rate': [], 'phred_mean': []}
            AUC_STEP = 10
            for k in range(0, len(df), AUC_STEP):
                sub_df = df.loc[k:, :]
                coords['fraction'].append(len(sub_df) / len(df))
                coords['match_rate'].append(np.mean(sub_df['Match_rate']))
                coords['phred_mean'].append(np.min(sub_df['Q-score']))
            
            auc = -integrate(coords['fraction'], coords['match_rate'])
            aucs.append([group_name, dataset_name, auc])
            coords = pd.DataFrame(coords)
            data.append(coords)

    auc_df = pd.DataFrame(aucs, columns=['Tool', 'Dataset', 'AUC'])

    # Prepare DataFrame for grouped boxplot
    def prepare_boxplot_data(metric_data, metric_name):
        data = []
        for i in range(groups):
            for j in range(3):
                group_name = group_labels[j]
                rates = metric_data[i * 3 + j]
                if isinstance(rates, np.ndarray) or isinstance(rates, list):
                    for rate in rates:
                        data.append([group_name, rate, experiment_labels[i]])
                else:
                    data.append([group_name, rates, experiment_labels[i]])
        df = pd.DataFrame(data, columns=['Group', 'Rate', 'Dataset'])
        return df

    align_df = prepare_boxplot_data(align_r_l, "Align ratio")
    match_df = prepare_boxplot_data(match_r_l, "Match rate")
    mismatch_df = prepare_boxplot_data(mismatch_r_l, "Mismatch rate")
    deletion_df = prepare_boxplot_data(del_r_l, "Deletion rate")
    insetion_df = prepare_boxplot_data(ins_r_l, "Insertion rate")
    
    # Prepare DataFrame for grouped bar plot
    pass_data = []
    for i in range(groups):
        for j in range(3):
            group_name = group_labels[j]
            rate = pass_l[i * 3 + j]  # Access the rate directly
            pass_data.append([group_name, rate, experiment_labels[i]])
    pass_df = pd.DataFrame(pass_data, columns=['Group', 'Rate', 'Dataset'])

    # Prepare data for Q-score violin plot
    q_data = []
    for i in range(groups):
        for j in range(3):
            group_name = group_labels[j]
            q_scores = q_score_l[i * 3 + j]
            for score in q_scores:
                q_data.append([group_name, score, experiment_labels[i]])
    q_df = pd.DataFrame(q_data, columns=['Group', 'Q Score', 'Dataset'])


    # Prepare data for read length
    read_len_data = []
    for i in range(groups):
        for j in range(3):
            group_name = group_labels[j]
            read_lens = read_len_l[i * 3 + j]
            for read_l in read_lens:
                read_len_data.append([group_name, read_l, experiment_labels[i]])
    readlength_df = pd.DataFrame(read_len_data, columns=['Group', 'Read length', 'Dataset'])

    # Plotting
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.5)
    
    fig, axes = plt.subplots(3, 2, figsize=(15, 18))
    # fig.suptitle('Comparison of Sequencing Metrics')

    sns.barplot(ax=axes[0, 0], x="Dataset", y="Rate", hue="Group", data=pass_df)
    axes[0, 0].set_title('Passed reads')
    axes[0, 0].text(-0.1, 1.1, 'A', transform=axes[0, 0].transAxes, 
                    size=20, weight='bold')
    axes[0, 0].set(xlabel='')

    sns.boxplot(ax=axes[0, 1], x='Dataset', y='Rate', hue='Group', data=align_df, showfliers=False)
    axes[0, 1].set_title('Align ratio')
    axes[0, 1].text(-0.1, 1.1, 'B', transform=axes[0, 1].transAxes, 
                    size=20, weight='bold')
    axes[0, 1].set(xlabel='')

    sns.boxplot(ax=axes[1, 0], x='Dataset', y='Rate', hue='Group', data=match_df, showfliers=False)
    axes[1, 0].set_title('Match rate')
    axes[1, 0].text(-0.1, 1.1, 'C', transform=axes[1, 0].transAxes, 
                    size=20, weight='bold')
    axes[1, 0].set(xlabel='')
    axes[1, 0].set_ylim(bottom=0)  

    sns.boxplot(ax=axes[1, 1], x='Dataset', y='Rate', hue='Group', data=mismatch_df, showfliers=False)
    axes[1, 1].set_title('Mismatch rate')
    axes[1, 1].text(-0.1, 1.1, 'D', transform=axes[1, 1].transAxes, 
                    size=20, weight='bold')
    axes[1, 1].set(xlabel='')

    sns.violinplot(ax=axes[2, 0], x='Dataset', y='Q Score', hue='Group', data=q_df, fill=True)
    axes[2, 0].set_title('PHRED Quality Score Distribution')
    axes[2, 0].text(-0.1, 1.1, 'E', transform=axes[2, 0].transAxes, 
                    size=20, weight='bold')
    axes[2, 0].set(xlabel='')

    sns.barplot(ax=axes[2, 1], x="Dataset", y="AUC", hue="Tool", data=auc_df)
    axes[2, 1].set_title('AUC Values')
    axes[2, 1].text(-0.1, 1.1, 'F', transform=axes[2, 1].transAxes, 
                    size=20, weight='bold')
    axes[2, 1].set(xlabel='')

    

    # Remove individual legends and add a shared legend
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=True, fontsize='large', bbox_to_anchor=(0.5, 0.005))
    for ax in axes.flatten():
        ax.get_legend().remove()

    # Adjust layout
    plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust the rect to leave space at the bottom
    plt.subplots_adjust(top=0.95, hspace=0.3)
    
    #plt.tight_layout()
    outfile = os.path.join(outdir, "Figure_03.png")
    plt.savefig(outfile, dpi=600)
    plt.close()


    # Plotting
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.5)

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    sns.boxplot(ax=axes[0, 0], x='Dataset', y='Rate', hue='Group', data=mismatch_df, showfliers=False)
    axes[0, 0].set_title('Mismatch rate')
    axes[0, 0].text(-0.1, 1.05, 'A', transform=axes[0, 0].transAxes, 
                    size=20, weight='bold')
    axes[0, 0].set(xlabel='')

    sns.boxplot(ax=axes[0, 1], x='Dataset', y='Rate', hue='Group', data=deletion_df, showfliers=False)
    axes[0, 1].set_title('Deletion rate')
    axes[0, 1].text(-0.1, 1.05, 'B', transform=axes[0, 1].transAxes, 
                    size=20, weight='bold')
    axes[0, 1].set(xlabel='')

    sns.boxplot(ax=axes[1, 0], x='Dataset', y='Rate', hue='Group', data=insetion_df, showfliers=False)
    axes[1, 0].set_title('Insertion rate')
    axes[1, 0].text(-0.1, 1.05, 'C', transform=axes[1, 0].transAxes, 
                    size=20, weight='bold')
    axes[1, 0].set(xlabel='')

    sns.violinplot(ax=axes[1, 1], x='Dataset', y='Read length', hue='Group', data=readlength_df, fill=True, log_scale=True)
    axes[1, 1].set_title('Read length distribution')
    axes[1, 1].text(-0.1, 1.05, 'D', transform=axes[1, 1].transAxes, 
                    size=20, weight='bold')
    axes[1, 1].set(xlabel='')

    
    # Remove individual legends and add a shared legend
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=True, fontsize='large', bbox_to_anchor=(0.5, 0.005))
    for ax in axes.flatten():
        ax.get_legend().remove()

    # Adjust layout
    plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust the rect to leave space at the bottom
    plt.subplots_adjust(top=0.95, hspace=0.3)
    
    outfile = os.path.join(outdir, "SuppFigure02.png")
    plt.savefig(outfile, dpi=600)
    plt.close()



    # Prepare data for facet grid
    data = []
    for i in range(groups):
        for j in range(3):
            index = i * 3 + j
            group_name = group_labels[j]
            dataset_name = experiment_labels[i]
            
            data_dict = {'Match_rate': match_r_l[index], 'Q-score': q_score_l[index]}
            df = pd.DataFrame(data_dict)
            df = df.sort_values('Q-score', ascending=True).reset_index(drop=True)
            
            coords = {
                'fraction': list(),
                'match_rate': list(),
                'phred_mean': list()
            }
            AUC_STEP = 10
            for k in range(0, len(df), AUC_STEP):
                sub_df = df.loc[k:, :]
                coords['fraction'].append(len(sub_df) / len(df))
                coords['match_rate'].append(np.mean(sub_df['Match_rate']))
                coords['phred_mean'].append(np.min(sub_df['Q-score']))
            
            auc = -integrate(coords['fraction'], coords['match_rate'])
            coords = pd.DataFrame(coords)
            coords['Group'] = group_name
            coords['Dataset'] = dataset_name
            data.append(coords)
    
    all_data = pd.concat(data)

    # Plotting
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.5)

    fig, axes = plt.subplots(nrows=groups, ncols=1, figsize=(10, 6 * groups))

    for i, dataset_name in enumerate(experiment_labels):
        ax = axes[i]
        sns.lineplot(ax=ax, data=all_data[all_data['Dataset'] == dataset_name], x='fraction', y='match_rate', hue='Group', linewidth=2.5)
        ax.set_title(dataset_name)
        ax.set_xlabel('Fraction of reads')
        ax.set_ylabel('Average match rate')
        ax.text(-0.1, 1.05, chr(65 + i), transform=ax.transAxes, weight='bold')

    # Remove individual legends and add a shared legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=True, bbox_to_anchor=(0.5, 0.005))
    for ax in axes:
        ax.get_legend().remove()

    # Adjust layout
    plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust the rect to leave space at the bottom
    plt.subplots_adjust(hspace=0.3)

    # Save plot
    outfile = os.path.join(outdir, "SuppFigure03.png")
    plt.savefig(outfile, dpi=600)
    plt.close()



    # Set up plotting parameters
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.5)
 
    x, y, hue = [], [], []
    # only take human GM mode
    npzfiles = npzfiles[3:6]
    names = names[3:6]

    # Choose a different color palette
    # colors = sns.color_palette('Set2', n_colors=4)

    # Create a figure and axes for subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    # fig.suptitle('Homopolymer Accuracy for different Kmer Lengths', fontsize=16)

    # Iterate over the kmer lengths 3, 4, and 5
    for i, homo_len in enumerate([3, 4, 5, 6]):
        x, y, hue = [], [], []
        for npzfile in npzfiles:
            prof = npzfile["homop_profile"].item()
            for kmer in sorted(prof.keys()):
                if len(kmer) == homo_len:
                    hue.append(kmer)
                    if prof[kmer][1] > 50:
                        y.append(prof[kmer][0] / prof[kmer][1] * 100)  # Convert to percentage
                    else:
                        y.append(0)
        # Flatten the x list
        x = [len(set(hue)) * [name] for name in names]
        x = [item for sublist in x for item in sublist]

        # Create a DataFrame from extracted data
        data = {
            'Dataset': x,
            'Accuracy': y,
            'Homopolymer': hue
        }
        df = pd.DataFrame(data)

        # Select the appropriate subplot
        ax = axes[i // 2, i % 2]
        
        # Plot using seaborn
        sns.barplot(ax=ax, x='Dataset', y='Accuracy', hue='Homopolymer', data=df)
        ax.set_title(f"Homopolymer accuracy (length = {homo_len})")
        ax.set_xticklabels(names, rotation=45, ha='right')
        ax.set_ylabel('Average match rate')
        ax.set_xlabel('')
        ax.set_ylim(0, 100)  # Set the same scale for all plots
        # Remove the individual legend
        ax.get_legend().remove()

    # Add a shared legend
    handles, labels = axes[0, 0].get_legend_handles_labels()  # Use the first subplot's legend
    fig.legend(handles, labels, title='Homopolymer', loc='lower center', bbox_to_anchor=(0.5, 0.02), ncol=4, frameon=True)

    plt.subplots_adjust(hspace=0.4, wspace=0.4)
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])  # Adjust the rect to leave space for the legend

    

    # Save the plot
    outfile = os.path.join(outdir, "SuppFigureHOMOPOLYMER.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
        



def integrate(x, y):
    sm = 0
    for i in range(1, len(x)):
        h = x[i] - x[i-1]
        sm += h * (y[i-1] + y[i]) / 2
    return sm




def setup_logging(verbosity):
    logging_levels = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
    }

    # Configure logging.
    logging.captureWarnings(True)
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    warnings_logger = logging.getLogger("py.warnings")

    # Formatters for file vs console:
    console_formatter = logging.Formatter(
        "{name} {levelname} {asctime}: {message}", style="{", datefmt="%H:%M:%S"
    )

    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(logging_levels[verbosity.lower()])
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)
    warnings_logger.addHandler(console_handler)

    warnings.filterwarnings("ignore")
    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("seaborn").setLevel(logging.ERROR)
    logging.getLogger("fsspec").setLevel(logging.WARNING)
    logging.getLogger("github").setLevel(logging.WARNING)
    logging.getLogger("h5py").setLevel(logging.WARNING)
    logging.getLogger("numba").setLevel(logging.WARNING)
    logging.getLogger("pytorch_lightning").setLevel(logging.WARNING)
    logging.getLogger("torch").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)


if __name__ == "__main__":
    main()
