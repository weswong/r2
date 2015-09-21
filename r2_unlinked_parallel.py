from pandas import DataFrame, read_csv
import itertools
import numpy as np
from collections import Counter
import sys
import multiprocessing as mp
import math

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)


def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]
    
    
def allele_count(df_column):
    counts = Counter(df_column[3:])
    counts.pop('.', 0)
    count = len(counts.keys())
    return count

def major_find(df_row):
    line = df_row[3:]
    counts = Counter(line)
    #print counts
    # removes and returns the value of '.'
    missing_count = counts.pop('.', 0)
    major = sorted(counts.items(), key=lambda x: x[1])[1][0]
    return major

def major_prop_find(df_row):
    line = df_row[3:]
    #print [x for x in df_row[0:3]]
    counts = Counter(line)
    # removes and returns the value of '.'
    missing_count = counts.pop('.', 0)
    total = float(sorted(counts.items(), key=lambda x: x[1])[0][1] + sorted(counts.items(), key=lambda x: x[1])[1][1])
    major_prop = sorted(counts.items(), key=lambda x: x[1])[1][1] / total
    return major_prop


def file_prep(file):
    df = DataFrame(read_csv(file, sep = '\t'))
    df.drop(df[df.apply(allele_count, axis = 1) != 2].index, inplace = True)
    major_freqs = df.apply(major_prop_find, axis = 1)
    major_alleles = df.apply(major_find, axis =1 )
    df.insert(3,'major_freqs', major_freqs)
    df.insert(3,'major_alleles', major_alleles)
    df = df.transpose()
    
    
    chrom, chrom_idx = np.unique(df.loc['chrom'], return_index=True)
    
    super_missing_df = df == '.'
    
    chromosome_dict = {}
    for number in np.unique(df.loc['chrom']):
        chromosome_dict[number] = df.loc['chrom'][df.loc['chrom'] == number].index
    return df, super_missing_df, chromosome_dict


def interlocus_r2(idx1, idx2, df, super_missing_df):
    haplotype= df[idx1] + df[idx2]
    missing_sample = haplotype[haplotype.str.contains('\.') == True]
    haplotype.drop(missing_sample.index, inplace=True)
    s = float(len(haplotype[5:]))
    
    loci1 = df[idx1].drop(missing_sample.index, inplace=False)[5:]
    loci2 = df[idx2].drop(missing_sample.index, inplace=False)[5:]

    missing_points = (len(df) - len(haplotype))/float(len(df))
    
    major_prop_idx1, major_prop_idx2 =Counter(loci1)[df[idx1].loc['major_alleles']] / s, Counter(loci2)[df[idx2].loc['major_alleles']] / s
    #major_prop_idx1, major_prop_idx2 =major_prop_find(loci1), major_prop_find(loci2)
    
    #major_prop_idx1, major_prop_idx2 = df[idx1].loc['major_freqs'], df[idx2].loc['major_freqs']
    major_hap = df[idx1].loc['major_alleles'] + df[idx2].loc['major_alleles']
    
    pos1 = str(df[idx1].loc['chrom']) + '.' + str(df[idx1].loc['pos'])
    pos2 = str(df[idx2].loc['chrom']) + '.' + str(df[idx2].loc['pos'])  
    
    comparison = pos1 + ':' + pos2
    print comparison
    #print major_find(haplotype_df[idx1])
    #print major_hap
    

    def linkage_disequillibria(haplotype, pA, pB, major_hap, s):
        expected_hap_freq =pA * pB
        #pAB
        observed_hap_freq = Counter(haplotype)[major_hap] / s
        
        #print 'expected hap freq = ' + str(expected_hap_freq)
        #print 'observed hap freq = ' + str(observed_hap_freq)
        
        #D = PAB - PaPb
        D = observed_hap_freq - expected_hap_freq
        
        if D < 0:
            Dmax = np.min([-pA*pB, -(1-pA)*(1-pB)])
        else: # D > 0
            Dmax = np.max([pA*(1-pB), (1-pA)*pB])

        Dprime = D/float(Dmax)
        #corrected r_square to handle allele frequency
        # D^2 / (Pa(1-Pa)Pb(1-Pb))
        r_square = (D**2) / (pA*pB*(1-pA)*(1-pB))
        
        #log.write(str(comparison) + '\t' + str(pA)  + '\t' + str(pB) + '\t' + str(expected_hap_freq) + '\t' + str(observed_hap_freq) + '\t' + str(r_square) + '\t' + str(D) + '\t' + str(Dprime) + '\n')
        
        return r_square, Dprime, observed_hap_freq
    
    r_square, Dprime, observed_hap_freq = linkage_disequillibria(haplotype, major_prop_idx1, major_prop_idx2,major_hap, s)
    
    #print 'pA = ' + str(major_prop_idx1) + ', ' + str(major_find(haplotype_df[idx1]))
    #print 'pB = ' + str(major_prop_idx2) + ', ' + str(major_find(haplotype_df[idx2]))
    #print 'AB = ' + major_hap
    return str(r_square), str(Dprime), comparison, str(major_prop_idx1), str(major_prop_idx2), str(observed_hap_freq), str(s)

def calculate(file, processes):
    fout = open('linkage.txt', 'w')
    fout.write('comparison\tr2\td\'pA\tpB\tobserved_hap_freq\ts\n')
    df, super_missing_df, chromosome_dict = file_prep(file)
    r_square, dprime, position = [],[], []
    chromosome_combos = list(itertools.combinations(chromosome_dict.keys(), 2))
    total = len(chromosome_combos)
    chunk_size = int(math.ceil(total/float(processes)))
    if chunk_size == 0:
        chunk_size = 1
    slices = chunks(chromosome_combos, chunk_size)
    
    results = [pool.apply_async(calculate_chunks, args=(chromosome_dict,df,super_missing_df, slice)) for slice in slices]
    output = []
    for p in results:
        print p.get()
        output += p.get()
    print len(output)
    
    for result in output:
        fout.write(('\t').join(result) + '\n')

def calculate_chunks(chromosome_dict,df,super_missing_df, slice):
    values = []
    print slice
    for x in slice:
        chr1, chr2 = x[0], x[1]
        for idx1, idx2 in itertools.product(chromosome_dict[chr1], chromosome_dict[chr2]):
            #print idx1, idx2
            r, d, p, pA, pB, observed_hap_freq, s = interlocus_r2(idx1, idx2, df, super_missing_df)
            values.append((p,r,d,pA, pB, observed_hap_freq, s))
    return values

if __name__ == '__main__':
    mp.freeze_support()
            
    file = sys.argv[1]
    processes = int(sys.argv[2])

    pool = mp.Pool(processes = processes)
    calculate(file, processes)

