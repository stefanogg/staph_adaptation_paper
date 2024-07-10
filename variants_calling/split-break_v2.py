#!/usr/bin/env python3

"""
Script to extract breakpoints from bam files containing
split reads
"""

import sys
import argparse
import re
import os
import shutil
import pandas as pd
import numpy as np
from collections import OrderedDict

def parse_input_file ( input_dict, input_file ):
    """
    Parse the tab separated input file and create
    a dictionary (key = isolate, value = splitters.bam)
    """
    with open( input_file , 'r' ) as f:
        for line in f:
            data = line.split()
            input_dict[data[0]] = data[1]


def make_directories( input_file, input_dict , refpath, bed_file = False):
    """
    Create the directory structure, inside a pre-existent main
    directory with the name of the analysis (or name of the pair):
    one directory for each isolate and a directory for the reference
    """
    input_file_name = os.path.basename(input_file)
    if not os.path.exists(input_file_name):
        shutil.copy( input_file, './isolates.tab')
    for i in input_dict:
        os.mkdir(i)
    os.mkdir('reference')
    os.chdir('reference')
    shutil.copy(refpath , 'ref.gbk')
    cmd = '~/perl5/bin/split-gbk.py ref.gbk'
    # os.system( cmd )
    splitters_bam_dirname = os.path.dirname( [ v for v in list(input_dict.values()) ][0] )
    print(splitters_bam_dirname)
    ref_fasta_dirname = os.path.dirname( splitters_bam_dirname )
    print(ref_fasta_dirname)
    #shutil.copy( ref_fasta_dirname + '/reference/ref.fa', '.' )
    # shutil.copy( ref_fasta_dirname + '/reference/ref.fa.fai', '.' )
    cmd = 'any2fasta ref.gbk > ref.fa'
    print(cmd)
    os.system(cmd)
    cmd = 'samtools faidx ref.fa'
    print(cmd)
    os.system(cmd)

    if bed_file:
    # generates regions
        cmd = 'bedtools makewindows -g ref.fa.fai -n 1 | bedtools subtract -a - -b ' + bed_file + ' > non_masked_regions.bed'
        print(cmd)
        os.system(cmd)
    os.chdir('..')
    return splitters_bam_dirname

def run_mpileup( isolate_name, mask_regions = False ):
    cmd = 'samtools mpileup -AB -C 0 ' + isolate_name + '.splitters.bam' + \
          ' > ' + isolate_name + '.splitters.pileup'
    if mask_regions:
        cmd = 'samtools mpileup -AB -C 0 -l ../reference/non_masked_regions.bed ' + \
              isolate_name + '.splitters.bam > ' + isolate_name + '.splitters.pileup'
    print(cmd)
    os.system( cmd )
    
def draw_graph( isolate_name ):
    with open( isolate_name + '.splitters.pileup', 'r') as f:
        with open( isolate_name + '.splitters.breaks.graph', 'w' ) as fout:
            for line in f:
                cols = line.split()
                if len(cols) >= 5:
                    breaks = ''
                    depth = ''
                    ratio_rounded = ''
                    if re.search( '[\^|\$]' , cols[4] ):
                        breaks = int( len( re.sub('[^\^|\$]','',cols[4]) ) )
                        depth = int ( cols[3] )
                        ratio = float( breaks ) / float ( depth )
                        ratio_rounded = round( ratio, 2 )
                    breaks_graph = re.sub('[^\^|\$]','', cols[4])
                    new_line = cols[0] + '\t' + cols[1] + '\t' + str( breaks ) + '\t' + str( ratio_rounded ) + '\t' + breaks_graph + '\n' 
                else:
                    new_line = line
                fout.writelines( new_line )
    print(isolate_name + '.splitters.break.graph saved.')

def draw_depth_graph ( isolate_name ):
    cmd = 'samtools depth -a ' + isolate_name + '.splitters.bam > ' + \
          isolate_name + '.splitters.depth'
    print(cmd)
    os.system ( cmd )
    with open( isolate_name + '.splitters.depth', 'r') as f:
        with open( isolate_name + '.splitters.depth.graph', 'w' ) as fout:
            for line in f:
                cols = line.split()
                chrom = cols[0]
                pos = cols[1]
                depth = int( cols[2] )
                graph = '*' * depth
                new_seq = [ chrom, pos, str( depth ), graph ]
                new_line = '\t'.join( map( str, new_seq ) ) + '\n'
                fout.writelines( new_line )

def write_table( isolate_name ):
    with open( isolate_name + '.breaks.tab', 'w' ) as fout:
        with open( isolate_name + '.splitters.pileup', 'r') as f:
            for line in f:
                cols = line.split()
                if re.search( '\^', line ):
                    chrom = cols[0]
                    pos = cols[1]
                    depth = int( cols[3])
                    breaks = int( len( re.sub('[^\^]','',cols[4]) ) )
                    ratio = float( breaks ) / float( depth )
#                    print ratio
                    new_seq = [ isolate_name, chrom, pos, str( depth ), 'start', str( breaks ), format( ratio, '.2f' ) ]
                    new_line = '\t'.join( map( str, new_seq ) ) + '\n'
                    fout.writelines( new_line )
                elif re.search( '\$', line ):
                    chrom = cols[0]
                    pos = cols[1]
                    depth = int( cols[3])
                    breaks = int( len( re.sub('[^\$]','',cols[4]) ) )
                    ratio = float (breaks) / float ( depth )
                    new_seq = [ isolate_name, chrom, pos, str( depth ), 'stop', str( breaks ), format( ratio, '.2f' ) ]
                    new_line = '\t'.join( map( str, new_seq ) ) + '\n'
                    fout.writelines( new_line )
    print(isolate_name + '.breaks.tab saved.')

def get_depth( isolate_name, input_dict ):
    with open( isolate_name + '.breaks.tab', 'r') as f, open( isolate_name + '.pos.tab', 'w') as fout:
        for line in f:
            new_line = line.split()[1:3]
            print(new_line)
            fout.writelines( '\t'.join( new_line ) + '\n' )
    splitters_bam_dirname = os.path.dirname( input_dict[ isolate_name ] )
    cmd = 'samtools depth -b ' + isolate_name + '.pos.tab ' + splitters_bam_dirname + '/' + isolate_name + \
          '.bam > ' + isolate_name + '.depth'
    print(cmd)
    os.system( cmd )

def filter_table( isolate_name, min_depth, min_ratio ):
    with open( isolate_name + '.breaks.filter.tab', 'w') as fout:
        with open( isolate_name + '.breaks.tab', 'r' ) as f:
            for line in f:
                break_depth = int( line.split()[5] )
                ratio = float( line.split()[6] )
                if break_depth >= min_depth and ratio >= min_ratio:
                    fout.writelines( line )
    print(isolate_name + '.breaks.filter.tab saved.')

def filter_breaks(isolate_name, min_depth):
    for b in ['start','stop']:
        with open(isolate_name + '.breaks.tab', 'r') as f:
            with open(isolate_name + '.' + b + '.bed', 'w') as fout:
                for line in f:
                    break_depth = int(line.split()[5])
                    break_type = line.split()[4]
                    if break_depth >= min_depth and break_type == b:
                        cols=line.split()
                        chrom=cols[1]
                        end=cols[2]
                        start=str(int(end) - 1)
                        depth=cols[3]
                        break_type=cols[4]
                        break_depth=cols[5]
                        ratio=cols[6]
                        new_line='\t'.join([chrom,start, end, depth, break_type, break_depth, ratio]) + '\n'
                        fout.writelines(new_line)

def get_insertion_pairs(isolate_name):
    a = isolate_name + '.start.bed'
    b = isolate_name + '.stop.bed'
    out = isolate_name + '.insertion.pairs.bed'
    cmd = 'bedtools window -a ' + a + ' -b ' + b + ' -r 10 -l 0 > ' + out
    print('Running: {}'.format(cmd))
    os.system(cmd)

def mask_internal_reference_positions(isolate_name, internal_ref):
    os.chdir(isolate_name)
    a = isolate_name + '.insertion.pairs.bed'
    b = '../' + internal_ref + '/' + internal_ref + '.insertion.pairs.bed'
    out = isolate_name + '.mask.insertion.pairs.bed'
    cmd = 'bedtools subtract -a ' + a + ' -b ' + b + ' > ' + out
    print('Running: {}'.format(cmd))
    os.system(cmd)
    os.chdir('..')

def intersect_tables( input_dict, name, add_splitters_depth = False ):
    '''
    merge filtered tables from all isolates (n_isolates) and
    keep positions that are found in < n_isolates.
    Write intersect table with columns: ISOLATE, chrom, pos, 
    depth (of split reads), break type (start or stop), breaks depth, 
    ratio breaks depth / split depth, n (number of isolates with position)
    '''
    n_isolates = len( input_dict) 
    print('n_isolates is: ' + str( n_isolates ))
    
    def read_table_pandas( filename ):
        names = ['ISOLATE','CHROM','POS','BREAK_DEPTH','BREAK_TYPE','BREAK_COUNT','RATIO']
        tabl = pd.read_table( filename, names = names )
        return tabl

    tabl_list = []
    for isolate_name in input_dict:
        df = read_table_pandas( isolate_name + '/' + isolate_name + '.breaks.tab' )
        tabl_list.append(df)
    breaks_df = pd.concat( tabl_list, axis = 0 )
    breaks_df['pair_name'] = name
    out = breaks_df
    print('{n} breaks.tab files concatenated'.format(n = n_isolates))

    if add_splitters_depth:
        print('Merging with splitters coverage data ...')
        # add depth
        tabl_list = []
        names = ['CHROM', 'POS', 'SPLITTERS_DEPTH']
        for isolate_name in input_dict:
            df = pd.read_table( isolate_name + '/' + isolate_name + '.splitters.depth', names = names)
            df['ISOLATE'] = isolate_name
            tabl_list.append(df)
        depth_df = pd.concat( tabl_list, axis = 0)
        print(depth_df)
        # merge
        breaks_df_select = breaks_df[['CHROM', 'POS']].drop_duplicates()
        depth_select = pd.merge(breaks_df_select, depth_df, how = 'left', on = ['CHROM', 'POS'])
        breaks_depth_df = pd.merge(depth_select, breaks_df, how = 'left', on = ['ISOLATE', 'CHROM', 'POS'])
        out = breaks_depth_df

    out.to_csv( name + '.breaks.tab', sep = '\t', index = None)
    


 #    merged.insert( 1, 'studyid', name )
 #    merged['N_ISOLATES']=merged.groupby(['CHROM','POS'])['ISOLATE'].transform('size')
 #    intersect = merged[ merged['N_ISOLATES'] < n_isolates ]
 #    if name == '.':
 #        name = 'split'
 #    intersect.to_csv( name + '.breaks.intersect.tab', sep = '\t', index = None )
 #    print name + '.breaks.intersect.tab saved.'

 #    # extract positions. 
 #    pos = intersect.reset_index()[['CHROM','POS']].values.astype(str).tolist()
 # #   print pos
 #    pos_for_search = []
 #    with open( name + '.pos.txt', 'w' ) as fout:
 #        for p in pos:
 #            chrom = p[0]
 #            pos = p[1]
 #            new_p = chrom + '\s+' + pos + '\n'
 #            fout.writelines( new_p )
 #    print name + '.pos.txt saved.'


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--outdir', required = True, 
                        help = 'Output directory')
    parser.add_argument('-n', '--name',
                        help = 'Pair name (default: same as output directory)')
    parser.add_argument('-i', '--input_file', required = True,
                        help = 'Tab-separated file ' 
                        'with isolate, splitters.bam filepath')
    parser.add_argument('-c', '--internal_reference', required = True,
                        help = 'Id of the internal reference')
    parser.add_argument('-r', '--reference', required = True,
                        help = 'Reference in Genbank format')
    parser.add_argument('-m', '--mask', 
                        help = 'BED file with regions to mask')
    parser.add_argument( '-d', '--min_depth',
                         help = 'Minimun breakpoint depth '
                         '(default 5)' )
    parser.add_argument( '-f', '--min_ratio',
                         help = 'Minimum ratio breakpoint depth /'
                         ' splitters depth (default 0.5)' )
    args = parser.parse_args()
    if not args.name:
        name = args.outdir
    else:
        name = args.name
    input_file = os.path.abspath(args.input_file)
    refpath = os.path.abspath(args.reference)
    internal_ref = args.internal_reference
    if args.min_depth:
        min_depth = int(args.min_depth)
    else:
        min_depth = 5
    if args.min_ratio:
        min_ratio = float(args.min_ratio)
    else:
        min_ratio = 0.5
    if args.mask:
        bed_file=os.path.abspath(args.mask)
    else:
        bed_file = None
    print('Outdir is: ' + args.outdir)
    print('Input file is: ' + input_file)
    print('Reference is: ' + refpath)
    print('Internal reference is {}'.format(internal_ref))
    print('Minimum breakpoint depth is: ' + str( min_depth ))
    print('Minimum ratio breakpoint depth / splitters depth is: ' + str( min_ratio ))
    input_dict = {}
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    os.chdir(args.outdir)
    parse_input_file( input_dict, input_file )
    if bed_file:
        make_directories( input_file, input_dict, refpath, bed_file )
    else:
        make_directories( input_file, input_dict, refpath )
    for isolate_name in input_dict:
        os.chdir( isolate_name )
        splitters_bam = input_dict [ isolate_name ]
        shutil.copy( splitters_bam, isolate_name + '.splitters.bam' )
        if args.mask:
            run_mpileup( isolate_name, mask_regions = True )
        else:
            run_mpileup( isolate_name )
#        write_allpos( isolate_name )
        draw_graph( isolate_name )
#        draw_depth_graph ( isolate_name )
        write_table ( isolate_name )
        filter_breaks(isolate_name, min_depth)
        get_insertion_pairs(isolate_name)
        os.chdir( '..' )
    for isolate_name in input_dict:
        mask_internal_reference_positions(isolate_name, internal_ref)

'''
Next steps: 1) extract start and stops at mindepth 2) use bedtools windows -r 10 -l 0 to find insertions 3) substract insertions of the internal reference
'''

# get_depth ( isolate_name, input_dict )
# filter_table( isolate_name, min_depth, min_ratio )
# intersect_tables( input_dict, name, add_splitters_depth = True )

if __name__ == "__main__":
    main(sys.argv[1:])

        
        
