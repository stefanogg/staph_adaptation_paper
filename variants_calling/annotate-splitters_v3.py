#!/usr/bin/env python3

'''
Script to summarise and annotate intervals with split reads that are unique to the isolate inside the pair.
Intervals are extracted from the list of unique splitter breakpoints generated by the split-break pipeline 
and verified manually.

This script needs improvements:
- fix the bug that generates duplicate intervals
- remove the time-consuming and often failing assembly step

'''

'''
We could add a break_start and break_stop attribute,
obtained with a method that uses pysam. Extact breaks with highest
count, remove reads with the break, extract breaks with highest count 
of the opposite type if count is >= min_count (e.g. 10),
In a next step, we use the breaks attributes and number of supp_intervals
        ( 1 or > 1 ) to tentatively define the splitter interval as 
        deletion type (left or right) or insertion type.
'''
'''
        This approach could be used also for the initial filtering of splits
        instead of split break. Define SplitterInterval class from the intervalList
        in pybedtools, filter, intervals with reads >= min_count, obtain break 1 and 
        2 (if present, a splitter interval can not have more than two breaks), 
        filter intervals with break count >= min_count (add other paramters,
        frac and especially frac as compared to total read coverage, as splits 
        may represent only a minority of reads). 
        Once we obtain a list of potential splitter intervals per isolate, 
        use R to filter breaks (not intervals) as done previously.
        
        Another idea would be to redefine the splitter intervals by selectin bam 
        files based on the breaks. The intervals could be smaller, more specific and
        we would exclude peripjeral splitters that are anyway difficult to interpret.
        
        we would have initial (provisional) splitter intervals and definitive 
        splitter intervals.
'''
'''
                An other idea: identify "true breaks" as opposed to normal start/end of read
'''


import sys
import os
import shutil
import argparse
import re
import pandas as pd
import pysam
import pybedtools
from pybedtools import BedTool
from BCBio import GFF

# Classes

class SplitterInterval( object ):
    class_counter = 1
    def __init__( self, type, isolate_name, intervalList = [], ):
        self.isolate = isolate_name
        self.chrom = intervalList[0]
        self.start = intervalList[1]
        self.stop = intervalList[2]
        # self.count = intervalList[3]
        self.length = str(int(self.stop) - int(self.start))
        # self.breakpoint = intervalList[5]
        self.type = type
        self.supp_intervals = []
        self.bamfiles = []
        self.id = isolate_name + '-SPLITTER' + str( SplitterInterval.class_counter )
        SplitterInterval.class_counter += 1

        
    
    def write_bed( self ):
        SplitterInterval = [ self.chrom, self.start, self.stop, self.id, 
            self.length, self.type]
        filename = self.id + '.bed'
        with open( filename, 'w') as fout:
            new_line = '\t'.join( SplitterInterval ) + '\n'
            fout.writelines( new_line )
        return filename
    
    def get_BAM( self, splitters_path ):
        '''
        Get splitter BAM files for the interval and save them in BAM format
        '''
        bedfile = self.id + '.bed'
        bamfile = self.id + '.bam'
        cmd = 'samtools view -b -L ' + bedfile + ' ' + \
        splitters_path + ' > ' + bamfile
        print('Running: ' + cmd) 
        os.system( cmd )
        self.bamfiles.append( bamfile )

        
    def get_suppAlignment( self ):
        '''
        Get intervals with supplementary alignment
        '''
        bampath = self.id + '.bam'
        suppAlignmentList = []
        with pysam.AlignmentFile( bampath, 'rb') as samfile:
            for i in samfile.fetch( until_eof = True ):
                sa_tag = i.get_tag( 'SA' )
                suppAlignmentList.append( sa_tag )
        # print 'samtools view ' + self.id + '.bam'
#         os.system('samtools view ' + self.id + '.bam')
        # wrap bedtools
        # a file
        a_filename = self.isolate + '.splitters.bed'
        a = pybedtools.BedTool( a_filename )
        # b file
        supp_posList = []
        for sa_tag in suppAlignmentList:
            columns = sa_tag.split(',')
            chrom = columns[0]
            start = columns[1]
            stop = columns[1]
            supp_pos = '\t'.join( [ chrom, start, stop ] ) + '\n'
            if supp_pos not in supp_posList:
                supp_posList.append( supp_pos )
        supp_pos_string = ''.join(supp_posList) 
        b = pybedtools.BedTool( supp_pos_string, from_string = True )
        # intersect
        a_and_b = a.intersect(b, wa = True)
        supp_intervals_list = []
        for i in a_and_b:
            if i not in supp_intervals_list:
                supp_intervals_list.append(i)
        self.supp_intervals = supp_intervals_list
        return supp_intervals_list
        
    def write_sa_bed( self ):
        '''
        Generates bed file for supplementary alignements
        '''
        filename = self.id + '.sa.bed'
        with open( filename, 'w') as fout:
            for i in self.supp_intervals:
                count = i[3]
                length = i[4]
                type = 'supplementary'
                id = self.id
                SplitterInterval = [ i.chrom, i.start, i.stop, id, length, type ]
                new_line = '\t'.join( map( str, SplitterInterval ) ) + '\n'
                fout.writelines( new_line )
        return filename
        
    def get_suppBAM ( self, splitters_path ):
        '''
        Get splitter BAM files for the intervals and save them in BAM format
        '''
        bedfile = self.id + '.sa.bed'
        bamfile = self.id + '.sa.bam'
        cmd = 'samtools view -b -L ' + bedfile + ' ' + \
            splitters_path + ' > ' + bamfile
        print('Running: ' + cmd)
        os.system( cmd )
        self.bamfiles.append( bamfile )
        return self.bamfiles
        
class annotatedInterval(object):
    def __init__( self, intervalList = [] ):
        self.chrom = intervalList[0]
        self.start = intervalList[1]
        self.stop = intervalList[2]
        self.id = intervalList[3]
        # self.count = intervalList[4]
        self.length = intervalList[4]
        # self.breakpoint = intervalList[6]
        self.type = intervalList[5]
        self.annotation = []

    def annotateInterval( self, reference ):
        # write interval to BedTool object
        interval = [ self.chrom, self.start, self.stop ]
        interval_string = '\t'.join( interval ) 
                # print('Interval string is {interval_string}'.format(interval_string = interval_string))
        a = reference
                # print('Argument -a is {a}'.format(a = reference))
        b = pybedtools.BedTool( interval_string, from_string = True )
                # print(b)
        gff_tmp = '/tmp/giulieris/gff.tmp'
                # print(gff_tmp)
        intervals_annotated = a.intersect( b, wa = True ).saveas( gff_tmp )
                # print('Annotated intervals are {i}'.format(i = intervals_annotated))
        # parse GFF
        with open( gff_tmp, 'r' ) as f:
            for rec in GFF.parse( f ):
                for feat in rec.features:
                    if feat.type == 'CDS':
                        self.annotation.append( feat )
                # print(self.annotation)
        annotated_entries = []
        common_variables = [self.chrom, self.start, self.stop, self.id,
                    self.length, self.type]
        for feat in self.annotation:



            type = feat.type
            product = feat.qualifiers['product'][0]
            locus_tag = feat.qualifiers['ID'][0]
            location = feat.location
            annotated_interval = common_variables + [type, product, locus_tag, str(location) ]
            annotated_entries.append( annotated_interval )
        return annotated_entries
        
        
# Functions

def createSplitterIntervals( isolate_name, split_dir ):
        splitters_path = split_dir + '/' + isolate_name + '/' + isolate_name + '.splitters.bam'
        filenames_list = []
        bamfiles_list = []
        bedfile = isolate_name + '.insertion.pairs.bed'
        bamfiles_file = isolate_name + '.bamfiles.txt'
        with open( bedfile, 'r') as f, open( bamfiles_file, 'w') as fout:
                for line in f:
                        intervalList = line.split()
                        interval = SplitterInterval( 'primary', isolate_name, intervalList )
                        filename = interval.write_bed()
                        filenames_list.append( filename )
                        interval.get_BAM( splitters_path )
                        interval.get_suppAlignment()
                        filename = interval.write_sa_bed()
                        filenames_list.append( filename )
                        bamfiles = interval.get_suppBAM( splitters_path )
                        for f in bamfiles:
                                bamfiles_list.append( f )
                        for f in bamfiles_list:
                                new_line = f + '\n'
                                fout.writelines( new_line )
        print('Splitter intervals created for {i}'.format(i = isolate_name))
        return filenames_list
            
def summarise_intervals( isolate_name, filenames_list ):
    # generate a list of intervals and write summary bed file 
        # intervals_summary = []
    summary_filename = isolate_name + '.primary.sa.bed'
    with open( summary_filename, 'w') as fout:
        for line in [open(f).read() for f in filenames_list]:
#             intervals_summary.append( line )
            fout.writelines( line )
    # annotate bed file with ref.gff
    reference = pybedtools.BedTool( '../reference/ref.gff' )
        # print('Reference in bedtools format is {r}'. format(r = reference))
    # list of intervals to annotate
    intervals_summary = pybedtools.BedTool( summary_filename )
    # annotate and write to file
    annotated_filename = isolate_name + '.annotated.bed'
    with open( annotated_filename, 'w' ) as fout:
        for i in intervals_summary:
                        intervalList = i.fields
                        annotated = annotatedInterval( intervalList)
                        annotated_entries = annotated.annotateInterval( reference )
                        if not annotated_entries:
                                new_line = '\t'.join( intervalList ) + '\n'
                                fout.writelines( new_line )
                        else:
                                for i in annotated_entries:
                                        new_line = '\t'.join( i ) + '\n'
                                        fout.writelines( new_line )

    # get fasta sequences. 
    intervals_newId = []
    for i in intervals_summary:
            i.name = i.name + '_' + i.fields[5]
            intervals_newId.append(i)
    intervals_newId = pybedtools.BedTool(intervals_newId)
    seq = intervals_newId.sequence( fi = '../reference/ref.fa', name = True )
    seq_file = isolate_name + '.seq.fa'
    with open( seq_file, 'w') as fout:
        fout.write( open( seq.seqfn).read())
    # run cd_hit
    out_file = isolate_name + '.out'
    cmd = 'cd-hit-est -d 0 -i ' + seq_file + ' -o ' + \
        out_file 
    print('Running: ' + cmd)
    os.system( cmd )
    # parse cd-hit
    clusterfile = out_file + '.clstr'
    print('Parsing cd-hit output ...')
    parse_outputfile = out_file + '.tab'
    cmd='clstr2txt.pl ' + clusterfile + ' > ' + parse_outputfile
    os.system( cmd )
        
    
def get_assembly( isolate_name ):
    # merge bam files and write to fastq format
    bamfiles_file = isolate_name + '.bamfiles.txt'
    out_bam = isolate_name + '.confirmed.bam'
    out_fq = isolate_name + '.confirmed.fastq'
    merge_cmd = 'samtools merge -f -b ' + bamfiles_file + ' ' + \
        out_bam
    bam2fq_cmd = 'samtools bam2fq ' + out_bam + ' > ' + \
        out_fq
    for cmd in [ merge_cmd, bam2fq_cmd ]:
        print('Running: ' + cmd)
        os.system( cmd )
    # assemble with spades
        # is it really necessary to assemble all the reads?
    if not os.path.exists( 'spades' ):
        os.mkdir( 'spades' )
    cmd = 'spades.py -o spades/all-confirmed -s ' + out_fq
    print('Running: ' + cmd)
    os.system( cmd )
    # splitters specific assemblies    
    # obtain the list of unique *SPLITTER* for the isolate
    summary_filename = isolate_name + '.primary.sa.bed'
    splitter_id_list = []
    with open( summary_filename, 'r') as f:
        for line in f:
            id = line.split()[3]
            if id not in splitter_id_list:
                splitter_id_list.append( id )
    # generate tmp file with *SPLITTER* bamfilenames and merge bam according to *SPLITTER*
    # run spades
    bamfiles_file_tmp = '/tmp/sgiulieri/bamfiles.tmp'
    for id in splitter_id_list:
        out_bam = id + '.bam'
        out_fq = id + '.fastq'
        with open( bamfiles_file, 'r' ) as f, open( bamfiles_file_tmp, 'w') as fout:
            for line in f:
                    if id in line:
                        fout.writelines( line )
        merge_cmd = 'samtools merge -f -b ' + bamfiles_file_tmp + ' ' + \
            out_bam
        bam2fq_cmd = 'samtools bam2fq ' + out_bam + ' > ' + \
            out_fq
        spades_cmd = 'spades.py -o spades/' + id + ' -s ' + out_fq
        for cmd in [ merge_cmd, bam2fq_cmd, spades_cmd ]:
            print('Running: ' + cmd)
            os.system( cmd )
        
            
        
    
            
def main():
    
    usage = """%prog -o <outdir> -i <file> -d <split_dir>
annotate-splitters.py
Author: Stefano Giulieri
Description: Annotate unique splitters intervals 
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir',
        help = 'Output directory')
    parser.add_argument('-i', '--input_file',
        help= 'Tab-separated file with chrom, pos, start, stop, isolate',
        metavar= 'FILE')
    parser.add_argument( '-d', '--split_dir',
        help= 'Directory with splitter bam files')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
        outdir = args.outdir
        input_file = os.path.abspath(args.input_file)
        split_dir = os.path.abspath(args.split_dir)
    os.chdir(outdir)
    isolates_list = []
    with open(input_file) as f:
        for line in f:
            iso=line.split()[3]
            isolates_list.append(iso)
    isolates_list=list(set(isolates_list))
    for iso in isolates_list:
        os.chdir(iso)
        filenames_list = createSplitterIntervals( isolate_name = iso, split_dir = split_dir )
        summarise_intervals( isolate_name = iso, filenames_list = filenames_list )
        # get_assembly( i )
        os.chdir('..')
            
if __name__ == "__main__":
    sys.exit(main()) 