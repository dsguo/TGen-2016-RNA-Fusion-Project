# Usage
# python Process_Discordant_SAM_File.py <Collapsed_TophatFusion_File> <Discordant_Reads.sam>

# requires python 3.4 with pandas, numpy, pysam packages

#Configure Enviroment
import pandas as pd
import numpy as np
import pysam
import sys


# Function to process SAM file from SAMBLASTER DISCORDANT EXPORT
def bam_to_df(bam, chr = None, start=None, stop = None):
    seq = []
    name = []
    r1_chr = []
    r1_pos = []
    r1_isRead1 = []
    r1_isReversed = []
    r1_cigar = []
    r1_mapq = []
    r2_isReversed = []
    r2_chr = []
    r2_pos = []
    frag_length = []
    r2_cigar = []
    r2_mapq = []
    for read in bam.fetch(chr, start, stop):
        seq.append(read.query_sequence)
        name.append(read.query_name)
        r1_chr.append(read.reference_name)
        r1_pos.append(read.reference_start)
        r1_isRead1.append(read.is_read1)
        r1_isReversed.append(read.is_reverse)
        r1_cigar.append(read.cigarstring)
        r1_mapq.append(read.mapping_quality)
        r2_isReversed.append(read.mate_is_reverse)
        r2_chr.append(read.next_reference_name)
        r2_pos.append(read.next_reference_start)
        frag_length.append(read.template_length)
        if read.has_tag('MC') == True:
            r2_cigar.append(read.get_tag('MC'))
        else:
            r2_cigar.append('NA')
        if read.has_tag('MQ') == True:
            r2_mapq.append(read.get_tag('MQ'))
        else:
            r2_mapq.append(255)
    return pd.DataFrame({'seq': seq,
                         'name': name,
                         'r1_pos': r1_pos,
                         'r1_chr': r1_chr,
                         'r1_isRead1': r1_isRead1,
                         'r1_isReversed': r1_isReversed,
                         'r1_cigar': r1_cigar,
                         'r1_mapq': r1_mapq,
                         'r2_isReversed': r2_isReversed,
                         'r2_chr': r2_chr,
                         'r2_pos': r2_pos,
                         'frag_length': frag_length,
                         'r2_cigar': r2_cigar,
                         'r2_mapq': r2_mapq})


# Import Collapsed Tophat-Fusion table
tf_table_file = tf_file = sys.argv[1]
tf_table = pd.read_csv(tf_table_file, sep="\t")

# Create python object for the SAM file
disc_sam_file = sys.argv[2]
samfile = pysam.AlignmentFile(disc_sam_file, "r")

# Call Function to convert SAM file into pandas dataframe, and add additional columns so original chromosome columns have integer values
# chromosomes are already in string format
chrList = range(23)
strChrList = ["{:01d}".format(x) for x in chrList]
table = bam_to_df(samfile)
table["r1_chr_int"] = np.nan
table["r2_chr_int"]= np.nan
for row in table.index:
    if table.at[row, "r1_chr"] in strChrList or table.at[row, "r1_chr"] in range(23):
        table.set_value(row, "r1_chr_int", table.at[row, "r1_chr"])
    elif table.at[row, "r1_chr"] == "X":
        table.set_value(row, "r1_chr_int", 23)
    elif table.at[row, "r1_chr"] == "Y":
        table.set_value(row, "r1_chr_int", 24)
    elif table.at[row, "r1_chr"] == "M": # or MT or mitochondrial chromosome calls
        table.set_value(row, "r1_chr_int", 25)

for row in table.index:
    if table.at[row, "r2_chr"] in strChrList or table.at[row, "r2_chr"] in range(23):
        table.set_value(row, "r2_chr_int", table.at[row, "r2_chr"])
    elif table.at[row, "r2_chr"] == "X":
        table.set_value(row, "r2_chr_int", 23)
    elif table.at[row, "r2_chr"] == "Y":
        table.set_value(row, "r2_chr_int", 24)
    elif table.at[row, "r2_chr"] == "M": # or MT or mitochondrial chromosome calls
        table.set_value(row, "r2_chr_int", 25)

# Extract just those rows that are the first read of the pair
firstReadTable = table[table.r1_isRead1 == True]

# Add needed columns to the collapsed tophat file
# Create new columns on the unique gene pairs table
tf_table['discordantFrag_Count'] = np.nan
tf_table['discordantFrag_For_Count'] = np.nan
tf_table['discordantFrag_Rev_Count'] = np.nan

# add columns to collapsed table so original chromosome columns have integer values
tf_table["firstGene_Chr_Int"] = np.nan
tf_table["secondGene_Chr_Int"]= np.nan
for row in tf_table.index:
    if tf_table.at[row, "firstGene_Chr"] in range(23) or tf_table.at[row, "firstGene_Chr"] in strChrList:
        tf_table.set_value(row, "firstGene_Chr_Int", tf_table.at[row, "firstGene_Chr"])
    elif tf_table.at[row, "firstGene_Chr"] == "X":
        tf_table.set_value(row, "firstGene_Chr_Int", 23)
    elif tf_table.at[row, "firstGene_Chr"] == "Y":
        tf_table.set_value(row, "firstGene_Chr_Int", 24)
    elif tf_table.at[row, "firstGene_Chr"] == "M": # or MT or mitochondrial chromosome calls
        tf_table.set_value(row, "firstGene_Chr_Int", 25)

for row in tf_table.index:
    if tf_table.at[row, "secondGene_Chr"] in range(23) or tf_table.at[row, "secondGene_Chr"] in strChrList:
        tf_table.set_value(row, "secondGene_Chr_Int", tf_table.at[row, "secondGene_Chr"])
    elif tf_table.at[row, "secondGene_Chr"] == "X":
        tf_table.set_value(row, "secondGene_Chr_Int", 23)
    elif tf_table.at[row, "secondGene_Chr"] == "Y":
        tf_table.set_value(row, "secondGene_Chr_Int", 24)
    elif tf_table.at[row, "secondGene_Chr"] == "M": # or MT or mitochondrial chromosome calls
        tf_table.set_value(row, "secondGene_Chr_Int", 25)


### THIS NEEDS TO BE A LOOP OF THE IMPORTED COLLAPPSED TOPHAT FUSION TABLE - updates the collapsed table :)
# Create test variables - THES ARE PRE_CALCULATED IN THE IMPORT
for row in tf_table.index:
    window1_chr=tf_table.at[row, 'firstGene_Chr_Int']
    window1_start=tf_table.at[row, 'firstGene_Window_Start']
    window1_end=tf_table.at[row, 'firstGene_Window_End']
    window2_chr=tf_table.at[row, 'secondGene_Chr_Int']
    window2_start=tf_table.at[row, 'secondGene_Window_Start']
    window2_end=tf_table.at[row, 'secondGene_Window_End']
    # extract possible pairs
    # because there can be two derivatives you need to test for both possible orrientations
    result_table = firstReadTable[((firstReadTable.r1_chr_int == window1_chr) &
                                   (firstReadTable.r1_pos >= window1_start) &
                                   (firstReadTable.r1_pos <= window1_end) &
                                   (firstReadTable.r2_chr_int == window2_chr) &
                                   (firstReadTable.r2_pos >= window2_start) &
                                   (firstReadTable.r2_pos <= window2_end)) |
                                  ((firstReadTable.r1_chr_int == window2_chr) &
                                   (firstReadTable.r1_pos >= window2_start) &
                                   (firstReadTable.r1_pos <= window2_end) &
                                   (firstReadTable.r2_chr_int == window1_chr) &
                                   (firstReadTable.r2_pos >= window1_start) &
                                   (firstReadTable.r2_pos <= window1_end))]
    # Get table length
    discordantFrag_CountNum = len(result_table.index)
    tf_table.set_value(row, 'discordantFrag_Count', discordantFrag_CountNum) # add value to collapsed table
    # Make tables for both possible derivatives
    # Read1 aligned to forward strand
    r1_for_result_table = result_table[result_table.r1_isReversed == False]
    r1_for_count = len(r1_for_result_table.index)
    tf_table.set_value(row, 'discordantFrag_For_Count', r1_for_count) # add value to collapsed table
    # Read1 aligned to reverse strand
    r1_rev_result_table = result_table[result_table.r1_isReversed == True]
    r1_rev_count = len(r1_rev_result_table.index)
    tf_table.set_value(row, 'discordantFrag_Rev_Count', r1_rev_count) # add value to collapsed table


# Create tables to figure out breakpoint locations
# Create forward sorted tables
# r1_for_result_table.sort_values(['r1_pos'], ascending=[1])

# Create reverse sorted tables
# r1_rev_result_table.sort_values(['r1_pos'], ascending=[1])

# Write out final table
# Save output to file
tf_table.to_csv("Tophat_Fusion_Results_Collapsed_Final.txt", sep="\t", index=False, float_format='%.0f')