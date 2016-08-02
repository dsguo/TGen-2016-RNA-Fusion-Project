# Usage Example
# python Process_TF_Output.py <TopHat_Fusion_Result.txt> <refGene.txt> <ensGene.txt>


#Configure Enviroment
import pandas as pd
import numpy as np
import sys

# Set a validation screening window
window = 60000

# set a distance-apart window for overlapping position windows (above)
dist = 2000

# Read in the original tophat-fusion result text file
#tf_file = "/Users/jkeats/MMRF_2000.txt"
tf_file = sys.argv[1]
# Define column header
col_headers = ['Specimen',
               'firstGene', 'firstGene_Chr', 'firstGene_Pos',
               'secondGene', 'secondGene_Chr', 'secondGene_Pos',
               'spanningReads', 'spanningMatePairs', 'spanningMPFusion', 'score']
tf_table = pd.read_csv(tf_file, sep="\t", header=None, names=col_headers)
#tf_table.head()

# Read in the refGene.txt and ensGene.txt files used for annotations by TopHat-Fusion
table_headers = ['uniq_id', 'transcript_id', 'chromosome', 'strand', 'transcript_start', 'transcript_end',
               'cds_start', 'cds_end', 'exon_count', 'exon_starts', 'exon_ends', 'blank', 'gene_id',
               'cds_start_stat', 'cds_end_stat', 'exon_frames']

refGene_file = sys.argv[2]
refGene_table = pd.read_csv(refGene_file, sep="\t", header=None, names=table_headers)
#refGene_table.head()

ensGene_file = sys.argv[3]
ensGene_table = pd.read_csv(ensGene_file, sep="\t", header=None, names=table_headers)
#ensGene_table.head()

# Read in the exon model file with ENSG and HUGO IDs
#exon_models_file = "/Users/jkeats/Ensembl_V74_ENSG_HUGO_ENSE_Table.txt"
#exon_models_file = sys.argv[2]
#exon_models_table = pd.read_csv(exon_models_file, sep="\t")
#exon_models_table.head()

# Create a list of unique first and second gene pairs
gene_pairs = tf_table.loc[:,['firstGene', 'secondGene']]
gene_pairs.drop_duplicates(inplace=True)
#gene_pairs.head()

# Create new columns on the unique gene pairs table
gene_pairs['ordered_genes'] = np.nan
gene_pairs['ordered_first'] = np.nan
gene_pairs['ordered_second'] = np.nan
gene_pairs['fusion_pairs'] = np.nan
gene_pairs['firstGene_Strand'] = np.nan
gene_pairs['firstGene_Chr'] = np.nan
gene_pairs['firstGene_Window_Start'] = np.nan
gene_pairs['firstGene_Window_End'] = np.nan
gene_pairs['firstGene_Pos_List'] = np.nan
gene_pairs['secondGene_Strand'] = np.nan
gene_pairs['secondGene_Chr'] = np.nan
gene_pairs['secondGene_Window_Start'] = np.nan
gene_pairs['secondGene_Window_End'] = np.nan
gene_pairs['secondGene_Pos_List'] = np.nan
gene_pairs['spanningReads_Sum'] = np.nan
gene_pairs['spanningReads_List'] = np.nan
gene_pairs['spanningMatePairs_Sum'] = np.nan
gene_pairs['spanningMatePairs_List'] = np.nan
gene_pairs['spanningMPFusion_Sum'] = np.nan
gene_pairs['spanningMPFusion_List'] = np.nan

# Loop through unique gene pair table to extract counts from tophat-fusion output table
for row in gene_pairs.index:
    # Create list with the two gene pairs on each row
    first = gene_pairs.at[row, 'firstGene']
    second = gene_pairs.at[row, 'secondGene']
    temp = [first, second]
    # Sort the temp list created for each row
    temp.sort()
    # Extract each of the sorted values and create a concatenation
    ordered_first = temp[0]
    ordered_second = temp[1]
    ordered_genes = ordered_first + "_" + ordered_second

    # Extract the strand for first gene
    if len(first) == 15 and first[:4] == "ENSG":
        # print(first + " is a ENSG_ID")
        first_Gene_table = ensGene_table[ensGene_table.gene_id == first]
    else:
        # print(first + " is a HUGO_ID")
        first_Gene_table = refGene_table[refGene_table.gene_id == first]
    loop = 0
    for i in first_Gene_table.index:
        if loop == 0:
            firstGene_Strand = first_Gene_table.at[i, 'strand']
            loop = 1
        elif loop == 1:
            firstGene_Strand2 = first_Gene_table.at[i, 'strand']

    # Extract the strand for second gene
    if len(second) == 15 and second[:4] == "ENSG":
        # print(first + " is a ENSG_ID")
        second_Gene_table = ensGene_table[ensGene_table.gene_id == second]
    else:
        # print(first + " is a HUGO_ID")
        second_Gene_table = refGene_table[refGene_table.gene_id == second]
    loop = 0
    for i in second_Gene_table.index:
        if loop == 0:
            secondGene_Strand = second_Gene_table.at[i, 'strand']
            loop = 1
        elif loop == 1:
            secondGene_Strand2 = second_Gene_table.at[i, 'strand']

    # Make a table from the input tophat table for each unique line
    pair_table = tf_table[(tf_table.firstGene == first) & (tf_table.secondGene == second)]
    fusion_pairs = len(pair_table.index)
    spanningReads_Sum = pair_table['spanningReads'].sum()
    spanningMatePairs_Sum = pair_table['spanningMatePairs'].sum()
    spanningMPFusion_Sum = pair_table['spanningMPFusion'].sum()

    loop = 0
    for line in pair_table.index:
        if loop == 0:
            firstGene_Chr = pair_table.at[line, 'firstGene_Chr']
            secondGene_Chr = pair_table.at[line, 'secondGene_Chr']
            firstGene_Pos_List = pair_table.at[line, 'firstGene_Pos']
            secondGene_Pos_List = pair_table.at[line, 'secondGene_Pos']
            spanningReads_List = pair_table.at[line, 'spanningReads']
            spanningMatePairs_List = pair_table.at[line, 'spanningMatePairs']
            spanningMPFusion_List = pair_table.at[line, 'spanningMPFusion']
            loop = 1
        elif loop == 1:
            firstGene_Chr2 = pair_table.at[line, 'firstGene_Chr']
            secondGene_Chr2 = pair_table.at[line, 'secondGene_Chr']
            # now test to ensure the chromosomes of each line are not changing
            if firstGene_Chr == firstGene_Chr2 and secondGene_Chr == secondGene_Chr2:
                print('Chromosomes Match')
            else:
                print('ERROR - ERROR')
            firstGene_Pos_List2 = pair_table.at[line, 'firstGene_Pos']
            firstGene_Pos_List = str(firstGene_Pos_List) + ";" + str(firstGene_Pos_List2)
            secondGene_Pos_List2 = pair_table.at[line, 'secondGene_Pos']
            secondGene_Pos_List = str(secondGene_Pos_List) + ";" + str(secondGene_Pos_List2)
            spanningReads_List2 = pair_table.at[line, 'spanningReads']
            spanningReads_List = str(spanningReads_List) + ";" + str(spanningReads_List2)
            spanningMatePairs_List2 = pair_table.at[line, 'spanningMatePairs']
            spanningMatePairs_List = str(spanningMatePairs_List) + ";" + str(spanningMatePairs_List2)
            spanningMPFusion_List2 = pair_table.at[line, 'spanningMPFusion']
            spanningMPFusion_List = str(spanningMPFusion_List) + ";" + str(spanningMPFusion_List2)

    # Determine the proper windows to target for WGS extraction/validation
    # first gene
    firstGene_Pos_Mean = pair_table['firstGene_Pos'].mean()
    firstGene_Window_Start = int(firstGene_Pos_Mean) - window
    firstGene_Window_End = int(firstGene_Pos_Mean) + window
    # second gene
    secondGene_Pos_Mean = pair_table['secondGene_Pos'].mean()
    secondGene_Window_Start = int(secondGene_Pos_Mean) - window
    secondGene_Window_End = int(secondGene_Pos_Mean) + window

    ### check if the two windows are overlapping: Not sure if other variables must be considered, such as the first gene being further along in the chromosome
    if firstGene_Chr == secondGene_Chr: # if the fusion is on the same chromosome
        if firstGene_Window_End < secondGene_Window_Start:
            pass
        elif firstGene_Window_End > secondGene_Window_Start:
            midpoint = int(round(sum([firstGene_Pos_Mean, secondGene_Pos_Mean])/2))
            firstGene_Window_End = midpoint - dist
            secondGene_Window_Start = midpoint + dist
    else: pass



    # Add the respective values to each respective row
    gene_pairs.at[[row], 'ordered_genes'] = ordered_genes
    gene_pairs.at[[row], 'ordered_first'] = ordered_first
    gene_pairs.at[[row], 'ordered_second'] = ordered_second
    gene_pairs.at[[row], 'fusion_pairs'] = fusion_pairs
    gene_pairs.at[[row], 'firstGene_Strand'] = firstGene_Strand
    gene_pairs.at[[row], 'firstGene_Chr'] = firstGene_Chr
    gene_pairs.at[[row], 'firstGene_Window_Start'] = firstGene_Window_Start
    gene_pairs.at[[row], 'firstGene_Window_End'] = firstGene_Window_End
    gene_pairs.at[[row], 'firstGene_Pos_List'] = firstGene_Pos_List
    gene_pairs.at[[row], 'secondGene_Strand'] = secondGene_Strand
    gene_pairs.at[[row], 'secondGene_Chr'] = secondGene_Chr
    gene_pairs.at[[row], 'secondGene_Window_Start'] = secondGene_Window_Start
    gene_pairs.at[[row], 'secondGene_Window_End'] = secondGene_Window_End
    gene_pairs.at[[row], 'secondGene_Pos_List'] = secondGene_Pos_List
    gene_pairs.at[[row], 'secondGene_Pos_List'] = secondGene_Pos_List
    gene_pairs.at[[row], 'spanningReads_Sum'] = spanningReads_Sum
    gene_pairs.at[[row], 'spanningReads_List'] = spanningReads_List
    gene_pairs.at[[row], 'spanningMatePairs_Sum'] = spanningMatePairs_Sum
    gene_pairs.at[[row], 'spanningMatePairs_List'] = spanningMatePairs_List
    gene_pairs.at[[row], 'spanningMPFusion_Sum'] = spanningMPFusion_Sum
    gene_pairs.at[[row], 'spanningMPFusion_List'] = spanningMPFusion_List

#gene_pairs.head()

# Save output to file
gene_pairs.to_csv("Tophat_Fusion_Results_Collapsed.txt", sep="\t", index=False, float_format='%.0f')

# Make a BED file with first and second gene Chr, Window_Start, and Window_End
firstGene_BED = gene_pairs.loc[:,['firstGene_Chr', 'firstGene_Window_Start', 'firstGene_Window_End']]
secondGene_BED = gene_pairs.loc[:,['secondGene_Chr', 'secondGene_Window_Start', 'secondGene_Window_End']]

# Concatentate and sort the two bed files
# Give each file a common header
col_headers2 = ['Chr', 'Start', 'End']
firstGene_BED.columns = col_headers2
secondGene_BED.columns = col_headers2

# Concatenate
query_BED = pd.concat([firstGene_BED, secondGene_BED], axis=0, join='outer', ignore_index=True)

# Sort the Query_BED
query_BED.sort_values(['Chr', 'Start', 'End'], ascending=[1, 1, 1], inplace=True)

# Save Output for manipulation with BEDTOOLS to collapse intervals
# Save output to file
query_BED.to_csv("Temp_Query_BED.bed", sep="\t", index=False, header=False, float_format='%.0f')