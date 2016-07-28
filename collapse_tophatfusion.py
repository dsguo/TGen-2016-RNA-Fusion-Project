#collapses tophatfusion output file, <filename>.results.txt

#import os
import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt

#read in text file in pandas to form csv table, adding header names
thpd = pd.read_csv("U138MG_ATCC.thFusion.result.txt", sep="\t", index_col=0, names=["Sample_name", "Left_gene", "Left_chr", "Left_pos_center", "Right_gene", "Right_chr", "Right_pos_center", "Spanning_sum", "Spanning_mate_pairs_sum", "End_spanning_fusion_sum", "Fusion_score" ])

#read in text file in numpy
#thnp = np.loadtxt("U138MG_ATCC.thFusion.result.txt", dtype='str')

#add header names to txt file
#thnp_names = np.savetxt("Headers_named", thnp, '''fmt= ", "''',  header="Sample_name, Left_gene, Left_chr, Left_coor, Right_gene, Right_chr, Right_coor, Spanning, Spanning_mate_pairs, End_spanning_fusion, Fusion_score")

#transpose panda table and save as another variable
thpd_trans = thpd.T

## set global variables, will be updated in needToCollapse function\
#nameL = []
# workTable = []
#tupleL = [] # contains tuples of each gene pair called in original table
#duplicatesL = [] # contains list of gene pairs that are duplicated and how many times they appear in the original table
#masterL = [] # contains list of gene pairs that are duplicated and which row indices they are in the original table
#posAverageL = [] # contains list of duplicated gene pairs as tuples, their averaged values as tuples (left, right), and lists of all their left and right values (two lists)
#readSumL = []
# shortTable = []

# Add unique number column and set index of table to numbers to avoid duplicate index problems
def addIndexL(table):
    '''Takes as input a pandas table. Adds a column of increasing whole numbers to the front of the table.'''
    nameL = [] # list of integers the length of the number of rows in the table
    for i in range(len(table.index)): # makes the python-numbered list (starts with 0) nameL
        nameL.append(i)
    global nameL
    table.insert(0, "Fusion_num", nameL) # adds the nameL list as a column into the table, at the front
    return

def setIndexL(table):
    '''Takes as input a pandas table that has a column of whole numbers added. Sets this column of numbers as the index of the pandas table, replacing the default. Returns an updated working table.'''
    workTable = table.set_index("Fusion_num") # sets the index as the newly added column Fusion_num, replacing the default index of duplicate Sample_names
    #global workTable
    return workTable


#see if there are duplicates
def needToCollapse(table):
    '''Takes as input a working table (with reset numbered index). Determines whether the table needs to be collapsed by detecting duplicates. Returns Boolean True if duplicates are detected, False if there are not duplicates and the table does not need to be collapsed.'''
    tupleL = [] # list of fusions, consisting of tuples of the gene pairs
    for row in table.index: # creates tupleL by looping through the rows and calling the gene names into tuples
        left_gene=table.at[row, "Left_gene"] # column Left_gene
        right_gene=table.at[row, "Right_gene"] # column Right_gene
        pair = (left_gene, right_gene)
        tupleL.append(pair)
    #print("tupleL:", tupleL)
    global tupleL # updates global(script) variable
    duplicatesL = [] # list of gene-pair duplicates and how many times they are duplicated
    for i in range(len(tupleL)): # for index of tuple in tupleL
        mypair = tupleL[i]
        mypaircount = 0
        for pairs in tupleL:
            if mypair == pairs:
                mypaircount += 1
        if mypaircount > 1:
            duplicate = (mypair, mypaircount)
            if duplicate not in duplicatesL: # removes repeats in duplicatesL
                duplicatesL.append(duplicate)
    global duplicatesL # updates global(script) variable
    if len(duplicatesL)>0:
        #status = "Yes, duplicates detected."
        #print(status)
        return True
    else:
        #status = "No, no duplicates detected."
        #print(status)
        return False


# identify rows with duplicates
def duplicateRowIndex(duplicatesL, tupleL):
    '''Must run needToCollapse function first. Takes as input two of needToCollapse's variables made global, duplicatesL and tupleL. Returns a master list of tuples, the duplicated gene pair tuple and the indices in the table where they occur. '''
    masterL = [] # list with duplicated gene pairs and their row indexes in the table
    for i in range(len(duplicatesL)):
        duppair = duplicatesL[i][0] # the gene pair tuple from duplicatesL
        indexL = [] # will store the gene pair's row index values in the table
        for j in range(len(tupleL)): # for row index in original table
            if duppair == tupleL[j]:
                indexL.append(j)
        duple = (duppair, indexL)
        masterL.append(duple)
        #global masterL
    return masterL

# find average of fusion position (left and right)
def averagePosition(table, masterL):
    '''Must run duplicateRowIndex first. Takes as input a pandas working table, and duplicateRowIndex output masterL. Returns a list, posAverageL, of lists containing gene pair name, average left and right positions in a tuple, as well as lists of all original given positions. '''
    posAverageL = [] # list with duplicated gene pairs, their averaged positions, and list of all of their left and right positions
    for i in range(len(masterL)): # for index in master list of duplicated gene pairs and their row indices
        pair = masterL[i][0] # gene pair as tuple
        indexL = masterL[i][1] # list of row indices
        leftValueL = [] # list of left gene positions
        rightValueL = [] # list of right gene positions
        for j in indexL:
            left_value = table.at[table.index[j], "Left_pos_center"] # position of left gene on a particular row
            right_value = table.at[table.index[j], "Right_pos_center"] # position of right gene on a particular row
            leftValueL.append(left_value)
            rightValueL.append(right_value)
        leftValueAvg = int(round(sum(leftValueL) / len(leftValueL))) # averages the left gene positions and rounds to nearest integer
        rightValueAvg = int(round(sum(rightValueL) / len(rightValueL))) # averages the right gene positions and rounds to nearest integer
        averageValues = (leftValueAvg, rightValueAvg) # tuple of average left and right gene positions
        averageInfo = [pair, averageValues, leftValueL, rightValueL] # list for one gene pair of the gene pair names, average positions, and lists of all the left and right posisions
        posAverageL.append(averageInfo)
    #global posAverageL
    return posAverageL

# find sums of supporting read counts
def sumReadCounts(table, masterL):
    '''Must run duplicateRowIndex first. Takes as input a pandas working table, and duplicateRowIndex output masterL. Returns a list, readSumL, of lists containing gene pair name, sum of supporting read counts in a list, as well as lists of all original given read counts.'''
    readSumL = [] # list with duplicated gene pairs, their sums of the three read counts as a list of length 3, and lists of all of their individual counts
    for i in range(len(masterL)): # for index in master list of duplicated gene pairs and their row indices
        pair = masterL[i][0] # gene pair as tuple
        indexL = masterL[i][1]  # list of row indices
        spanningValueL = []
        spanningMateValueL = []
        endSpanningValueL = []
        for j in indexL:
            spanning_value = table.at[table.index[j], "Spanning_sum"]
            spanningMate_value = table.at[table.index[j], "Spanning_mate_pairs_sum"]
            endSpanning_value = table.at[table.index[j], "End_spanning_fusion_sum"]
            spanningValueL.append(spanning_value)
            spanningMateValueL.append(spanningMate_value)
            endSpanningValueL.append(endSpanning_value)
        spanningSum = sum(spanningValueL)
        spanningMateSum = sum(spanningMateValueL)
        endSpanningSum = sum(endSpanningValueL)
        readSums = [spanningSum, spanningMateSum, endSpanningSum]
        sumInfo = [pair, readSums, spanningValueL, spanningMateValueL, endSpanningValueL]
        readSumL.append(sumInfo)
    #global readSumL
    return readSumL

# make new table with duplicates removed
# table.drop(table.index[list_of_rows_to_remove]
def removeDuplicateRows(table, masterL):
    '''Must run duplicateRowIndex first. Takes as input a pandas working table and duplicateRowIndex output masterL. Removes all but one of the duplicate rows for each duplicated gene pair. Returns a new shortened table which does not contain any average or summed values.'''
    removeRowL = []
    for i in range(len(masterL)):
        pair = masterL[i][0]  # gene pair as tuple
        indexL = masterL[i][1] # list of row indices
        newIndexL = indexL[1:] # list of row indices to remove (keeping the first one)
        removeRowL.extend(newIndexL)
    shortTable = table.drop(table.index[removeRowL])
    #global shortTable
    return shortTable

# make lists of revised column names for insertion, and replace duplicate values for position and read count with correct averages and sums

def columnLists(table, posAverageL, readSumL):
    '''Must run removeDuplicateRows, averagePosition, and sumReadCounts first. Takes as input a pandas shortened table, and outputs from averagePosition and sumReadCounts. Updates position and read count information for previously duplicated rows. Returns lists of values necessary for the insertion of new columns into the final pandas table.'''
    MleftPosArrayL = []
    MrightPosArrayL = []
    MspanningArrayL = []
    MspanningMateArrayL = []
    MendSpanningArrayL = []
    dupGenePairL = []
    for i in range(len(posAverageL)): # get list of gene pair tuples
        dupGenePair = posAverageL[i][0]
        dupGenePairL.append(dupGenePair)
    for row in table.index:
        leftPosArrayL = []
        rightPosArrayL = []
        spanningArrayL = []
        spanningMateArrayL = []
        endSpanningArrayL = []
        if (table.at[row, "Left_gene"], table.at[row, "Right_gene"]) not in dupGenePairL: # if the gene pair is not duplicated
            leftPosArrayL.append(table.at[row, "Left_pos_center"])
            rightPosArrayL.append(table.at[row, "Right_pos_center"])
            spanningArrayL.append(table.at[row, "Spanning_sum"])
            spanningMateArrayL.append(table.at[row, "Spanning_mate_pairs_sum"])
            endSpanningArrayL.append(table.at[row, "End_spanning_fusion_sum"])
            MleftPosArrayL.append(leftPosArrayL) # append short list to master lists that will become column values
            MrightPosArrayL.append(rightPosArrayL)
            MspanningArrayL.append(spanningArrayL)
            MspanningMateArrayL.append(spanningMateArrayL)
            MendSpanningArrayL.append(endSpanningArrayL)
        else:
            j = dupGenePairL.index((table.at[row, "Left_gene"], table.at[row, "Right_gene"])) # finds the index of this gene pair in teh dupGenePair list
            left_pos_center = posAverageL[j][1][0]
            right_pos_center = posAverageL[j][1][1]
            spanning_sum = readSumL[j][1][0]
            spanning_mate_sum = readSumL[j][1][1]
            end_spanning_sum = readSumL[j][1][2]
            table.set_value(row, "Left_pos_center", left_pos_center) # replace with average and summed values
            table.set_value(row, "Right_pos_center", right_pos_center)
            table.set_value(row, "Spanning_sum", spanning_sum)
            table.set_value(row, "Spanning_mate_pairs_sum", spanning_mate_sum)
            table.set_value(row, "End_spanning_fusion_sum", end_spanning_sum)
            left_pos_array = posAverageL[j][2]
            right_pos_array = posAverageL[j][3]
            spanning_array = readSumL[j][2]
            spanning_mate_array = readSumL[j][3]
            end_spanning_array = readSumL[j][4]
            leftPosArrayL.extend(left_pos_array) # assign individual arrays
            rightPosArrayL.extend(right_pos_array)
            spanningArrayL.extend(spanning_array)
            spanningMateArrayL.extend(spanning_mate_array)
            endSpanningArrayL.extend(end_spanning_array)
            MleftPosArrayL.append(leftPosArrayL) # append array list to master lists that will become column values
            MrightPosArrayL.append(rightPosArrayL)
            MspanningArrayL.append(spanningArrayL)
            MspanningMateArrayL.append(spanningMateArrayL)
            MendSpanningArrayL.append(endSpanningArrayL)
    #global MleftPosArrayL, MrightPosArrayL, MspanningArrayL, MspanningMateArrayL, MendSpanningArrayL
    return MleftPosArrayL, MrightPosArrayL, MspanningArrayL, MspanningMateArrayL, MendSpanningArrayL

# add necessary additional columns to pandas table
# thpd.insert(newColumnIndex, column_Name, list_of_values)

def insertNewColumns(table, MleftPosArrayL, MrightPosArrayL, MspanningArrayL, MspanningMateArrayL, MendSpanningArrayL):
    '''Must run columnLists and all other funcitons first. Takes as input a pandas shortened table and the column value lists from columnLists. Inserts and names new columns starting from the end of the table for ease of positioning. Returns a revised, updated table with appropriate columns and column values.'''
    # add from reverse to not mess up order
    table.insert(9, "End_spanning_fusion_array", MendSpanningArrayL)
    table.insert(8, "Spanning_mate_pairs_array", MspanningMateArrayL)
    table.insert(7, "Spanning_array", MspanningArrayL)
    table.insert(6, "Right_pos_array", MrightPosArrayL)
    table.insert(3, "Left_pos_array", MleftPosArrayL)
    return table

