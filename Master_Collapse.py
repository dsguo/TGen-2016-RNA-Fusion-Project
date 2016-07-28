import collapse_tophatfusion as cthf
import pandas as pd

#inputFile
#headerNameList = ["Sample_name", "Left_gene", "Left_chr", "Left_pos_center", "Right_gene", "Right_chr", "Right_pos_center", "Spanning_sum", "Spanning_mate_pairs_sum", "End_spanning_fusion_sum", "Fusion_score"]

def masterCollapse(inputFile, headerNameList):
    '''Takes as input a TopHat-Fuions output .txt file and a list of header names (of length 11) and makes a table. Looks for duplicate gene fusions in the table rows. If found, duplicates will be collapsed into row, their positions averaged and supporting reads summed. Returns the final modified table, with has extra columns for lists of collapsed values.'''
    startTable = pd.read_csv(inputFile, sep="\t", index_col=0, names=headerNameList) # reads in .txt file and converts to a pandas table
    cthf.addIndexL(startTable) # adds a numbered column to front of table
    workTable = cthf.setIndexL(startTable) # sets the index of the working table to the numbered column, replacing default index
    if cthf.needToCollapse(workTable) == False: # if there are no duplicate gene pairs in the table
        posAverageL = [] # empty list needed for columnLists function
        readSumL = [] # empty list needed for columnLists function
        MleftPosArrayL, MrightPosArrayL, MspanningArrayL, MspanningMateArrayL, MendSpanningArrayL = cthf.columnLists(workTable, posAverageL, readSumL) # create value lists for new columns
        revisedTable = cthf.insertNewColumns(workTable, MleftPosArrayL, MrightPosArrayL, MspanningArrayL, MspanningMateArrayL, MendSpanningArrayL) # insert new columns
        print("This table does not contain duplicates. Exporting a header-inclusive .csv file.")
        return revisedTable
    else:
        masterL = cthf.duplicateRowIndex(cthf.duplicatesL, cthf.tupleL) # creates list of duplicate gene pairs and their row indices within the table
        posAverageL = cthf.averagePosition(workTable, masterL) # creates list of information on average positions and corresponding lists of original positions called
        readSumL = cthf.sumReadCounts(workTable, masterL) # creates list of information on summed read counts and corresponding lists of original read counts called
        shortTable = cthf.removeDuplicateRows(workTable, masterL) # creates shortened table with duplicate rows removed
        MleftPosArrayL, MrightPosArrayL, MspanningArrayL, MspanningMateArrayL, MendSpanningArrayL = cthf.columnLists(shortTable, posAverageL, readSumL) # creates value lists for new columns
        revisedTable = cthf.insertNewColumns(shortTable, MleftPosArrayL, MrightPosArrayL, MspanningArrayL, MspanningMateArrayL, MendSpanningArrayL) # insert new columns
        print("Duplicates detected. Exporting collapsed, header-inclusive .csv file.")
        return revisedTable

def exportAsCSV(inputFile, headerNameList):
    '''Takes as input a TopHat-Fuions output .txt file and a list of header names (of length 11) and makes a table. Runs the function masterCollapse and exports final, edited / collapsed table as a csv with file name {QC}.thFusion.result.collapsed.csv.'''
    finalTable = masterCollapse(inputFile, headerNameList)
    finalTable.to_csv(inputFile[:-4] + ".collapsed.csv", sep="\t")