# This script outputs results for BRACNAV

import xlsxwriter as xls

def createOutputFile(outFile,
                     colNames,
                     amplNumUserSortedToAmplName,
                     amplNameToAmplNumSortedByCoord,
                     amplToExon,
                     amplToIntron):
    # Create file for output
    wb=xls.Workbook(outFile)
    ws=wb.add_worksheet('CNVs')
    f0=wb.add_format({'font_name':'Times New Roman',
                      'bold':True,
                      'font_size':12,
                      'align':'center'})
    f1=wb.add_format({'font_name':'Times New Roman',
                      'bold':False,
                      'font_size':12,
                      'align':'left'})
    # colNames starts from sampleNum
    # from 1 - numbers of amplicons
    newColNamesOrderDict={}
    for i,colName in enumerate(colNames):
        if i==0:
            ws.set_column(i,i,10)
            newColNamesOrderDict[0]=colName
        elif 1<=i<=len(colNames):
            ws.set_column(i,i,15)
            # amplToExon and amplToIntron contain amplicon numbers
            # of sorted amplicons by coordinate and strand
            # So we need to convert it from user-defined to sorted
            # amplNumUserSortedToAmplName's keys start from 0
            amplName=amplNumUserSortedToAmplName[i-1]
            if amplName not in amplNameToAmplNumSortedByCoord.keys():
                continue
            sortedAmplNum=amplNameToAmplNumSortedByCoord[amplName]
            newColNamesOrderDict[sortedAmplNum+1]=colName
            if sortedAmplNum in amplToExon.keys():
                if 100<amplToExon[sortedAmplNum]<200:
                    ws.write(0,sortedAmplNum+1,
                             'brca1_exon'+str(amplToExon[sortedAmplNum]-100),f0)
                elif 200<amplToExon[sortedAmplNum]:
                    ws.write(0,sortedAmplNum+1,
                             'brca2_exon'+str(amplToExon[sortedAmplNum]-200),f0)
            elif sortedAmplNum in amplToIntron.keys():
                if 100<amplToIntron[sortedAmplNum]<200:
                    ws.write(0,sortedAmplNum+1,
                             'brca1_exon'+str(amplToIntron[sortedAmplNum]-100),f0)
                elif 200<amplToIntron[sortedAmplNum]:
                    ws.write(0,sortedAmplNum+1,
                             'brca2_exon'+str(amplToIntron[sortedAmplNum]-200),f0)
    # Set width of columns with results
    for i in range(8):
        if i<2:
            ws.set_column(len(colNames)+i,len(colNames)+i,15)
        else:
            ws.set_column(len(colNames)+i,len(colNames)+i,30)
    newColNamesOrdered=[]
    for amplNum,amplName in sorted(newColNamesOrderDict.items()):
        newColNamesOrdered.append(amplName)
    # Write column names for results
    ws.write_row(1,0,[*newColNamesOrdered,'Patient_Name','Median_Coverage',
                      'Variance_of_Values','Deleted_Exons#',
                      'Deleted_Exons_Scores','Deleted_Exons_P_values',
                      'Duplicated_Exons#','Duplicated_Exons_Scores',
                      'Duplicated_Exons_P_values'],
                 f0)
    return(wb,ws,f0,f1)
