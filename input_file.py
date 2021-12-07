# This script reads input file for BRACNAV

import os
import statistics as stat
import numpy as np
from copy import deepcopy
from scipy.cluster.hierarchy import linkage

# Our modules
from clusterization import showDendrogram,getClusters

# Reads input file with coverages
# As an input it takes (inputFile):
# - file name with coverage values with the following columns:
#   + patient number
#   + amplicon 1 coverages
#   + amplicon 2 coverages
#   + etc...
# - output file name (outFile, for saving file with sample dendrogram)
# - dictionary that converts multiplexes 
# - minimal acceptable median coverage for sample (minCov)
def readInputFile(inputFile,outFile,
                  multToAmpls,amplToMults,
                  amplNumUserSortedToAmplName,
                  amplNumCoordSortedToAmplName,
                  amplNameToAmplNumSortedByCoord,
                  amplNameToAmplNum,
                  minCov,notClust):
    # Read input file with coverage values
    file=open(inputFile)
    # Stores column names
    colNames=[]
    # Stores sample numbers
    patNames=[]
    # Stores matrix of coverage values
    data=[]
    # Stores, if the input file contains patient info
    withPatIDs=False
    for string in file:
        cols=string.replace('\n','').split('\t')
        # If it is the 1st string
        if ('Patient#' in string or
            'Patient_Num' in string):
            # Save column names
            if 'Patient_ID' in cols:
                colNames=cols[0:1]+cols[5:]
                withPatIDs=True
            else:
                colNames=cols[:]
            continue
        # Save sample name
        patNames.append(cols[0])
        # Save all coverage values into the matrix
        if withPatIDs:
            data.append(list(map(float,cols[5:])))
        else:
            data.append(list(map(float,cols[1:])))
    file.close()
    # Convert matrix as list to numpy array
    data0=np.array(data)
    # If some coverage values are less than 1,
    ## increase it by 1. It's necessary for not dividing by 0
    data0[data0<1]+=1
    # Create new array that will store normalized values
    data1=np.array(data0)
    # Save number of samples (maxi) and amplicons (maxj)
    maxi,maxj=data0.shape
    # Stores sample numbers that have normal median coverage
    # (not less than one that was defined by user
    normCovPatients=[]
    # Stores median coverages for each sample
    medianCovs=[]
    for i in range(maxi):
        try:
            medianCovs.append(round(stat.median(data0[i]),3))
        except stat.StatisticsError:
            print('ERROR (3)! No data in the following row:')
            print('Row number:',i)
            print(data0)
            exit(3)
        # If median coverage for this sample is not less than user defined
        if stat.median(data0[i])>=minCov:
            # We save it to list of samples with normal coverage
            normCovPatients.append(i)
        # Go through all amplicons
        for j in range(maxj):
            # We change order of values to order of exons
            # We need to convert amplicon numbers from user-defined to sorted
            # amplNumUserSortedToAmplName's keys start from 0
            try:
                amplName=amplNumUserSortedToAmplName[j]
            except KeyError:
                print('ERROR (1)! Amplicons file does not correcpond to coverage file by amplicon set')
                print('Amplicons have the following IDs:')
                print(amplNumUserSortedToAmplName.keys())
                print('But the following was searched:')
                print(j)
                exit(1)
            if amplName in amplNameToAmplNumSortedByCoord.keys():
                sortedAmplNum=amplNameToAmplNumSortedByCoord[amplName]
                # Normalize the value by the median value of coverage
                # for all amplicons of this sample
                data1[i,sortedAmplNum]=data0[i,j]/stat.median(data0[i])
    
    # Cluster all samples to normalize them in each cluster
    # We determine output directory name to save dendrogram plot later
    dirName=outFile[:outFile.rfind('/')+1]
    outFilePart=outFile[outFile.rfind('/')+1:-4]
    # Perform clusterization
    z=linkage(data1,'ward')
    # Create directory for outputing plots
    if not os.path.exists(dirName+outFilePart+'_plots'):
        os.mkdir(dirName+outFilePart+'_plots')
    # Draw dendrogram with all samples
    showDendrogram(z,''.join([dirName,
                              outFilePart,
                              '_plots/',
                              outFilePart,
                              '_sample_clusters.png']))
    clusters0=getClusters(z)
    clusters=deepcopy(clusters0)
    patToCluster={}
    for key,item in clusters0.items():
        for it in item:
            patToCluster[it]=key
            if it not in normCovPatients:
                clusters[key].remove(it)
    # Get numbers of amplicons that has low coverage. Further we will exclude them from the analysis
    lowCovAmplicons=[]
    for j in range(maxj):
        # We need to convert amplicon numbers from user-defined to sorted
        # amplNumUserSortedToAmplName's keys start from 0
        amplName=amplNumUserSortedToAmplName[j]
        if amplName in amplNameToAmplNumSortedByCoord.keys():
            sortedAmplNum=amplNameToAmplNumSortedByCoord[amplName]
            if stat.median(data0[:,j])<minCov:
                lowCovAmplicons.append(sortedAmplNum)
    data2=np.array(data1)
    # Normalize by number of reads for the amplicon
    for i in range(maxi):
        if i in normCovPatients:
            # If user chose not to cluster samples
            if notClust:
                for j in range(maxj):
                    data2[i,j]=data1[i,j]/stat.median(data1[:,j])
            else:
                clusterKeys=clusters[patToCluster[i]]
                for j in range(maxj):
                    data2[i,j]=data1[i,j]/stat.median(data1[clusterKeys,j])
    data3=np.array(data2)
    # Values of deviation from 2
    devs=np.array(data2)
    # Normalize by number of reads for the multiplexes
    for i in range(maxi):
        if i in normCovPatients:
            multMedianValues=[]
            for k in range(len(multToAmpls)):
                # Collect coverage values for all amplicons
                # of current multiplex (primer pool)
                multValues=[]
                for m in multToAmpls[k]:
                    # m is a user-sorted number of amplicon
                    # We need to convert amplicon numbers from user-defined to sorted
                    # amplNumUserSortedToAmplName's keys start from 0
                    try:
                        amplName=amplNumUserSortedToAmplName[int(m)]
                    except KeyError:
                        print('ERROR (2)! Amplicons file does not correcpond to coverage file by amplicon set')
                        print('Amplicons have the following IDs:')
                        print(amplNumUserSortedToAmplName.keys())
                        print('But the following was searched:')
                        print(m,type(m),
                              m in amplNumUserSortedToAmplName.keys(),
                              int(m) in amplNumUserSortedToAmplName.keys(),
                              str(m) in amplNumUserSortedToAmplName.keys())
                        exit(2)
                    if amplName in amplNameToAmplNumSortedByCoord.keys():
                        sortedAmplNum=amplNameToAmplNumSortedByCoord[amplName]
                        multValues.append(float(data2[i,sortedAmplNum]))
                # Calculate median value for current multiplex
                multMedianValues.append(stat.median(multValues))
            for j in range(maxj):
                # We need to convert amplicon numbers
                # from sorted by coordinate to
                # sorted by user in the amplicons file
                # because multMedianValues contains such numbers
                if j in amplNumCoordSortedToAmplName.keys():
                    amplName=amplNumCoordSortedToAmplName[j]
                    userSortedAmplNum=amplNameToAmplNum[amplName]
                    data3[i,j]=float(data2[i,j])*2/multMedianValues[amplToMults[userSortedAmplNum]]
##                devs[i,j]=abs(data3[i,j]-2)
    data3[data3>4]=4
    # Get percentile values for each amplicon for plotting
    minYs=[]
    maxYs=[]
    for j in range(maxj):
        perc1=np.percentile(data3[:,j],20)
        perc2=np.percentile(data3[:,j],80)
        minYs.append(perc1)
        maxYs.append(perc2)

    # Calculate differences between current amplicon
    # and previous one
    allValVars=[]
    for i in range(maxi):
        valVars=[]
        for j in range(maxj)[1:]:
            valVars.append(abs(data3[i,j]-data3[i,j-1]))
        allValVars.append(stat.median(valVars))
    return(data3,normCovPatients,
           medianCovs,lowCovAmplicons,
           maxi,maxj,colNames,patNames,
           dirName,outFilePart,
           minYs,maxYs,
           allValVars)


