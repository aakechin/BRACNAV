# This script contains function for clusterization of value arrays

import numpy as np
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt

# Draws and saves dendrogram of clusterized samples
# As an input it takes:
# - result of linage function from scipy.cluster.hierarchy (z)
# - file name for output plot of dendrogram (fileName)
def showDendrogram(z,fileName):
    # calculate full dendrogram
    plt.figure(figsize=(25, 10))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    dendrogram(
        z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    )
    plt.savefig(fileName)
    plt.close()

# Split all samples into clusters
# As input it takes:
# - result of linage function from scipy.cluster.hierarchy (z)
def getClusters(z):
    # Stores size of z
    size=len(z)+1
    # Get maximal distance between branches
    maxDist=z[-1][2]+10
    clusters={}
    bigClusters={}
    bigCl=False
    zarr=np.array(z)
    zStd=np.std(zarr[:,2])
    for i in range(len(z))[::-1]:
        if z[i][2]-z[i-1][2]>=zStd*1.5:
            maxDist=z[i][2]
    for i in range(len(z)):
        if z[i][2]>=maxDist and not bigCl:
            bigCl=True
        if not bigCl: clusters[size+i]=[]
        # If all big clusters has been already formed, we do not join them
        if ((int(z[i][0]) in clusters.keys() or int(z[i][0]) in bigClusters.keys())
            and (int(z[i][1]) in clusters.keys() or int(z[i][1]) in bigClusters.keys()) and bigCl):
            bigClusters[size+i]=[int(z[i][0]),int(z[i][1])]
            continue
        if int(z[i][0]) in clusters.keys():
            if bigCl:
                clusters[size+i]=[]
            clusters[size+i].extend(clusters[int(z[i][0])])
            clusters.pop(int(z[i][0]))
        elif int(z[i][0]) in bigClusters.keys() and int(z[i][1]) not in bigClusters.keys() and int(z[i][1]) not in clusters.keys():
            clusters[size+i].extend(clusters[bigClusters[int(z[i][0])][0]])
            clusters.pop(bigClusters[int(z[i][0])][0])
        elif int(z[i][0]) not in clusters.keys():
            if bigCl:
                clusters[size+i]=[]
            clusters[size+i].append(int(z[i][0]))
        if int(z[i][1]) in clusters.keys():
            clusters[size+i].extend(clusters[int(z[i][1])])
            clusters.pop(int(z[i][1]))
        elif int(z[i][1]) in bigClusters.keys()  and int(z[i][0]) not in bigClusters.keys() and int(z[i][0]) not in clusters.keys():
            clusters[size+i].extend(clusters[bigClusters[int(z[i][1])][0]])
            clusters.pop(bigClusters[int(z[i][1])][0])
        elif int(z[i][1]) not in clusters.keys():
            clusters[size+i].append(int(z[i][1]))
    if len(clusters)==0:
        clusters[0]=list(range(size))
##    print(z)
##    print(zStd)
##    print(clusters)
    return(clusters)
