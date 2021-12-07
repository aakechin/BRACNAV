# This script reads file with amplicon coordinates
# and determines which amplicons match to which exons

import os

thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

def matchAmpliconsToExonsAndIntrons(amplicons,posToExon,prefix=100,
                                    amplToExon={},exonToAmpls={},
                                    amplToIntron={},intronToAmpls={},
                                    amplNameToAmplNumSortedByCoord={},
                                    amplNumCoordSortedToAmplName={}):
    # Amplicon number
    if len(amplToExon)==0 and len(amplToIntron)==0:
        i=0
    elif len(amplToExon)>0 and len(amplToIntron)>0:
        i=max(max(amplToExon.keys()),max(amplToIntron.keys()))+1
    elif len(amplToExon)>0:
        i=max(amplToExon.keys())+1
    else:
        i=max(amplToIntron.keys())+1
    # If gene is located at plus strand
    if posToExon[min(posToExon.keys())]<posToExon[max(posToExon.keys())]:
        reversedOrder=False
    # If gene is located at minus strand
    else:
        reversedOrder=True
    amplNums=[]
    allExonPoses=list(posToExon.keys())
    for amplName,coordinates in sorted(amplicons.items(),
                                       key=lambda item:item[1],
                                       reverse=reversedOrder):
        thisAmpliconOverlapsExon=False
##        print(amplName,coordinates)
        for pos in range(coordinates[1],coordinates[2]+1):
            if pos in posToExon.keys():
                thisAmpliconOverlapsExon=True
                amplToExon[i]=prefix+posToExon[pos]
                if prefix+posToExon[pos] not in exonToAmpls.keys():
                    exonToAmpls[prefix+posToExon[pos]]=[i]
                else:
                    exonToAmpls[prefix+posToExon[pos]].append(i)
                break
##        print(amplToExon)
        if not thisAmpliconOverlapsExon:
            allExonPosesPlusThis=sorted(allExonPoses+[coordinates[1]])
            thisNum=allExonPosesPlusThis.index(coordinates[1])
            if thisNum==0 or thisNum==len(allExonPosesPlusThis)-1:
                intronNum=0
            else:
                # If current gene is on plus strand
                if posToExon[min(posToExon.keys())]<posToExon[max(posToExon.keys())]:
                    intronNum=posToExon[allExonPosesPlusThis[thisNum-1]]
                else:
                    intronNum=posToExon[allExonPosesPlusThis[thisNum+1]]
            amplToIntron[i]=prefix+intronNum
            if intronNum not in intronToAmpls.keys():
                intronToAmpls[prefix+intronNum]=[i]
            else:
                intronToAmpls[prefix+intronNum].append(i)
##        print(amplToIntron)
        amplNums.append(i)
##        print(amplNums)
        amplNameToAmplNumSortedByCoord[amplName]=i
        amplNumCoordSortedToAmplName[i]=amplName
##        print(amplNameToAmplNumSortedByCoord)
##        print(amplNumCoordSortedToAmplName)
        i+=1
    return(amplToExon,exonToAmpls,
           amplToIntron,intronToAmpls,
           amplNums,amplNameToAmplNumSortedByCoord,
           amplNumCoordSortedToAmplName)

def getAmpliconCoordinates(ampliconFile,posToExon1,posToExon2):
    if ampliconFile:
        file=open(ampliconFile)
    else:
        file=thisDir+'default_amplicons.csv'
    # Converts number of amplicon to number of multiplex
    amplToMults={}
    # Converts number of multiplex to number of amplicon
    multToAmpls={}
    # Converts number of amplicon to number of exon
    amplToExon={}
    # Converts number of amplicon to number of intron
    amplToIntron={}
    # Contains numbers of amplicons for each exon
    exonToAmpls={0:[]}
    # Contains numbers of amplicons for some introns
    intronToAmpls={}
    # Converts name of amplicon to number of amplicon among
    # amplicons sorted by user
    amplNameToAmplNum={}
    # Converts number of amplicon among all amplicons
    # sorted by user to name of amplicon
    amplNumUserSortedToAmplName={}
    # Converts number of amplicon 
    # sorted by coordinate to name of amplicon
    amplNumCoordSortedToAmplName={}
    # BRCA1 amplicons defined by user
    brca1_amplicons={}
    # BRCA2 amplicons defined by user
    brca2_amplicons={}
    # Other amplicons
    other_amplicons={}
    # Amplicon number among amplicons sorted by user (from 0)
    i=0
    for string in file:
        cols=string.replace('\n','').split('\t')
        if len(cols)>4:
            amplToMults[i]=int(cols[4])-1
            if int(cols[4])-1 not in multToAmpls.keys():
                multToAmpls[int(cols[4])-1]=[i]
            else:
                multToAmpls[int(cols[4])-1].append(i)
        else:
            amplToMults[i]=0
            if 0 not in multToAmpls.keys():
                multToAmpls[0]=[i]
            else:
                multToAmpls[0].append(i)
        if cols[1]=='17' or cols[1]=='chr17':
            brca1_amplicons[cols[0]]=[cols[1],int(cols[2]),int(cols[3])]
        elif cols[1]=='13' or cols[1]=='chr13':
            brca2_amplicons[cols[0]]=[cols[1],int(cols[2]),int(cols[3])]
        else:
            other_amplicons[cols[0]]=[cols[1],int(cols[2]),int(cols[3])]
##            print('ERROR! Unknown chromosome number in the amplicon file:\n',cols[1])
##            exit(0)
        if cols[0] in amplNameToAmplNum.keys():
            print('ERROR! The following amplicon name is repeated in the amplicons file:')
            print(cols[0])
            exit(1)
        amplNameToAmplNum[cols[0]]=i
        amplNumUserSortedToAmplName[i]=cols[0]
        i+=1
    file.close()
    # Amplicons order. It it necessary for drawing plot, because we need to sort amplicons
    # in the following order: BRCA1 exons 1-23, BRCA2 exons 1-27
    ampliconsOrder=[]
    amplExonIntronMatches=matchAmpliconsToExonsAndIntrons(brca1_amplicons,posToExon1,
                                                          prefix=100,amplToExon={},exonToAmpls={},
                                                          amplToIntron={},intronToAmpls={},
                                                          amplNameToAmplNumSortedByCoord={},
                                                          amplNumCoordSortedToAmplName={})
    amplToExon,exonToAmpls,amplToIntron,intronToAmpls=amplExonIntronMatches[:4]
    amplNums,amplNameToAmplNumSortedByCoord=amplExonIntronMatches[4:6]
    amplNumCoordSortedToAmplName=amplExonIntronMatches[6]
##    print(amplNums)
##    # If gene is located at plus strand
##    if posToExon1[min(posToExon1.keys())]<posToExon1[max(posToExon1.keys())]:
    ampliconsOrder.extend(sorted(amplNums,reverse=False))
##    else:
##        ampliconsOrder.extend(sorted(amplNums,reverse=True))
##    print(ampliconsOrder)
##    input()
    amplExonIntronMatches=matchAmpliconsToExonsAndIntrons(brca2_amplicons,
                                                          posToExon2,
                                                          200,
                                                          amplToExon,
                                                          exonToAmpls,
                                                          amplToIntron,
                                                          intronToAmpls,
                                                          amplNameToAmplNumSortedByCoord,
                                                          amplNumCoordSortedToAmplName)
    amplToExon,exonToAmpls,amplToIntron,intronToAmpls=amplExonIntronMatches[:4]
    amplNums,amplNameToAmplNumSortedByCoord=amplExonIntronMatches[4:6]
    amplNumCoordSortedToAmplName=amplExonIntronMatches[6]
##    # If gene is located at plus strand
##    if posToExon2[min(posToExon2.keys())]<posToExon2[max(posToExon2.keys())]:
    ampliconsOrder.extend(sorted(amplNums,reverse=False))
##    else:
##        ampliconsOrder.extend(sorted(amplNums,reverse=True))
    return(amplToExon,
           exonToAmpls,
           amplToIntron,
           intronToAmpls,
           ampliconsOrder,
           multToAmpls,
           amplToMults,
           amplNumUserSortedToAmplName,
           amplNumCoordSortedToAmplName,
           amplNameToAmplNumSortedByCoord,
           amplNameToAmplNum,
           len(brca1_amplicons))
##            # We have 3 variants of covering exon by an amplicon:
##            # ----[------------]----
##            #    <-------------->
##            # <------->  or <------->
##            #         <---->
##            # Go through all BRCA1 exons
##            for j,exon in enumerate(brca1_exons):
##                # If current amplicon's start is less than start of current exon
##                # and amplicon's end is more than end of current exon
##                # and it is not the last exon and current amplicon's end is less than start of the next one
##                if int(cols[2])<=exon[0]+EXON_COVERAGE_BUFFER and int(cols[3])>=exon[0]-EXON_COVERAGE_BUFFER and (j<len(brca1_exons)-1 and int(cols[3])<=brca1_exons[j+1][0]):
##                    amplToExon[i]=prefix+j+1
##                    if prefix+j+1 not in exonToAmpls.keys():
##                        exonToAmpls[prefix+j+1]=[i]
##                    else:
##                        exonToAmpls[prefix+j+1].append(i)
##                    break
##                # If current amplicon's end is more than end of current exon minus buffer
##                # and the amplicon's start is less than start of current exon plus buffer
##                # and it is not the last exon and current amplicon's end is less than start of the next one
##                if int(cols[3])>=exon[1]-EXON_COVERAGE_BUFFER and int(cols[2])<=exon[1]+EXON_COVERAGE_BUFFER and (j<len(brca1_exons)-1 and int(cols[3])<=brca1_exons[j+1][0]):
##                    amplToExon[i]=prefix+j+1
##                    if prefix+j+1 not in exonToAmpls.keys():
##                        exonToAmpls[prefix+j+1]=[i]
##                    elif i not in exonToAmpls[prefix+j+1]:
##                        exonToAmpls[prefix+j+1].append(i)
##                    break
##                # If current amplicon's start is more than start of current exon minus buffer
##                # and the amplicon's end is less than end of current exon plus buffer
##                if int(cols[2])>=exon[0]-EXON_COVERAGE_BUFFER and int(cols[3])<=exon[1]+EXON_COVERAGE_BUFFER:
##                    amplToExon[i]=prefix+j+1
##                    if prefix+j+1 not in exonToAmpls.keys():
##                        exonToAmpls[prefix+j+1]=[i]
##                    else:
##                        exonToAmpls[prefix+j+1].append(i)
##                    break
##                # Else we save it as intron
##                amplToExon[i]=0
##                exonToAmpls[0].append(i)
##                amplToIntron
