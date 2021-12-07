# This script reads information about each sample

def getSampleInfo(sampleFile):
    file=open(sampleFile)
    patIds={}
    for string in file:
        if 'Patient_num' in string:
            continue
        cols=string.replace('\n','').split('\t')
        patIds[cols[0]]=cols[1]
    file.close()
    return(patIds)
