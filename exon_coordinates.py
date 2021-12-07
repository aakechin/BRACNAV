# This script gets coordinates of BRCA1 and BRCA2 exons

def getExonCoordinates(version='hg19'):
    exons={'hg19':
           ['41196312..41197819','41199660..41199720',
           '41201138..41201211','41203080..41203134','41209069..41209152',
           '41215350..41215390','41215891..41215968','41219625..41219712',
           '41222945..41223255','41226348..41226538','41228505..41228631',
           '41234421..41234592','41242961..41243049','41243452..41246877',
           '41247863..41247939','41249261..41249306','41251792..41251897',
           '41256139..41256278','41256885..41256973','41258473..41258550',
           '41267743..41267796','41276034..41276132','41277288..41277500'],
           'hg38':
           ['40131959..40132171','40133327..40133425',
            '40141663..40141716','40150909..40150986','40152486..40152574',
            '40153181..40153320','40157562..40157667','40160153..40160198',
            '40161520..40161596','40162582..40166007','40166410..40166498',
            '40174867..40175038','40180828..40180954','40182921..40183111',
            '40186204..40186514','40189747..40189834','40193491..40193568',
            '40194069..40194109','40200307..40200390','40206325..40206379',
            '40208248..40208321','40209739..40209799','40211640..40213147']}
    brca1_exons=[]
    for ex in exons[version]:
        start,end=ex.split('..')
        brca1_exons.append([int(start),int(end)])
    if version=='hg19':
        brca1_exons=brca1_exons[::-1]
    posToExon1={}
    for i,exon in enumerate(brca1_exons):
        for pos in range(exon[0],exon[1]+1):
            posToExon1[pos]=i+1
    exons={'hg19':
           ['32889617..32889804','32890559..32890664',
           '32893214..32893462','32899213..32899321','32900238..32900287',
           '32900379..32900419','32900636..32900750','32903580..32903629',
           '32905056..32905167','32906409..32907524','32910402..32915333',
           '32918695..32918790','32920964..32921033','32928998..32929425',
           '32930565..32930746','32931879..32932066','32936660..32936830',
           '32937316..32937670','32944539..32944694','32945093..32945237',
           '32950807..32950928','32953454..32953652','32953887..32954050',
           '32954144..32954282','32968826..32969070','32971035..32971181',
           '32972299..32973809'],
           'hg38':
           ['32315480..32315667','32316422..32316527',
            '32319077..32319325','32325076..32325184','32326101..32326150',
            '32326242..32326282','32326499..32326613','32329443..32329492',
            '32330919..32331030','32332272..32333387','32336265..32341196',
            '32344558..32344653','32346827..32346896','32354861..32355288',
            '32356428..32356609','32357742..32357929','32362523..32362693',
            '32363179..32363533','32370402..32370557','32370956..32371100',
            '32376670..32376791','32379317..32379515','32379750..32379913',
            '32380007..32380145','32394689..32394933','32396898..32397044',
            '32398162..32399672']}
    brca2_exons=[]
    for ex in exons[version]:
        start,end=ex.split('..')
        brca2_exons.append([int(start),int(end)])
    posToExon2={}
    for i,exon in enumerate(brca2_exons):
        for pos in range(exon[0],exon[1]+1):
            posToExon2[pos]=i+1
    return(brca1_exons,brca2_exons,posToExon1,posToExon2)
