########################### Amplicon analysis - Avi Levy Lab #######################
####################################################################################
#################################### Tal Dahan - Meir ##############################
####################################################################################
########################### Last updated - April 15th 2018 #########################
####################################################################################


import pandas
import glob, os
import re


path = '/home/labs/alevy/tald/'
path2 = "/home/labs/alevy/tald/"
outputdir='/home/labs/alevy/tald/'

myFastqs = pandas.read_excel(path+"libraries.xlsx")
fnames=list();

for columns,j in myFastqs.iterrows():
    filename=list()
    currentSerial=j['Serial']
    currentlib=j['library']
    currentBarcode=j['Barcode']
    gene=j['gene']
    currentexperiment=j['experiment']
    currenttimepoint=j['timepoint']
    currentgenotype=j['genotype']
#    currentTemp=j['Temperature']
#    currentVolt=j['Voltage']
#    currentHarvest=j['Harvest_time (h)']
#    currentPulse=j['Pulse']
#    currentCas=j['Cas9']



    treatID=str(currentlib)+'_'+'Barcode_'+currentBarcode+'_'+str(gene)+'_'+str(currentexperiment)+'_'+str(currenttimepoint)
    print (treatID)

    for filename in glob.glob(os.path.join(path,'*.fastq')):
        print(path2+str(currentlib))
        if filename.startswith(path2+str(currentlib)+"_"):
            print("nice!")
            print(filename)
            fnames.append(filename)
    print ("working on these files: "+str(fnames))
    newname=outputdir+str(currentSerial)+"_"+treatID+".txt"
    f= open(newname, 'w')
    count = 0
    f.write("treatID\treadID\ttile\tX_coor\tY_coor\tread_type\trbarcode\tbarcode_30Qual\tgene\tstate\tmut_class\tmut_len\tfull_seq\n")


    if gene == 'CRTISO':
        primerF = 'ATCTATAAAAGACAGC'
        dsb = 'GCATTCTGGG'
        WT_regexp='AACCCAGGAT([A,C,T,G]*?)TAAATTAAGA'
        WT_regexp_length=92;
        print (gene)
        print (currentlib)



    if gene == 'PSY1':
        primerF = 'TCGCCCCTGAATCAAAG'
        dsb = 'TGCTGCTTTG'
        WT_regexp='GCAACAACAGAG([A,C,T,G]*?)AAGAGTAAGT'
        WT_regexp_length=92;



    for i in fnames:
        name=i
        fhand = open(name)

        for line1 in fhand:
            read_identifiers=line1
            read_seq=fhand.next()
            plusLine=fhand.next()
            read_qual=fhand.next()
            #print(read_seq)
            read_tile=read_identifiers.split(":")[4]
            read_x=read_identifiers.split(":")[5]
            read_y_and_orientation=read_identifiers.split(":")[6]
            read_y=read_y_and_orientation.split()[0]
            read_ori=read_y_and_orientation.split()[1]

            markerQ_score = read_qual[:4]
            qual_4bp=0;
            for i in markerQ_score:
                if ord(i)>=63:
                    qual_4bp+=1

            if read_seq.startswith(currentBarcode):

                if read_seq.find(primerF) != -1:
                    if not read_seq.find(dsb) != -1:
                        count = count + 1

                        mut_search=re.search(WT_regexp,read_seq)

                        if mut_search:
                            mut_fragment=mut_search.group(0)
                            indel_size=len(mut_fragment)-WT_regexp_length
                            if indel_size==0:
                                mut_class="SNP"
                            elif indel_size>=1:
                                mut_class="insertion"
                            elif indel_size<=-1:
                                mut_class="deletion"
                        else:
                            mut_class="unDef"
                            indel_size="none"


                        f.write(treatID+"\t"+str(count)+"\t"+str(read_tile)+"\t"+str(read_x)+"\t"+str(read_y)+"\t"+str(read_ori)+"\t"+currentBarcode+"\t"+str(qual_4bp)+"\t"+gene+"\tmut\t"+mut_class+"\t"+str(indel_size)+"\t"+read_seq)

                    else:
                        count+=1
                        f.write(treatID+"\t"+str(count)+"\t"+str(read_tile)+"\t"+str(read_x)+"\t"+str(read_y)+"\t"+str(read_ori)+"\t"+currentBarcode+"\t"+str(qual_4bp)+"\t"+gene+"\tWT\t"+"none"+"\t"+"none"+"\t"+read_seq)


    #print (fnames)
    fnames=list();



########################### Amplicon analysis - Avi Levy Lab #######################
####################################################################################
#################################### Tal Dahan - Meir ##############################
####################################################################################
########################### Last updated - April 15th 2018 #########################
####################################################################################
