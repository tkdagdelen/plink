from tempfile import mkstemp
from shutil import move
from os import remove, close
import csv

def getIIDsAndWrite(recodeADfilename, phenofilename):
	fh, output_file = mkstemp()
	with open(recodeADfilename) as ids_toread:
		with open(phenofilename) as phenos_toread:
			with open(output_file, "wb") as tmp_file:
				idreader = csv.reader(ids_toread, delimiter = "\t")
				phenoreader = csv.reader(phenos_read, delimiter = "\t")
				writer = csv.writer(tmp_file, delimiter = " ")
				isfirst = True
				idDict = set{}
				for idrow in idreader:     # read one row at a time
					if (len(idrow) >= 2) and (not idrow[0].startswith('#')):
						if isFirst: 
							isFirst = False
						else:
							idDict[idrow[0]]=idrow[1]
				for phenorow in phenoreader:
					if (len(phenorow) >= 2) and (not phenorow[0].startswith('#')):
						if isFirst: 
							myColumn = ['FID', 'IID'] + list(phenorow[i].replace(" ","") for i in range(2,len(phenorow)))
							writer.writerow(myColumn)
							isfirst = False
						else:
							myColumn = [phenorow[1], idDict[phenorow[1]]]  + list(phenorow[i] for i in range(2,len(row)))
							writer.writerow(myColumn) # write it
	close(fh)
	move(output_file,"fixed_with_iids_"+filename)

phenofilename = input("Pheno filename:")
recodeADfilename = input("RecodeAD filename:")
getIIDsAndWrite(recodeADfilename, phenofilename)
