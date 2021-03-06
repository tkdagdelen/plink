from tempfile import mkstemp
from shutil import move
from os import remove, close
import csv
import sys

def getIIDsAndWrite(imissFilename, phenofilename):
	fh, output_file = mkstemp()
	with open(imissFilename) as ids_toread:
		with open(phenofilename) as phenos_toread:
			with open(output_file, "wb") as tmp_file:
				csv.field_size_limit(sys.maxsize)
				idreader = csv.reader(ids_toread, delimiter = " ")
				phenoreader = csv.reader(phenos_toread, delimiter = "\t")
				writer = csv.writer(tmp_file, delimiter = "\t")
				isfirst = True
				idDict = {}
				for idrow in idreader:     # read one row at a time
					idrow = filter(None,idrow)
					if (len(idrow) >= 2) and (not idrow[0].startswith('#')):
						if isfirst: 
							isfirst = False
						else:
							idDict[idrow[0]]=idrow[1]
				isfirst = True
				for phenorow in phenoreader:
					if (len(phenorow) >= 2) and (not phenorow[0].startswith('#')):
						if isfirst: 
							myColumn = ['FID', 'IID'] + list(phenorow[i].replace(" ","") for i in range(2,len(phenorow)))
							writer.writerow(myColumn)
							isfirst = False
						else:
							try:
							 	genderStatus = [0]
							 	if (phenorow[5]=="M"): 
							 		genderStatus = [1]
							 	for p in range(0, len(phenorow)): 
							 		if phenorow[p]=="Yes":
							 			phenorow[p] = 1
							 		if phenorow[p] == "No": 
							 			phenorow[p] = 0
								myColumn = [phenorow[1], idDict[phenorow[1]]]  + list(phenorow[i] for i in range(2,5)) + genderStatus + list(phenorow[j] for j in range(6,len(phenorow)))
							except KeyError: 
								# do nothing
								do=None
							writer.writerow(myColumn) # write it
	close(fh)
	move(output_file,"fixed_with_iids_"+phenofilename)

phenofilename = input("Pheno filename:")
imissFilename = input("imiss filename:")
getIIDsAndWrite(imissFilename, phenofilename)
