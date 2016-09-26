from tempfile import mkstemp
from shutil import move
from os import remove, close
import csv


""" Takes the phenotype file given by dbGaP and changes all of the garnet accession numbers to 1"""
def fixPhenoFile(filename):
        fh, output_file = mkstemp()
        with open(filename) as to_read:
                with open(output_file, "wb") as tmp_file:
                        reader = csv.reader(to_read, delimiter = "\t")
                        writer = csv.writer(tmp_file, delimiter = " ")
                        isTrue = False
                        isfirst = True
                        for row in reader:     # read one row at a time
                                if len(row) >= 2:
                                        if not row[0].startswith('#'):#'dbGaP SubjID':
                                                isTrue = True
                                if isTrue and not isfirst:
                                        myColumn = [row[1], row[1]] + list(row[i] for i in range(2,len(row)))
                                        writer.writerow(myColumn) # write it
                                if isfirst and isTrue:
                                        myColumn = ['FID', 'IID'] + list(row[i].replace(" ","") for i in range(2,len(row)))
                                        writer.writerow(myColumn)
                                        isfirst = False

        close(fh)
        move(output_file,"fixed_"+filename)

name = input("filename:")
fixPhenoFile(name)