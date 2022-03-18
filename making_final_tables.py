import glob
import os
import csv
import  subprocess

dir=['01/', '02/', '03/', '04/', '05/', '06/', '07/', '08/', '09/', '10/', '10_12_17/', '11/', '12/', '13/', '14/', '15/', '16/', '17/', '18/', '19/', '20/', '21/', '22/', '23/', '24/', '25/', '26/', '27/']

results_consesus=open("merge_consensus.txt","w")
results_reference=open("merge_reference.txt","w")

for file in dir:
	directory=file.replace("/","")
	os.chdir(file)
	os.chdir("reference/")
	read_csv = list(csv.reader(open("covrin_reference_"+directory+".txt"), delimiter="\t"))
	line=directory+"\t".join(read_csv[-1])
	line+="\n"
	results_reference.write(line)
	os.chdir("../..")

for file in dir:
	directory=file.replace("/","")
	os.chdir(file)
	os.chdir("consensus/")
	read_csv = list(csv.reader(open("covrin_consensus_"+directory+".txt"), delimiter="\t"))
	line=directory+"\t".join(read_csv[-1])
	line+="\n"
	results_consesus.write(line)
	os.chdir("../..")




