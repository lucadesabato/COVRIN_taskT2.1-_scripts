import csv
from Bio import SeqIO
import operator
import subprocess
import sys
import os
import glob
import statistics

dir=sys.argv[1].replace("/","")
reference_sequence=sys.argv[2]
os.chdir(dir)

#count the occurence of each element in a list and return the frequence of the most repeated 
def most_frequent(lists):
	counter = 0
	num = lists[0]

	for i in lists:
		curr_frequency = lists.count(i)
		if(curr_frequency> counter):
			counter = curr_frequency
			num = i
	return num

#find the Pangolin lineage of a sequence in a Pangoline results file (csv)
def find_lineage_in_Pangolin(name_sequence, pangolin):
	lineages=''
	read_csv = list(csv.reader(open(pangolin), delimiter=","))
	for line in read_csv[1:]:
		if line[0]==name_sequence:
			lineages+=line[1]
	return lineages

#choose the degeneration based on the nucleotides in a list
def find_degeneration(nucleotide):
	deg=''
	if "a" in nucleotide and "g" in nucleotide and "t" not in nucleotide:
		deg+="r"
	elif "a" in nucleotide and "g" in nucleotide and "t" in nucleotide:
		deg+="d"
	elif "a" in nucleotide and "c" in nucleotide and "t" not in nucleotide:
		deg+="m"
	elif "c" in nucleotide and "g" in nucleotide and "t" not in nucleotide:
		deg+="s"
	elif "g" in nucleotide and "t" in nucleotide:
		deg+="k"
	elif "a" in nucleotide and "t" in nucleotide:
		deg+="w"
	elif "a" in nucleotide and "c" in nucleotide and "t" in nucleotide:
		deg+="h"
	elif "c" in nucleotide and "g" in nucleotide and "t" in nucleotide:
		deg+="b"
	elif "a" in nucleotide and "c" in nucleotide and "g" in nucleotide:
		deg+="v"
	elif "c" in nucleotide and "t" in nucleotide:
		deg+="y"
	return deg

#find deletions at 5' and 3' of a sequence
def find_del_at_beg_end(seq):
	k=0
	l__begin=0
	l_end=0
	while seq[k]=="-":
		l__begin+=1
		k+=1
	j=len(seq)-1
	while seq[j]=="-":
		l_end+=1
		j-=1
	return [l__begin, l_end]

#based on alignment of sequences with a reference, with the reference in the first position, change the N of the refence
#put the nucleotide with the highest frequency

def change_reference(allineamento):
	filename=allineamento.replace(".fasta","")
	alignment = list(SeqIO.parse(allineamento, "fasta"))
	file=open(filename+"_FINE.fasta", "w")
	k=0
	reference=alignment[0].seq
	consensus=''
	while k<len(reference):
		if alignment[0].seq[k] in ["n", "N"]:
			v=[]
			for seq in alignment[1:]:
				if seq.seq[k] in ["a", "g", "c", "t", "A", "G", "C", "T"] and seq.seq[k] not in ["n","N"]:
					v.append(seq.seq[k])
			#print(k)
			#print(v)
			if v and len(v)>1:
				common=most_frequent(v)
				#print(common)
				consensus+=common
			elif v and len(v)==1:
				consensus+=v[0]
			else:
				consensus+="n"
			k+=1
		else:
			consensus+=alignment[0].seq[k]
			k+=1


	file.write(">"+str(alignment[0].id)+"\n")
	file.write(consensus+"\n")
	for seq in alignment[1:]:
		file.write(">"+str(seq.id)+"\n")
		file.write(str(seq.seq)+"\n")
	consesus_only=open(dir+"_New_consensus_only.fasta","w")
	consesus_only.write(">"+str(alignment[0].id)+"\n")
	consesus_only.write(consensus+"\n")
	return True

def count_nucleotide(seq):
	count=0
	for nucleotide in seq:
		if nucleotide in ["a", "g", "c", "t", "A", "G", "C", "T"]:
			count+=1
	return count

def count_deg(seq):
	count=0
	for nucleotide in seq:
		if nucleotide in ["R", "D", "M", "S", "K", "W", "H", "B", "V", "Y",
							  "r", "d", "m", "s", "k", "w", "h", "b", "v", "y"]:
			count+=1
	return count

def count_N(seq):
	count=0
	for nucleotide in seq:
		if nucleotide in ["n", "N"]:
			count+=1
	return count

def count_gap(seq):
	count=0
	for nucleotide in seq:
		if nucleotide=="-":
			count+=1
	return count

#return the SNP considering only AGCT bases
#return a dictionary with position of the reference, the nucleotide reference and the nucleotide of sequence
def find_SNP(reference, seq):
	k=0
	differences={}
	while k<len(reference):
		if reference[k] in ["a", "g", "c", "t", "A", "G", "C", "T"] and seq[k]!="-" and seq[k] not in ["n", "N", "R", "D", "M", "S", "K", "W", "H", "B", "V", "Y",
							  "r", "d", "m", "s", "k", "w", "h", "b", "v", "y"] and reference[k]!=seq[k]:
			differences[k]=[reference[k], seq[k]]
			k+=1
		else:
			k+=1
	return differences

#return a dictionary with the start position of the insertion  and the lenght of the insertion
def find_insertions(reference, seq):
	k=0
	ins_long={}
	lun=0
	while k<len(reference):
		if reference[k] == "-" and seq[k]!="-" and seq[k] not in ["n", "N"] and reference[k]!=seq[k]:
			lun+=1
			k+=1
		else:
			if lun>0:
				ins_long[k-lun+1]=lun
				k+=1
				lun=0
			else:
				lun=0
				k+=1
	if lun>0:
			ins_long[k-lun+1]=lun
	return ins_long

#return a dictionary with the position of the deletion, the nucleotide reference and '-' of the sequence
def find_deletions_and_count(reference, seq):
	k=0
	differences={}
	while k<len(reference):
		if reference[k] in ["a", "g", "c", "t", "A", "G", "C", "T", "R", "D", "M", "S", "K", "W", "H", "B", "V", "Y",
							  "r", "d", "m", "s", "k", "w", "h", "b", "v", "y"] and seq[k]=="-" and reference[k]!=seq[k]:
			differences[k]=[reference[k], seq[k]]
			k+=1
		else:
			k+=1
	return differences

#return a dictionary with the start position of the deletion and the lenght of the deletion
def find_deletions(seq):
	del_long={}
	k=0
	lun=0
	for n in seq:
		if n == "-":
			lun+=1
			k+=1
		else:
			if lun>0:
				del_long[k-lun+1]=lun
				k+=1
				lun=0
			else:
				lun=0
				k+=1
	if lun>0:
			del_long[k-lun+1]=lun
	return del_long

#create a fasta  file for each sequence in a fasta file
def extract_seq():
	records = list(SeqIO.parse(dir.replace("/","")+".fasta", "fasta"))
	for seq in records:
		file=open(str(seq.id)+".fasta", "w")
		file.write(">"+str(seq.id)+"\n")
		file.write(str(seq.seq))
	return True

#return minum, maximum and median value from a list of number
def extract_min_max_median(list):
	results=[]
	results.append(min(list))
	results.append(max(list))
	results.append(statistics.median(list))
	return results

def collect_results(directory):
	update=list(SeqIO.parse("../UPDATE_alignment_samples_01_to_22_ref_SARS_CoV_2_MN908947.fasta", "fasta"))
	if directory=="reference/":
		os.chdir(directory)
		vs_reference=open("covrin_reference_"+dir+".txt", "w")
		vs_reference.write("Sample\tlenght\tdegenerations\tSNP\tNNN_size\tinsertions\tdeletions\n")
		alignment=sorted(glob.glob("*_ref_ALL.fasta"))

		consensus_size = []
		NNN_size = []
		insertions_size = []
		degenerations_size = []
		lineages = []
		lenght_to_compare=[]

		deletions_samples={}
		insertions_samples={}
		for file_sequence in alignment:
			riga = file_sequence.replace("_ref_ALL.fasta", "") + "\t"
			r = list(SeqIO.parse(file_sequence, "fasta"))
			reference = r[0].seq
			to_analyze = r[1].seq
			for updated in update:
				if updated.id==r[1].id:
					lenght_seq=str(count_nucleotide(updated.seq))
					lenght_to_compare.append(int(lenght_seq))
					break
			consensus_size.append(int(lenght_seq))
			n_differences = find_SNP(reference, to_analyze)
			n_degenerations = str(count_deg(to_analyze))
			degenerations_size.append(int(n_degenerations))
			read_distance = list(csv.reader(open("../distance/"+dir+"_pdistance.txt"), delimiter="\t"))
			number_SNP=""
			for line in read_distance[2:]:
				if line[0]==file_sequence.replace("_ref_ALL.fasta", ""):
					number_SNP=line[1]
	#print("-----")
	#print(n_differences)
	#print("-----")
			number_NNN = str(count_N(to_analyze))
			NNN_size.append(int(number_NNN))
			# number_gap=str(count_gap(to_analyze))

			n_insertions = find_insertions(reference, to_analyze)
			number_ins = len(n_insertions)
			total_ins=0
			for key in n_insertions:
				total_ins+=n_insertions[key]
			insertions_samples[file_sequence]=[total_ins, number_ins]

	# print(n_insertions)

			n_deletions = find_deletions(to_analyze)
			total_del = count_gap(to_analyze) - find_del_at_beg_end(to_analyze)[0] - find_del_at_beg_end(to_analyze)[1]
			number_del = len(n_deletions)
			deletions_samples[file_sequence]=[total_del, number_del]

	#print("DEELEZIONI")
	#print(n_deletions)
	#print("del--------")
			riga += lenght_seq + "\t" + n_degenerations + "\t" + number_SNP + "\t" + number_NNN + "\t" + str(total_ins) + "\t" + str(total_del) + "\n"
			vs_reference.write(riga)
			lineages.append(find_lineage_in_Pangolin(file_sequence.replace("_ref_ALL.fasta", ""), "../../pangolin_covrin.csv"))

			#scrivo tabella dei singoli per grafici
			tab.write(file_sequence.replace("_ref_ALL.fasta", "") + "\t"+lenght_seq+"\t"+str(number_ins)+"\t"+str(total_ins) +"\t"+str(number_del)+"\t"+str(total_del)+"\t"+number_SNP+"\t"+number_NNN+"\t"+str(n_degenerations)+"\n")

		risultati="Final_results"+dir+"\tSNP_min\tSNP_max\tSNP_median\tpd_min\tpd_max\tpd_median\tCons_size_min\tCons_size_max\tCons_size_diff\tNNN_size_min\tNNN_size_max" \
				  "\tNNNsize_noMinMax\tNumb_ins\tins_min\tins_max\tins_median\tNumb_del\tdel_min\tdel_max\tdel_median\tDeg_min\tDeg_max\tlineage\n"
		vs_reference.write(risultati)
		vs_reference.write("\t")
#write snp with reference
		read_distance = list(csv.reader(open("../distance/"+dir+"_pdistance.txt"), delimiter="\t"))
		list_distanze=[]
		for line in read_distance[1:]:
			if line[0]=="min_ref":
					list_distanze.append(line[1])
			elif line[0]=="max_ref":
					list_distanze.append(line[1])
			elif line[0]=="median_ref":
					list_distanze.append(line[1])
		vs_reference.write(list_distanze[0]+"\t"+list_distanze[1]+"\t"+list_distanze[2]+"\t")
#write snp with consensus sequences
		read_distance = list(csv.reader(open("../distance/"+dir+"_pdistance.txt"), delimiter="\t"))
		list_distanze=[]
		for line in read_distance[1:]:
			if line[0]=="min_cons":
					list_distanze.append(line[1])
			elif line[0]=="max_cons":
					list_distanze.append(line[1])
			elif line[0]=="median_cons":
					list_distanze.append(line[1])
		vs_reference.write(list_distanze[0]+"\t"+list_distanze[1]+"\t"+list_distanze[2]+"\t")
#write consensus size
		diff_lunghezze=extract_min_max_median(consensus_size)[1]-extract_min_max_median(consensus_size)[0]
		vs_reference.write(str(extract_min_max_median(consensus_size)[0])+"\t"+str(extract_min_max_median(consensus_size)[1])+"\t"+str(diff_lunghezze)+"\t")
#write NNNsize
		item_nominmax_list = [e for e in NNN_size if e not in (extract_min_max_median(NNN_size)[0], extract_min_max_median(NNN_size)[1])]
		if item_nominmax_list:
			vs_reference.write(str(extract_min_max_median(NNN_size)[0]) + "\t" + str(extract_min_max_median(NNN_size)[1]) + "\t"+str(extract_min_max_median(item_nominmax_list)[2])+"\t")
		else:
			vs_reference.write(str(extract_min_max_median(NNN_size)[0]) + "\t" + str(extract_min_max_median(NNN_size)[1]) + "\t0\t")

#write insertions
		sorted_insertions_samples = dict(sorted(insertions_samples.items(), key=operator.itemgetter(1), reverse=True))
		conto=0
		insertions_size=[]
		insertions_numbers=[]
		for key in sorted_insertions_samples:
			insertions_size.append(sorted_insertions_samples[key][0])
			insertions_numbers.append(sorted_insertions_samples[key][1])
			if sorted_insertions_samples[key][0]!=0:
				conto+=1
		median_insertions=""
		if extract_min_max_median(insertions_size)[2]==0:
			median_insertions+="0"
		else:
			median_insertions+=str(extract_min_max_median(insertions_size)[2])+"("+str(extract_min_max_median(insertions_numbers)[2])+")"

		if sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]==0 and sorted_insertions_samples[list(sorted_insertions_samples)[0]][0]!=0:
			vs_reference.write(str(conto) + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]) + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[0]][0])+"("+str(sorted_insertions_samples[list(sorted_insertions_samples)[0]][1])+")" + "\t" + median_insertions + "\t")
		elif sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]!=0 and sorted_insertions_samples[list(sorted_insertions_samples)[0]][0]!=0:
			vs_reference.write(str(conto) + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]) +"("+str(sorted_insertions_samples[list(sorted_insertions_samples)[-1]][1])+")" + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[0]][0])+"("+str(sorted_insertions_samples[list(sorted_insertions_samples)[0]][1])+")" + "\t" + median_insertions + "\t")
		elif sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]==0 and sorted_insertions_samples[list(sorted_insertions_samples)[0]][0]==0:
			vs_reference.write(str(conto) + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]) + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[0]][0])+"\t"+ str(extract_min_max_median(insertions_size)[2]) + "\t")
#write deletions
		sorted_deletions_samples = dict(sorted(deletions_samples.items(), key=operator.itemgetter(1), reverse=True))
		conto=0
		deletions_size=[]
		deletions_numbers=[]
		for key in sorted_deletions_samples:
			deletions_size.append(sorted_deletions_samples[key][0])
			deletions_numbers.append(sorted_deletions_samples[key][1])
			if sorted_deletions_samples[key][0]!=0:
				conto+=1
		median_deletions=""
		if extract_min_max_median(deletions_size)[2]==0:
			median_deletions+="0"
		else:
			median_deletions+=str(extract_min_max_median(deletions_size)[2])+"("+str(extract_min_max_median(deletions_numbers)[2])+")"
		if sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]==0 and sorted_deletions_samples[list(sorted_deletions_samples)[0]][0]!=0:
			vs_reference.write(str(conto) + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]) + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[0]][0])+"("+str(sorted_deletions_samples[list(sorted_deletions_samples)[0]][1])+")" + "\t" + median_deletions +"\t")
		elif sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]!=0 and sorted_deletions_samples[list(sorted_deletions_samples)[0]][0]!=0:
			vs_reference.write(str(conto) + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]) +"("+str(sorted_deletions_samples[list(sorted_deletions_samples)[-1]][1])+")" + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[0]][0])+"("+str(sorted_deletions_samples[list(sorted_deletions_samples)[0]][1])+")" + "\t" + median_deletions +"\t")
		elif sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]==0 and sorted_deletions_samples[list(sorted_deletions_samples)[0]][0]==0:
			vs_reference.write(str(conto) + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]) + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[0]][0])+"\t"+ str(extract_min_max_median(deletions_size)[2]) + "\t")
#write degenerations
		vs_reference.write(str(extract_min_max_median(degenerations_size)[0]) + "\t" + str(
			extract_min_max_median(degenerations_size)[1]) + "\t")
#write lineages
		lineages_set = set(lineages)
		count_lineages = {}
		for l in lineages_set:
			count_lineages[l] = lineages.count(l)
		sorted(count_lineages.items(), key=lambda item: item[1])
		for key in count_lineages:
			vs_reference.write(key + "(" + str(count_lineages[key]) + "), ")

		os.chdir("../")
	elif directory=="consensus/":
		os.chdir(directory)
		vs_consensus=open("covrin_consensus_"+dir+".txt", "w")
		vs_consensus.write("Sample\tlenght\tdegenerations\tSNP\tNNN_size\tinsertions\tdeletions\n")
		alignment=sorted(glob.glob("*_compare.fasta"))

		consensus_size=[]
		NNN_size=[]
		insertions_size=[]
		deletions_size=[]
		degenerations_size=[]
		lineages=[]
		lenght_to_compare=[]
		insertions_samples={}
		deletions_samples={}
		for file_sequence in alignment:
			riga=file_sequence.replace("_cons_compare.fasta","")+"\t"
			r = list(SeqIO.parse(file_sequence, "fasta"))
			reference=r[0].seq
			to_analyze=r[1].seq
			for updated in update:
				if updated.id==r[1].id:
					lenght_seq=str(count_nucleotide(updated.seq))
					lenght_to_compare.append(int(lenght_seq))
					break
			consensus_size.append(int(lenght_seq))
			n_differences=find_SNP(reference, to_analyze)
			n_degenerations=str(count_deg(to_analyze))
			degenerations_size.append(int(n_degenerations))
			read_distance = list(csv.reader(open("../distance/"+dir+"_pdistance.txt"), delimiter="\t"))
			number_SNP=""
			for line in read_distance[2:]:
				if line[0]==file_sequence.replace("_ref_ALL.fasta", ""):
					number_SNP=line[1]
	#print("-----")
	#print(n_differences)
	#print("-----")
			number_NNN=str(count_N(to_analyze))
			NNN_size.append(int(number_NNN))
			n_insertions = find_insertions(reference, to_analyze)
			number_ins = len(n_insertions)
			total_ins=0
			for key in n_insertions:
				total_ins+=n_insertions[key]
			insertions_samples[file_sequence]=[total_ins, number_ins]

	#print(n_insertions)
			n_deletions = find_deletions(to_analyze)
			total_del = count_gap(to_analyze) - find_del_at_beg_end(to_analyze)[0] - find_del_at_beg_end(to_analyze)[1]
			number_del = len(n_deletions)
			deletions_samples[file_sequence]=[total_del, number_del]
	#print("DEELEZIONI")
	#print(n_deletions)
	#print("del--------")
			riga += lenght_seq + "\t" + n_degenerations + "\t" + number_SNP + "\t" + number_NNN + "\t" + str(total_ins) + "\t" + str(total_del) + "\n"
			vs_consensus.write(riga)
			lineages.append(find_lineage_in_Pangolin(file_sequence.replace("_cons_compare.fasta",""),"../../pangolin_covrin.csv"))

		risultati="Final_results"+dir+"\tSNP_min\tSNP_max\tSNP_median\tpd_min\tpd_max\tpd_median\tCons_size_min\tCons_size_max\tCons_size_diff\tNNN_size_min\tNNN_size_max" \
				  "\tNNNsize_noMinMax\tNumb_ins\tins_min\tins_max\tins_median\tNumb_del\tdel_min\tdel_max\tdel_median\tDeg_min\tDeg_max\tlineage\n"
		vs_consensus.write(risultati)
		vs_consensus.write("\t")
#write snp with reference
		read_distance = list(csv.reader(open("../distance/"+dir+"_pdistance.txt"), delimiter="\t"))
		list_distanze=[]
		for line in read_distance[1:]:
			if line[0]=="min_ref":
					list_distanze.append(line[1])
			elif line[0]=="max_ref":
					list_distanze.append(line[1])
			elif line[0]=="median_ref":
					list_distanze.append(line[1])
		vs_consensus.write(list_distanze[0]+"\t"+list_distanze[1]+"\t"+list_distanze[2]+"\t")
#write snp with consensus  sequences
		read_distance = list(csv.reader(open("../distance/"+dir+"_pdistance.txt"), delimiter="\t"))
		list_distanze=[]
		for line in read_distance[1:]:
			if line[0]=="min_cons":
					list_distanze.append(line[1])
			elif line[0]=="max_cons":
					list_distanze.append(line[1])
			elif line[0]=="median_cons":
					list_distanze.append(line[1])
		vs_consensus.write(list_distanze[0]+"\t"+list_distanze[1]+"\t"+list_distanze[2]+"\t")
#write consensus size
		diff_lunghezze=extract_min_max_median(consensus_size)[1]-extract_min_max_median(consensus_size)[0]
		vs_consensus.write(str(extract_min_max_median(consensus_size)[0])+"\t"+str(extract_min_max_median(consensus_size)[1])+"\t"+str(diff_lunghezze)+"\t")
#write NNNsize
		item_nominmax_list = [e for e in NNN_size if e not in (extract_min_max_median(NNN_size)[0], extract_min_max_median(NNN_size)[1])]
		if item_nominmax_list:
			vs_consensus.write(str(extract_min_max_median(NNN_size)[0]) + "\t" + str(extract_min_max_median(NNN_size)[1]) + "\t"+str(extract_min_max_median(item_nominmax_list)[2])+"\t")
		else:
			vs_consensus.write(str(extract_min_max_median(NNN_size)[0]) + "\t" + str(extract_min_max_median(NNN_size)[1]) + "\t0\t")
#write insertions
		sorted_insertions_samples = dict(sorted(insertions_samples.items(), key=operator.itemgetter(1), reverse=True))
		conto=0
		insertions_size=[]
		insertions_numbers=[]
		for key in sorted_insertions_samples:
			insertions_size.append(sorted_insertions_samples[key][0])
			insertions_numbers.append(sorted_insertions_samples[key][1])
			if sorted_insertions_samples[key][0]!=0:
				conto+=1
		median_insertions=""
		if extract_min_max_median(insertions_size)[2]==0:
			median_insertions+="0"
		else:
			median_insertions+=str(extract_min_max_median(insertions_size)[2])+"("+str(extract_min_max_median(insertions_numbers)[2])+")"
		if sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]==0 and sorted_insertions_samples[list(sorted_insertions_samples)[0]][0]!=0:
			vs_consensus.write(str(conto) + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]) + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[0]][0])+"("+str(sorted_insertions_samples[list(sorted_insertions_samples)[0]][1])+")" + "\t" + median_insertions+"\t")
		elif sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]!=0 and sorted_insertions_samples[list(sorted_insertions_samples)[0]][0]!=0:
			vs_consensus.write(str(conto) + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]) +"("+str(sorted_insertions_samples[list(sorted_insertions_samples)[-1]][1])+")" + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[0]][0])+"("+str(sorted_insertions_samples[list(sorted_insertions_samples)[0]][1])+")" + "\t" + median_insertions +"\t")
		elif sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]==0 and sorted_insertions_samples[list(sorted_insertions_samples)[0]][0]==0:
			vs_consensus.write(str(conto) + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[-1]][0]) + "\t" + str(sorted_insertions_samples[list(sorted_insertions_samples)[0]][0])+"\t"+ str(extract_min_max_median(insertions_size)[2]) + "\t")
#write deletions
		sorted_deletions_samples = dict(sorted(deletions_samples.items(), key=operator.itemgetter(1), reverse=True))
		conto=0
		deletions_size=[]
		deletions_numbers=[]
		for key in sorted_deletions_samples:
			deletions_size.append(sorted_deletions_samples[key][0])
			deletions_numbers.append(sorted_deletions_samples[key][1])
			if sorted_deletions_samples[key][0]!=0:
				conto+=1
		median_deletions=""
		if extract_min_max_median(deletions_size)[2]==0:
			median_deletions+="0"
		else:
			median_deletions+=str(extract_min_max_median(deletions_size)[2])+"("+str(extract_min_max_median(deletions_numbers)[2])+")"
		if sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]==0 and sorted_deletions_samples[list(sorted_deletions_samples)[0]][0]!=0:
			vs_consensus.write(str(conto) + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]) + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[0]][0])+"("+str(sorted_deletions_samples[list(sorted_deletions_samples)[0]][1])+")" + "\t" + median_deletions +"\t")
		elif sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]!=0 and sorted_deletions_samples[list(sorted_deletions_samples)[0]][0]!=0:
			vs_consensus.write(str(conto) + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]) +"("+str(sorted_deletions_samples[list(sorted_deletions_samples)[-1]][1])+")" + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[0]][0])+"("+str(sorted_deletions_samples[list(sorted_deletions_samples)[0]][1])+")" + "\t" + median_deletions +"\t")
		elif sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]==0 and sorted_deletions_samples[list(sorted_deletions_samples)[0]][0]==0:
			vs_consensus.write(str(conto) + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[-1]][0]) + "\t" + str(sorted_deletions_samples[list(sorted_deletions_samples)[0]][0])+"\t"+ str(extract_min_max_median(deletions_size)[2]) + "\t")

#write degenerations
		vs_consensus.write(str(extract_min_max_median(degenerations_size)[0])+"\t"+str(extract_min_max_median(degenerations_size)[1])+"\t")
#write lineages
		lineages_set=set(lineages)
		count_lineages={}
		for l in lineages_set:
			count_lineages[l]=lineages.count(l)
		sorted(count_lineages.items(), key=lambda item: item[1])
		for key in count_lineages:
			vs_consensus.write(key+"("+str(count_lineages[key])+"), ")
		os.chdir("../")
	return True

extract_seq()

subprocess.call("mkdir reference", shell=True)
subprocess.call("mkdir consensus", shell=True)
subprocess.call("mkdir distance", shell=True)
subprocess.call("cp "+dir+".fasta distance/", shell=True)

os.chdir("distance/")
subprocess.call("cat ../../sars2/"+reference_sequence+".fasta "+dir+".fasta > intermidiate.fasta", shell=True)
subprocess.call("mafft --localpair --maxiterate 1000 intermidiate.fasta > "+dir+"_ALL.fasta", shell=True)
os.chdir("../")

sequences=[]
for i in os.listdir("./"):
	if os.path.isfile(os.path.join("./",i)) and 'COVRIN' in i:
		subprocess.call("cp "+i+" reference/", shell=True)
		subprocess.call("cp "+i+" consensus/", shell=True)
		sequences.append(i)
subprocess.call("rm COVRIN*", shell=True)

os.chdir("reference/")
for file in sequences:
	subprocess.call("cat ../../sars2/"+reference_sequence+".fasta "+file+" > "+file.replace(".fasta","")+"_ref.fasta", shell=True)
	subprocess.call("mafft --localpair --maxiterate 1000 "+file.replace(".fasta","")+"_ref.fasta > "+file.replace(".fasta","")+"_ref_ALL.fasta", shell=True)
subprocess.call("rm *_ref.fasta", shell=True)
os.chdir("../")

subprocess.call("cp "+dir+".fasta consensus/"+dir+"_for_consensus.fasta", shell=True)
os.chdir("consensus/")
subprocess.call("./../../minimap2/minimap2 -a ../../sars2/"+reference_sequence+".fasta "+dir+"_for_consensus.fasta > "+dir+"_filtered.sam", shell=True)
subprocess.call("samtools view -bt ../../sars2/"+reference_sequence+".fasta -o "+dir+"_filtered.bam "+dir+"_filtered.sam", shell=True)
subprocess.call("samtools sort "+dir+"_filtered.bam -o "+dir+"_sort_filtered.bam", shell=True)
subprocess.call("samtools index "+dir+"_sort_filtered.bam", shell=True)
subprocess.call("samtools mpileup -aa -A -d 0 -Q 0 "+dir+"_sort_filtered.bam | ivar consensus -q 1 -t 0.2 -m 1 -n 1 -p "+dir+"_consensus", shell=True)
subprocess.call("rm *filtered.sam", shell=True)
subprocess.call("rm *filtered.bam", shell=True)
subprocess.call("cat "+dir+"_consensus.fa "+dir+"_for_consensus.fasta > correct"+dir+"_consensus.fa", shell=True)
subprocess.call("mafft --localpair --maxiterate 1000 correct"+dir+"_consensus.fa > correct"+dir+"_ALL.fasta", shell=True)
change_reference("correct"+dir+"_ALL.fasta")

for file in sequences:
	subprocess.call("cat "+dir+"_New_consensus_only.fasta "+file+" > "+file.replace(".fasta","")+"_ref.fasta", shell=True)
	subprocess.call("mafft --localpair --maxiterate 1000 "+file.replace(".fasta","")+"_ref.fasta > "+file.replace(".fasta","")+"_cons_compare.fasta", shell=True)
subprocess.call("rm *_ref.fasta", shell=True)
os.chdir("../")
#collect results for table graph
tab=open("../table_graph.txt","a")
tab.write("Sample\tconsensus_size\tn_insertions\tsum_size_insertions (bp)\tn_deletions\tsum_size_deletions (bp)\tSNPs to REF\tNNs_size\tDegenerations\n")

collect_results("reference/")
collect_results("consensus/")