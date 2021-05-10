import sys
import os
import math
import numpy as np
import random
from os import path
import argparse
from Bio.Data import CodonTable
table = CodonTable.ambiguous_dna_by_id[1]
from Bio.Seq import _translate_str
import time
from ete3 import Tree

#©EMBL-European Bioinformatics Institues, 2021

#run phylogeographic parsimony SARS-CoV-2 analyses
createDiffFile=False
reduceDiffFile=True
fromDiffToFasta=True
runFastTree=True
prepareTrees=True
parsimonyUK=True
combineAllVOCs=True

VOCs=["A.23.1","B.1.1.318","B.1.1.7","B.1.617","B.1.351","B.1.525","B.1.526","P.1","P.2","R.1"]
VOCs=["A.23.1","B.1.1.318","B.1.525","B.1.617","B.1.351"]
VOCs=["B.1.617.1","B.1.617.2","B.1.617.3"]
countsVOCs={}
countsVOCs2={}
moritzSamples={}
gisaidSamples={}
for v in VOCs:
	countsVOCs[v]=[0,0]
	countsVOCs2[v]=[0,0]
	moritzSamples[v]={}
	gisaidSamples[v]={}

lineageDict={}
#fileMoritz=open("/Users/demaio/Desktop/Moritz_analysis/cog_global_2021-04-14_mutations.is_rapid_seq.is_surveillance.csv")
fileMoritz=open("/Users/demaio/Desktop/Moritz_analysis/cog_global_2021-05-04_mutations.is_rapid_seq.is_surveillance.csv")
headers=fileMoritz.readline()
line=fileMoritz.readline()
while line!="" and line!="\n":
	linelist=line.split(",")
	if (linelist[-1].replace("\n",""))=="True" and linelist[-2]=="True" and (linelist[2] in VOCs):
		name=linelist[0]
		date=linelist[1]
		lineage=linelist[2]
		variants=linelist[4]
		e484=linelist[10]
		lineageDict[linelist[-3]]=lineage
		if e484=="E":
			countsVOCs[lineage][0]+=1
		else:
			countsVOCs[lineage][1]+=1
		moritzSamples[lineage][name]=[e484,date,variants]
		#moritzSamples[linelist[-3]]=[lineage,e484,name,date,variants]
	line=fileMoritz.readline()
print(countsVOCs)

lineageDict2={}
#fileMeta=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/metadata_tsv_2021_04_15/metadata.tsv")
fileMeta=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/metadata_tsv_2021_05_04/metadata.tsv")
headers=fileMeta.readline()
line=fileMeta.readline()
while line!="" and line!="\n":
	linelist=line.split("\t")
	name=linelist[0].replace("hCoV-19/","")
	ID=linelist[2]
	date=linelist[3]
	location=linelist[4]
	lineage=linelist[11]
	variants=linelist[14]
	if lineage in VOCs:
		lineageDict2[ID]=lineage
		if "E484K" in variants:
			e484="K"
			countsVOCs2[lineage][1]+=1
		elif "E484Q" in variants:
			e484="Q"
			countsVOCs2[lineage][1]+=1
		else:
			e484="E"
			countsVOCs2[lineage][0]+=1
		gisaidSamples[lineage][name]=[e484,date,variants,ID,location]
	line=fileMeta.readline()
print(countsVOCs2)






alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,1.0,0.0,1.0],"r":[1.0,0.0,1.0,0.0],"w":[1.0,0.0,0.0,1.0],"s":[0.0,1.0,1.0,0.0],"k":[0.0,0.0,1.0,1.0],"m":[1.0,1.0,0.0,0.0],"d":[1.0,0.0,1.0,1.0],"v":[1.0,1.0,1.0,0.0],"h":[1.0,1.0,0.0,1.0],"b":[0.0,1.0,1.0,1.0]}

#collect reference
def collectReference(fileName):
	file=open(fileName)
	line=file.readline()
	ref=""
	while line!="":
		line=file.readline()
		ref+=line.replace("\n","")
	lRef=len(ref)
	print("Ref genome length: "+str(lRef))
	file.close()
	#vector to count how many bases of each type are cumulatively in the reference genome up to a certain position
	cumulativeBases=np.zeros((lRef+1,4))
	for i in range(lRef):
		for k in range(4):
			cumulativeBases[i+1][k]=cumulativeBases[i][k]
		cumulativeBases[i+1][allelesLow[ref[i]]]+=1
	#print(cumulativeBases)
	print(cumulativeBases[-1])
	rootFreqs=np.zeros(4)
	rootFreqsLog=np.zeros(4)
	for i in range(4):
		rootFreqs[i]=cumulativeBases[-1][i]/float(lRef)
		rootFreqsLog[i]=math.log(rootFreqs[i])
	#print(rootFreqs)
	#print(rootFreqsLog)
	return ref, cumulativeBases, rootFreqs, rootFreqsLog


ref, cumulativeBases, rootFreqs, rootFreqsLog = collectReference("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/EPI_ISL_402124_lowercase.fasta")
lRef=len(ref)

# if createDiffFile:
# 	fileI=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/msa_2021-05-06/2021-05-06_unmasked.fa")
# 	#fileI=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked.fa")
# 	line=fileI.readline()
# 	nSeqs=0
# 	while line!="" and line!="\n":
# 		nSeqs+=1
# 		seq=""
# 		name=line.replace(">","").replace("\n","")
# 		line=fileI.readline(10000000)
# 		if len(line)>10*lRef:
# 			print("name line")
# 			print(line)
# 			print(len(line))
# 			exit()
# 		while line!="" and line!="\n" and line[0]!=">":
# 			seq+=line.replace("\n","")
# 			line=fileI.readline(10000000)
# 			if len(seq)>100*lRef:
# 				print("seq line")
# 				for i in range(int(len(seq)/1000)):
# 					print(i)
# 					print(seq[i*1000:(i+1)*1000])
# 				print(seq)
# 				print(len(seq))
# 				exit()
# 		if (nSeqs%1000)==0:
# 			print(nSeqs)
# exit()

if createDiffFile:
	start = time.time()
	#collect alignment and translate into diff file
	#fileI=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked.fa")
	#fileI=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/msa_2021-05-06/2021-05-06_unmasked.fa")
	fileI=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked.fa")
	#fileI=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/msa_0501/msa_0501.fasta")
	#fileO=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_differences.txt","w")
	fileO=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_differences_new.txt","w")
	line=fileI.readline()
	nSeqs=0
	while line!="" and line!="\n":
		nSeqs+=1
		seq=""
		name=line.replace(">","").replace("\n","")
		fileO.write(line)
		line=fileI.readline()
		while line!="" and line!="\n" and line[0]!=">":
			seq+=line.replace("\n","")
			line=fileI.readline()
			if len(seq)>10*lRef:
				print(seq)
				exit()
		lSeq=len(seq)
		if lSeq!=lRef:
			print("Seq "+name+" has length "+str(len(seq))+" while ref is "+str(lRef))
			print(seq[:200000])
			#exit()
			if lSeq>lRef:
				lSeq=lRef
		
		# state 0=ref; 1=N; 2=-; 
		state=0
		seqList=[]
		length=0
		for i in range(lSeq):
			if state==1:
				if seq[i]=="n":
					length+=1
				else:
					seqList.append(("n",i+1-length,length))
					length=0
					state=0
			elif state==2:
				if seq[i]=="-":
					length+=1
				else:
					seqList.append(("-",i+1-length,length))
					length=0
					state=0
			elif seq[i]=="n" and state!=1:
				length=1
				state=1
			elif seq[i]=="-" and state!=2:
				length=1
				state=2
			if seq[i]!=ref[i] and seq[i]!="-" and seq[i]!="n":
				seqList.append((seq[i],i+1))
		if lSeq<lRef:
			seqList.append(("n",lSeq+1-length,length+lRef-lSeq))
		else:
			if state==1:
					seqList.append(("n",lRef+1-length,length))
			elif state==2:
					seqList.append(("-",lRef+1-length,length))
		
		for m in seqList:
			if len(m)==2:
				fileO.write(m[0]+"\t"+str(m[1])+"\n")
			else:
				fileO.write(m[0]+"\t"+str(m[1])+"\t"+str(m[2])+"\n")
		if (nSeqs%1000)==0:
			print(nSeqs)
	fileI.close()
	fileO.close()

	time2 = time.time() - start
	print("Time to convert alignment file: "+str(time2))
	print(str(nSeqs)+" sequences converted.")
	#Time take to convert alignment file: 15593.673340320587
	#915508 sequences converted.



def readConciseAlignment(fileName):
	start = time.time()
	fileI=open(fileName)
	line=fileI.readline()
	nSeqs=0
	data={}
	while line!="" and line!="\n":
		nSeqs+=1
		seqList=[]
		name=line.replace(">","").replace("\n","")
		line=fileI.readline()
		while line!="" and line!="\n" and line[0]!=">":
			linelist=line.split()
			if len(linelist)>2:
				entry=(linelist[0],int(linelist[1]),int(linelist[2]))
			else:
				entry=(linelist[0],int(linelist[1]))
			seqList.append(entry)
			line=fileI.readline()
		data[name]=seqList
	fileI.close()
	time2 = time.time() - start
	print("Time to read DNA reduced data file: "+str(time2))
	print(str(nSeqs)+" sequences in file.")
	return data


if reduceDiffFile:
	#read sequence data from file
	#data=readConciseAlignment("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_differences.txt")
	#data=readConciseAlignment("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/msa_2021-05-06/2021-05-06_unmasked_differences.txt")
	data=readConciseAlignment("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_differences_new.txt")
	dataReduced={}
	for v in VOCs:
		dataReduced[v]={}
		print(v)
		nameDict={}
		names1=moritzSamples[v].keys()
		names2=gisaidSamples[v].keys()
		IDs=[]
		for n in names2:
			IDs.append(gisaidSamples[v][n][3])
			nameDict[gisaidSamples[v][n][3]]=n
		IDs2=[]
		for n in names1:
			if n in names2:
				IDs2.append(gisaidSamples[v][n][3])

		count=0
		count2=0
		for sample in data:
			if sample in IDs:
				count+=1
				dataReduced[v][sample]=data[sample]
				if sample in IDs2:
					count2+=1
		
		print("Samples from Moritz in the alignment")
		print(count2)
		print("Samples from GISAID metadata in the alignment")
		print(count)
		print("Total samples from GISAID metadata for this lineage:")
		print(len(IDs))
		print("")


		fileO=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_differences_"+v+".txt","w")
		for k in dataReduced[v].keys():
			fileO.write(">"+k+"\n")
			for m in dataReduced[v][k]:
				if len(m)==2:
					fileO.write(m[0]+"\t"+str(m[1])+"\n")
				else:
					fileO.write(m[0]+"\t"+str(m[1])+"\t"+str(m[2])+"\n")
		fileO.close()
	#data=dataReduced

	


pos484=21563+483*3
print("Position of E484K")
print(pos484)


def sequenceFromDiff(diffs,ref):
	seq=""
	pos=1
	for m in diffs:
			currPos=m[1]
			if currPos>pos:
				#region identical to the ref.
				seq+=ref[pos-1:currPos-1]
				pos=currPos
			if m[0]=="n" or m[0]=="-":
				length=m[2]
				#region with ambiguity.
				seq+=(m[0]*length)
				pos=currPos+length
			else:
				seq+=m[0]
				pos=currPos+1
	if pos<=lRef:
		seq+=ref[pos-1:]
	# print("Final length:")
	# print(len(seq))
	# count=0
	# for i in range(lRef):
	# 	if ref[i]!=seq[i] and seq[i]!="n" and seq[i]!="-":
	# 		count+=1
	# print(count)
	if len(seq)!=len(ref):
		print("Problem with sequence length!")
		print(diffs)
		print(len(seq))
		exit()
	return seq


if fromDiffToFasta:
	if not combineAllVOCs:
		for v in VOCs:
			#data=readConciseAlignment("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_differences_"+v+".txt")
			data=readConciseAlignment("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_differences_"+v+".txt")
			#fileO=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_"+v+".fa","w")
			fileO=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_"+v+".fa","w")
			for sample in data:
				fileO.write(">"+sample+"\n"+sequenceFromDiff(data[sample],ref)+"\n")
			fileO.close()
	if combineAllVOCs:
		fileO=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_allVOCs.fa","w")
		for v in VOCs:
			#data=readConciseAlignment("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_differences_"+v+".txt")
			data=readConciseAlignment("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_differences_"+v+".txt")
			#fileO=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_"+v+".fa","w")
			#fileO=open("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/msa_2021-05-06/2021-05-06_unmasked_allVOCs.fa","w")
			for sample in data:
				fileO.write(">"+sample+"\n"+sequenceFromDiff(data[sample],ref)+"\n")
		fileO.close()

if runFastTree:
	if combineAllVOCs:
		os.system("/Applications/FastTree/FastTree -quiet -nosupport -nt -gtr -nocat < /Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_allVOCs.fa > /Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_allVOCs_fastTree.tree")
	else:
		for v in VOCs:
			#os.system("/Applications/FastTree/FastTree -quiet -nosupport -nt -gtr -nocat < /Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_"+v+".fa > /Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_"+v+"_fastTree.tree")
			os.system("/Applications/FastTree/FastTree -quiet -nosupport -nt -gtr -nocat < /Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_"+v+".fa > /Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_"+v+"_fastTree.tree")

	

#function to run likelihood calculation in a subsampled tree to compare it to full likelihood calculation in other methods.
def countLeaves(phylo,newSamples):
		if len(phylo.children)==0:
			if phylo.name in newSamples:
				phylo.descendants=1
			else:
				phylo.descendants=0
		else:
			tot=0
			for c in phylo.children:
				countLeaves(c,newSamples)
				tot+=c.descendants
			phylo.descendants=tot

def extractSubtree(phylo,t1,minBlen=0.0000001):
		if phylo.descendants>0:
			nChildren=0
			t1.dist+=phylo.dist
			if len(phylo.children)==0:
				t1.name=phylo.name
				#t1.dist+=phylo.dist
			newChildren=[]
			for c in phylo.children:
				if c.descendants>0:
					nChildren+=1
					newChildren.append(c)
			if nChildren==1:
				#t1.dist+=phylo.dist
				extractSubtree(newChildren[0],t1,minBlen=minBlen)
			if nChildren>1:
				for i in range(nChildren-1):
					t1.add_child(dist=minBlen)
					extractSubtree(newChildren[i],t1.children[0],minBlen=minBlen)
					if i!=(nChildren-2):
						t1.add_child(dist=minBlen)
						t1=t1.children[1]
				t1.add_child(dist=minBlen)
				extractSubtree(newChildren[-1],t1.children[1],minBlen=minBlen)




minDist=0.000001
def parsimonyTraverse(node,IDs2,gisaidSamples,nameDict):
	printExit=False
	#0: min into surv UK conditional on ancestral UK; 1: max into surv UK conditional on ancestral UK; 2: total number of state changes conditional on UK; 3-5: conditional on non-UK
	newCounts=np.zeros(6, dtype=int)
	node.scores=[0,0]
	node.childrenPars=[]
	node.cladeChildren=[]
	if len(node.children)==0:
		if node.name in IDs2:
			#UK surveillance
			state=0
			newCounts[3]=1
			newCounts[4]=1
			newCounts[5]=100
			node.scores[1]=100
			cladeSize=[1,0]
			pastClades=[[],[1]]
		else:
			country=gisaidSamples[nameDict[node.name]][-1].split("/")[1].replace(" ","")
			if country=="UnitedKingdom":
				#UK non-surveillance
				state=1
				newCounts[5]=100
				node.scores[1]=100
				cladeSize=[0,0]
				pastClades=[[],[]]
			else:
				#non-uK
				state=2
				newCounts[2]=100
				node.scores[0]=100
				cladeSize=[0,0]
				pastClades=[[],[]]
		return [newCounts],[cladeSize],[pastClades],[node]
	else:
		nodeList=[]
		currentClade=[]
		cladeSizes=[]
		for c in node.children:
			if c.name=='Agamennon':
			#if c.name=='EPI_ISL_1440292':
			#if c.name=='EPI_ISL_1099865':
				printExit=True
			countsChild,cladeSize,cladesChild,nodes=parsimonyTraverse(c,IDs2,gisaidSamples,nameDict)
			for ch in range(len(countsChild)):
				nodeList.append(countsChild[ch])
				currentClade.append(cladeSize[ch])
				cladeSizes.append(cladesChild[ch])
				#node.scores=[0,0]
				node.childrenPars.append([0,0])
				node.cladeChildren.append(nodes[ch])
		if node.dist<minDist and (not node.is_root()):
			#print("internal node 0 top branch length")
			#print(nodeList)
			return nodeList, currentClade, cladeSizes, node.cladeChildren
		else:
			if len(nodeList)==1:
				print("Strange, one-childed parent?")
				return nodeList, currentClade, cladeSizes, node.cladeChildren
			found01=0 #does it make sense to have at least a 0->1 transition among one of the two final child branches of the current node?
			found10=0
			#cladeSizeNew=[[0,[]],[0,[]]]
			currentCladeNew=[0,0]
			cladeListNew=[[],[]]
			for n in nodeList:
				if n[5]<n[2]:
					found01+=1
					if n[5]+1<n[2]:
						found01+=1
				if n[2]<n[5]:
					found10+=1
					if n[2]+1<n[5]:
						found10+=1
			if found01==0: # no 0->1 transition happened among children if we have 0 ancestor
				chIndex=0
				for n in nodeList:
					newCounts[0]+=n[0]
					newCounts[1]+=n[1]
					newCounts[2]+=n[2]
					currentCladeNew[0]+=currentClade[chIndex][0]
					for cl in cladeSizes[chIndex][0]:
						cladeListNew[0].append(cl)
					node.scores[0]+=node.cladeChildren[chIndex].scores[0]
					#node.childrenPars.append([0,0])
					#node.cladeChildren.append(nodes[ch])	
					chIndex+=1
			else:
				ifNoMig=[0,0,0]
				for n in nodeList:
					ifNoMig[0]+=n[0]
					ifNoMig[1]+=n[1]
					ifNoMig[2]+=n[2]
				ifMig=[0,0,1]
				for n in nodeList:
					if n[2]<n[5]:
						ifMig[0]+=n[0]
						ifMig[1]+=n[1]
						ifMig[2]+=n[2]
					elif n[5]<n[2]:
						ifMig[0]+=n[3]
						ifMig[1]+=n[4]
						ifMig[2]+=n[5]
					else:
						ifMig[2]+=n[5]
						if n[3]>n[0]:
							ifMig[0]+=n[0]
						else:
							ifMig[0]+=n[3]
						if n[4]>n[1]:
							ifMig[1]+=n[4]
						else:
							ifMig[1]+=n[1]
				if ifNoMig[2]<ifMig[2]:
					newCounts[0]=ifNoMig[0]
					newCounts[1]=ifNoMig[1]
					newCounts[2]=ifNoMig[2]
				elif ifNoMig[2]>ifMig[2]:
					newCounts[0]=ifMig[0]
					newCounts[1]=ifMig[1]
					newCounts[2]=ifMig[2]
				else:
					newCounts[2]=ifMig[2]
					if ifNoMig[0]<ifMig[0]:
						newCounts[0]=ifNoMig[0]
					else:
						newCounts[0]=ifMig[0]
					if ifNoMig[1]<ifMig[1]:
						newCounts[1]=ifMig[1]
					else:
						newCounts[1]=ifNoMig[1]

				# chIndex=0
				# for n in nodeList:
				# 	if "1589857" in node.cladeChildren[chIndex].name:
				# 		print(ifNoMig)
				# 		print(ifMig)
				# 		print(n)
				# 		exit()
				# 	chIndex+=1
				chIndex=0
				if ifNoMig[2]<=ifMig[2]:
					for n in nodeList:
						currentCladeNew[0]+=currentClade[chIndex][0]
						for cl in cladeSizes[chIndex][0]:
							cladeListNew[0].append(cl)
						node.scores[0]+=node.cladeChildren[chIndex].scores[0]
						#node.childrenPars.append([0,0])
						#node.cladeChildren.append(nodes[ch])	
						chIndex+=1
				else:
					#migratedCladeSize=0
					node.scores[0]+=1
					for n in nodeList:
						if n[2]<n[5]:
							currentCladeNew[0]+=currentClade[chIndex][0]
							for cl in cladeSizes[chIndex][0]:
								cladeListNew[0].append(cl)
							node.scores[0]+=node.cladeChildren[chIndex].scores[0]
							#node.childrenPars.append([0,0])
							#node.cladeChildren.append(nodes[ch])	
						else:
							#migratedCladeSize+=currentClade[chIndex][0]
							for cl in cladeSizes[chIndex][1]:
								cladeListNew[0].append(cl)
							node.scores[0]+=node.cladeChildren[chIndex].scores[1]
							node.childrenPars[chIndex][0]=1
							#node.cladeChildren.append(nodes[ch])	
						chIndex+=1

			#REPEAT THE ABOVE for ancestral state 1 - in this case we probably don't need to add the extra migration to min and max counts.
			if found10==0: # no 1->0 transition happened among children if we have 1 ancestor
				chIndex=0
				for n in nodeList:
					newCounts[3]+=n[3]
					newCounts[4]+=n[4]
					newCounts[5]+=n[5]
					for cl in cladeSizes[chIndex][1]:
						cladeListNew[1].append(cl)
					node.scores[1]+=node.cladeChildren[chIndex].scores[1]
					node.childrenPars[chIndex][1]=1
					#node.cladeChildren.append(nodes[ch])	
					chIndex+=1

			else:
				ifNoMig=[0,0,0]
				for n in nodeList:
					ifNoMig[0]+=n[3]
					ifNoMig[1]+=n[4]
					ifNoMig[2]+=n[5]
				ifMig=[0,0,1]
				foundDesc=False
				foundDescMin=False
				foundDescMax=False
				for n in nodeList:
					if n[2]<n[5]:
						ifMig[0]+=n[0]
						ifMig[1]+=n[1]
						ifMig[2]+=n[2]
						if not foundDesc:
							if n[4]>0:
								foundDesc=True
								ifMig[0]+=1
								ifMig[1]+=1
					elif n[5]<n[2]:
						ifMig[0]+=n[3]
						ifMig[1]+=n[4]
						ifMig[2]+=n[5]
						
					else:
						ifMig[2]+=n[5]
						if n[3]>n[0]:
							ifMig[0]+=n[0]
							if (not foundDesc) and (not foundDescMin):
								if n[4]>0:
									foundDescMin=True
									ifMig[0]+=1
						else:
							ifMig[0]+=n[3]
							
						if n[1]>n[4]:
							ifMig[1]+=n[1]
							if (not foundDesc) and (not foundDescMax):
								if n[4]>0:
									foundDescMax=True
									ifMig[1]+=1
						else:
							ifMig[1]+=n[4]

				if ifNoMig[2]<ifMig[2]:
					newCounts[3]=ifNoMig[0]
					newCounts[4]=ifNoMig[1]
					newCounts[5]=ifNoMig[2]
				elif ifNoMig[2]>ifMig[2]:
					newCounts[3]=ifMig[0]
					newCounts[4]=ifMig[1]
					newCounts[5]=ifMig[2]
				else:
					newCounts[5]=ifMig[2]
					if ifNoMig[0]<ifMig[0]:
						newCounts[3]=ifNoMig[0]
					else:
						newCounts[3]=ifMig[0]
					if ifNoMig[1]<ifMig[1]:
						newCounts[4]=ifMig[1]
					else:
						newCounts[4]=ifNoMig[1]

				chIndex=0
				if ifNoMig[2]<ifMig[2]:
					for n in nodeList:
						#currentCladeNew[0]+=currentClade[chIndex][0]
						for cl in cladeSizes[chIndex][1]:
							cladeListNew[1].append(cl)
						node.scores[1]+=node.cladeChildren[chIndex].scores[1]
						node.childrenPars[chIndex][1]=1
						#node.cladeChildren.append(nodes[ch])
						chIndex+=1
				else:
					node.scores[1]+=1
					migratedCladeSize=0
					for n in nodeList:
						if n[2]<=n[5]:
							migratedCladeSize+=currentClade[chIndex][0]
							#currentCladeNew[0]+=currentClade[chIndex][0]
							for cl in cladeSizes[chIndex][0]:
								cladeListNew[1].append(cl)
							node.scores[1]+=node.cladeChildren[chIndex].scores[0]
							node.childrenPars[chIndex][1]=0
							#node.cladeChildren.append(nodes[ch])	
						else:
							#migratedCladeSize+=currentClade[chIndex][0]
							for cl in cladeSizes[chIndex][1]:
								cladeListNew[1].append(cl)
							node.scores[1]+=node.cladeChildren[chIndex].scores[1]
							node.childrenPars[chIndex][1]=1
							#node.cladeChildren.append(nodes[ch])
						chIndex+=1
					if migratedCladeSize>0:
						cladeListNew[1].append(migratedCladeSize)

			#if newCounts[2]>1:
			if printExit:
				print("internal node non-0 top branch length")
				print(nodeList)
				print(newCounts)
				exit()
			return [newCounts], [currentCladeNew], [cladeListNew], [node]


leafDict={}
cladeCount=[1]
def assignClade(phylo,state,clade):
	newClade=clade
	printing=False
	if "1589857" in phylo.name:
		printing=True
	for n in range(len(phylo.cladeChildren)):
		if "1589857" in phylo.cladeChildren[n].name:
			printing=True
	if printing:
		print("printing ")
		print(state)
		print(clade)
		print(len(phylo.cladeChildren))
	foundMig=False
	if state==1:
		for n in range(len(phylo.cladeChildren)):
			chState=phylo.childrenPars[n][state]
			if chState==0:
				foundMig=True
				break
		if foundMig:
			newClade=cladeCount[0]
			cladeCount[0]+=1
	if printing:
		print(foundMig)
		print(newClade)
		print("children")
	for n in range(len(phylo.cladeChildren)):
		chState=phylo.childrenPars[n][state]
		chNode=phylo.cladeChildren[n]
		if printing:
			print(phylo.childrenPars[n])
		if state==0 and chState==1:
			assignClade(chNode,chState,-1)
			chNode.cladeNum=-1
		elif state==1 and chState==0:
			assignClade(chNode,chState,newClade)
			chNode.cladeNum=newClade
			#cladeCount[0]+=1
		elif state==0 and chState==0:
			assignClade(chNode,chState,clade)
			chNode.cladeNum=clade
		else:
			assignClade(chNode,chState,-1)
			chNode.cladeNum=-1
		if printing:
			print("child "+str(n+1))
			print(chNode.cladeNum)
	phylo.cladeNum=clade
	if printing:
		print(phylo.cladeNum)
	if phylo.is_leaf:
		leafDict[phylo.name]=phylo






#infer number of UK introductions with parsimony
if parsimonyUK:
	if not combineAllVOCs:
		for v in VOCs:
			print(v)
			cladeCount[0]=1
			nameDict={}
			#phylo = Tree("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_"+v+"_fastTree.tree")
			phylo = Tree("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_"+v+"_fastTree.tree")
			names1=moritzSamples[v].keys()
			names2=gisaidSamples[v].keys()
			IDs=[]
			for n in names2:
				IDs.append(gisaidSamples[v][n][3])
				nameDict[gisaidSamples[v][n][3]]=n
			IDs2=[]
			for n in names1:
				if n in names2:
					IDs2.append(gisaidSamples[v][n][3])

			count=0
			count2=0
			for leaf in phylo:
				if leaf.name in IDs:
					count+=1
					if leaf.name in IDs2:
						count2+=1
			print(count2)
			print(count)

			parsimonyCounts, currClade, cladeList, rootChildren =parsimonyTraverse(phylo,IDs2,gisaidSamples[v],nameDict)
			if phylo.scores[0]>phylo.scores[1]:
				state=0
				phylo.cladeNum=0
				assignClade(phylo,state,0)
			else:
				state=1
				assignClade(phylo,state,-1)

			print(parsimonyCounts[0])
			print(currClade)
			print(cladeList)
			print(phylo.scores)
			print(phylo.childrenPars)

	if combineAllVOCs:
		print("Now combining all VOCs above")
		cladeCount[0]=1
		nameDict={}
		#phylo = Tree("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_"+v+"_fastTree.tree")
		phylo = Tree("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_allVOCs_fastTree.tree")
		names1=[]
		names2=[]
		IDs=[]
		IDs2=[]
		gisaidSamplesAll={}
		for v in VOCs:
			n2=gisaidSamples[v].keys()
			n1=moritzSamples[v].keys()
			names1.extend(n1)
			names2.extend(n2)
			for n in n2:
				IDs.append(gisaidSamples[v][n][3])
				nameDict[gisaidSamples[v][n][3]]=n
				gisaidSamplesAll[n]=gisaidSamples[v][n]
			for n in n1:
				if n in n2:
					IDs2.append(gisaidSamples[v][n][3])

		count=0
		count2=0
		for leaf in phylo:
			if leaf.name in IDs:
				count+=1
				if leaf.name in IDs2:
					count2+=1
		print(count2)
		print(count)

		parsimonyCounts, currClade, cladeList, rootChildren =parsimonyTraverse(phylo,IDs2,gisaidSamplesAll,nameDict)
		if phylo.scores[0]>phylo.scores[1]:
			state=0
			phylo.cladeNum=0
			assignClade(phylo,state,0)
		else:
			state=1
			assignClade(phylo,state,-1)

		print(parsimonyCounts[0])
		print(currClade)
		print(cladeList)
		print(phylo.scores)
		print(phylo.childrenPars)
		#exit()







if prepareTrees:
	if not combineAllVOCs:
		for v in VOCs:
			print(v)
			nameDict={}
			#phylo = Tree("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_"+v+"_fastTree.tree")
			phylo = Tree("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_"+v+"_fastTree.tree")
			names1=moritzSamples[v].keys()
			names2=gisaidSamples[v].keys()
			IDs=[]
			for n in names2:
				IDs.append(gisaidSamples[v][n][3])
				nameDict[gisaidSamples[v][n][3]]=n
			IDs2=[]
			for n in names1:
				if n in names2:
					IDs2.append(gisaidSamples[v][n][3])

			count=0
			count2=0
			for leaf in phylo:
				if leaf.name in IDs:
					count+=1
					if leaf.name in IDs2:
						count2+=1
			print(count2)
			print(count)
			print(len(IDs))
			print("")

			#extract subtree containing these samples and replace 0 lengths with non-zero lengths and multifurcations with bifurcations.
			#countLeaves(phylo,IDs)
			#t1 = Tree()	
			#extractSubtree(phylo,t1)
			t1=phylo
			print(str(len(t1))+" sequences in the new tree.")
			#t1.write(outfile="/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/"+v+"_fastTree_total.tree")
			t1.write(outfile="/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/"+v+"_fastTree_total_May.tree")
			file=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/"+v+"_fastTree_total_May.tree")
			t1String=file.readline()
			file.close()
			file=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/"+v+"_fastTree_total_May_features_new.tree","w")
			file.write("#NEXUS\n")
			file.write("begin taxa;\n")
			file.write("	dimensions ntax="+str(len(t1))+";\n")
			file.write("	taxlabels\n")
			UKsamples=0
			UKsamplesk=0
			nonSurvUK=0
			nonUK=0
			nonSurvUKk=0
			nonUKk=0
			for leaf in t1:
				file.write("	\'"+leaf.name+"\'\n")
				newName=""
				newName+="\'"+leaf.name+"\'"
				country=gisaidSamples[v][nameDict[leaf.name]][-1].split("/")[1].replace(" ","")
				isUK=1
				isUK2=4
				if leaf.name in IDs2:
					#newName+="[&type=\"UK_surveillance_"+gisaidSamples[v][nameDict[leaf.name]][0]+"\","
					newName+="[&type=\"UK_surveillance\","
					UKsamples+=1
					if gisaidSamples[v][nameDict[leaf.name]][0]=="K":
						UKsamplesk+=1
				elif country=="UnitedKingdom":
					newName+="[&type=\"UK\","
					nonSurvUK+=1
					if gisaidSamples[v][nameDict[leaf.name]][0]=="K":
						nonSurvUKk+=1
				else:
					newName+="[&type=\"non-UK\","
					nonUK+=1
					if gisaidSamples[v][nameDict[leaf.name]][0]=="K":
						nonUKk+=1
						#print(leaf.name)
						#print(country)
					isUK=0
					isUK2=3
				newName+="eek=\""+gisaidSamples[v][nameDict[leaf.name]][0]+"\","
				newName+="cladeNum="+str(leafDict[leaf.name].cladeNum)+","
				
				#if gisaidSamples[v][nameDict[leaf.name]][0]=="K":
				#	print(leaf.name)
				#	print(country)
				#	newName+="mutation=1,"
				#else:
				#	newName+="mutation=0,"
				newName+="isUK="+str(isUK)+","
				newName+="isUK2="+str(isUK2)+","
				newName+="location=\""+country+"\"]"
				t1String=t1String.replace(leaf.name,newName)
			print("UK samples from surveillance: "+str(UKsamples)+" of which K: "+str(UKsamplesk))
			print("UK samples not from surveillance: "+str(nonSurvUK)+" of which K: "+str(nonSurvUKk))
			print("non-UK samples: "+str(nonUK)+" of which K: "+str(nonUKk))
			file.write("	;\n")
			file.write("end;\n")
			file.write("\n")
			file.write("begin trees;\n")
			file.write("tree TREE1 = [&R] [&R=true]")
			file.write(t1String)
			file.write("\n"+"end;\n")
			file.close()

			countLeaves(phylo,IDs2)
			t1 = Tree()	
			extractSubtree(phylo,t1)
			t1.write(outfile="/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/"+v+"_fastTree_onlyUK_May.tree")
			file=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/"+v+"_fastTree_onlyUK_May.tree")
			t1String=file.readline()
			file.close()
			file=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/"+v+"_fastTree_onlyUK_May_features_new.tree","w")
			file.write("#NEXUS\n")
			file.write("begin taxa;\n")
			file.write("	dimensions ntax="+str(len(t1))+";\n")
			file.write("	taxlabels\n")
			UKsamples=0
			UKsamplesk=0
			nonSurvUK=0
			nonUK=0
			nonSurvUKk=0
			nonUKk=0
			for leaf in t1:
				file.write("	\'"+leaf.name+"\'\n")
				newName=""
				newName+="\'"+leaf.name+"\'"
				country=gisaidSamples[v][nameDict[leaf.name]][-1].split("/")[1].replace(" ","")
				isUK=1
				isUK2=4
				if leaf.name in IDs2:
					#newName+="[&type=\"UK_surveillance_"+gisaidSamples[v][nameDict[leaf.name]][0]+"\","
					newName+="[&type=\"UK_surveillance\","
					UKsamples+=1
					if gisaidSamples[v][nameDict[leaf.name]][0]=="K":
						UKsamplesk+=1
				elif country=="UnitedKingdom":
					newName+="[&type=\"UK\","
					nonSurvUK+=1
					if gisaidSamples[v][nameDict[leaf.name]][0]=="K":
						nonSurvUKk+=1
				else:
					newName+="[&type=\"non-UK\","
					nonUK+=1
					if gisaidSamples[v][nameDict[leaf.name]][0]=="K":
						nonUKk+=1
						#print(leaf.name)
						#print(country)
					isUK=0
					isUK2=3
				newName+="eek=\""+gisaidSamples[v][nameDict[leaf.name]][0]+"\","
				newName+="cladeNum="+str(leafDict[leaf.name].cladeNum)+","
				
				#if gisaidSamples[v][nameDict[leaf.name]][0]=="K":
				#	print(leaf.name)
				#	print(country)
				#	newName+="mutation=1,"
				#else:
				#	newName+="mutation=0,"
				newName+="isUK="+str(isUK)+","
				newName+="isUK2="+str(isUK2)+","
				newName+="location=\""+country+"\"]"
				t1String=t1String.replace(leaf.name,newName)
			print("UK samples from surveillance: "+str(UKsamples)+" of which K: "+str(UKsamplesk))
			print("UK samples not from surveillance: "+str(nonSurvUK)+" of which K: "+str(nonSurvUKk))
			print("non-UK samples: "+str(nonUK)+" of which K: "+str(nonUKk))
			file.write("	;\n")
			file.write("end;\n")
			file.write("\n")
			file.write("begin trees;\n")
			file.write("tree TREE1 = [&R] [&R=true]")
			file.write(t1String)
			file.write("\n"+"end;\n")
			file.close()


			#exit()

	if combineAllVOCs:
		#for v in VOCs:
		nameDict={}
		#phylo = Tree("/Users/demaio/Desktop/Moritz_analysis/GISAID/sequences_fasta_2021_04_15/msa_2021-04-17/2021-04-17_unmasked_"+v+"_fastTree.tree")
		phylo = Tree("/Users/demaio/Desktop/Moritz_analysis/GISAID_5-5-2021/sequences_fasta_2021_05_04/mmsa_2021-05-06/2021-05-06_masked_allVOCs_fastTree.tree")
		names1=[]
		names2=[]
		IDs=[]
		IDs2=[]
		gisaidSamplesAll={}
		for v in VOCs:
			n2=gisaidSamples[v].keys()
			n1=moritzSamples[v].keys()
			names1.extend(n1)
			names2.extend(n2)
			for n in n2:
				IDs.append(gisaidSamples[v][n][3])
				nameDict[gisaidSamples[v][n][3]]=n
				gisaidSamplesAll[n]=gisaidSamples[v][n]
			for n in n1:
				if n in n2:
					IDs2.append(gisaidSamples[v][n][3])

		count=0
		count2=0
		for leaf in phylo:
			if leaf.name in IDs:
				count+=1
				if leaf.name in IDs2:
					count2+=1
		print(count2)
		print(count)
		print(len(IDs))
		print("")

		#extract subtree containing these samples and replace 0 lengths with non-zero lengths and multifurcations with bifurcations.
		#countLeaves(phylo,IDs)
		#t1 = Tree()	
		#extractSubtree(phylo,t1)
		t1=phylo
		print(str(len(t1))+" sequences in the new tree.")
		#t1.write(outfile="/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/"+v+"_fastTree_total.tree")
		#t1.write(outfile="/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/"+v+"_fastTree_total_May.tree")
		t1.write(outfile="/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/allVOCs_fastTree_total_May.tree")
		file=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/allVOCs_fastTree_total_May.tree")
		t1String=file.readline()
		file.close()
		file=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/allVOCs_fastTree_total_May_features_new.tree","w")
		file.write("#NEXUS\n")
		file.write("begin taxa;\n")
		file.write("	dimensions ntax="+str(len(t1))+";\n")
		file.write("	taxlabels\n")
		UKsamples=0
		UKsamplesk=0
		nonSurvUK=0
		nonUK=0
		nonSurvUKk=0
		nonUKk=0
		for leaf in t1:
			file.write("	\'"+leaf.name+"\'\n")
			newName=""
			newName+="\'"+leaf.name+"\'"
			country=gisaidSamplesAll[nameDict[leaf.name]][-1].split("/")[1].replace(" ","")
			isUK=1
			isUK2=4
			if leaf.name in IDs2:
				#newName+="[&type=\"UK_surveillance_"+gisaidSamples[v][nameDict[leaf.name]][0]+"\","
				newName+="[&type=\"UK_surveillance\","
				UKsamples+=1
				if gisaidSamplesAll[nameDict[leaf.name]][0]=="K":
					UKsamplesk+=1
			elif country=="UnitedKingdom":
				newName+="[&type=\"UK\","
				nonSurvUK+=1
				if gisaidSamplesAll[nameDict[leaf.name]][0]=="K":
					nonSurvUKk+=1
			else:
				newName+="[&type=\"non-UK\","
				nonUK+=1
				if gisaidSamplesAll[nameDict[leaf.name]][0]=="K":
					nonUKk+=1
					#print(leaf.name)
					#print(country)
				isUK=0
				isUK2=3
			newName+="eek=\""+gisaidSamplesAll[nameDict[leaf.name]][0]+"\","
			newName+="cladeNum="+str(leafDict[leaf.name].cladeNum)+","
			
			#if gisaidSamples[v][nameDict[leaf.name]][0]=="K":
			#	print(leaf.name)
			#	print(country)
			#	newName+="mutation=1,"
			#else:
			#	newName+="mutation=0,"
			newName+="isUK="+str(isUK)+","
			newName+="isUK2="+str(isUK2)+","
			newName+="location=\""+country+"\"]"
			t1String=t1String.replace(leaf.name,newName)
		print("UK samples from surveillance: "+str(UKsamples)+" of which K: "+str(UKsamplesk))
		print("UK samples not from surveillance: "+str(nonSurvUK)+" of which K: "+str(nonSurvUKk))
		print("non-UK samples: "+str(nonUK)+" of which K: "+str(nonUKk))
		file.write("	;\n")
		file.write("end;\n")
		file.write("\n")
		file.write("begin trees;\n")
		file.write("tree TREE1 = [&R] [&R=true]")
		file.write(t1String)
		file.write("\n"+"end;\n")
		file.close()

		countLeaves(phylo,IDs2)
		t1 = Tree()	
		extractSubtree(phylo,t1)
		t1.write(outfile="/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/allVOCs_fastTree_onlyUK_May.tree")
		file=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/allVOCs_fastTree_onlyUK_May.tree")
		t1String=file.readline()
		file.close()
		file=open("/Users/demaio/Desktop/Moritz_analysis/GISAID/GISAID-hCoV-19-phylogeny-2021-04-10/allVOCs_fastTree_onlyUK_May_features_new.tree","w")
		file.write("#NEXUS\n")
		file.write("begin taxa;\n")
		file.write("	dimensions ntax="+str(len(t1))+";\n")
		file.write("	taxlabels\n")
		UKsamples=0
		UKsamplesk=0
		nonSurvUK=0
		nonUK=0
		nonSurvUKk=0
		nonUKk=0
		for leaf in t1:
			file.write("	\'"+leaf.name+"\'\n")
			newName=""
			newName+="\'"+leaf.name+"\'"
			country=gisaidSamplesAll[nameDict[leaf.name]][-1].split("/")[1].replace(" ","")
			isUK=1
			isUK2=4
			if leaf.name in IDs2:
				#newName+="[&type=\"UK_surveillance_"+gisaidSamples[v][nameDict[leaf.name]][0]+"\","
				newName+="[&type=\"UK_surveillance\","
				UKsamples+=1
				if gisaidSamplesAll[nameDict[leaf.name]][0]=="K":
					UKsamplesk+=1
			elif country=="UnitedKingdom":
				newName+="[&type=\"UK\","
				nonSurvUK+=1
				if gisaidSamplesAll[nameDict[leaf.name]][0]=="K":
					nonSurvUKk+=1
			else:
				newName+="[&type=\"non-UK\","
				nonUK+=1
				if gisaidSamplesAll[nameDict[leaf.name]][0]=="K":
					nonUKk+=1
					#print(leaf.name)
					#print(country)
				isUK=0
				isUK2=3
			newName+="eek=\""+gisaidSamplesAll[nameDict[leaf.name]][0]+"\","
			newName+="cladeNum="+str(leafDict[leaf.name].cladeNum)+","
			
			#if gisaidSamples[v][nameDict[leaf.name]][0]=="K":
			#	print(leaf.name)
			#	print(country)
			#	newName+="mutation=1,"
			#else:
			#	newName+="mutation=0,"
			newName+="isUK="+str(isUK)+","
			newName+="isUK2="+str(isUK2)+","
			newName+="location=\""+country+"\"]"
			t1String=t1String.replace(leaf.name,newName)
		print("UK samples from surveillance: "+str(UKsamples)+" of which K: "+str(UKsamplesk))
		print("UK samples not from surveillance: "+str(nonSurvUK)+" of which K: "+str(nonSurvUKk))
		print("non-UK samples: "+str(nonUK)+" of which K: "+str(nonUKk))
		file.write("	;\n")
		file.write("end;\n")
		file.write("\n")
		file.write("begin trees;\n")
		file.write("tree TREE1 = [&R] [&R=true]")
		file.write(t1String)
		file.write("\n"+"end;\n")
		file.close()










exit()




