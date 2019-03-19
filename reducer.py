#!/usr/bin/python
import sys
import os

def process(argc, argv):
	if (argc < 3):
		print "Usage: reducer.py name.prematrix name.prepartition name.prephy"
		return

	matrixFile = argv[1]
	partitionFile = argv[2]
	phyFile = argv[3]

	#read in binary matrixFile
	matLines = [line.strip() for line in open(matrixFile) if len(line.strip()) > 0]
	numTaxa = matLines[0].split()[0]

	#find columns to be removed in order to have a comprehensive taxon (find maximum of 1s)
	maximum = 0
	posZero = []
	for line in matLines:
		counter = countOnes(line)
		if (counter > maximum):
			posZero = getPositionOfZeros(line)
			maximum = counter

	#reduction only needed if no comprehensive taxon available
	if (len(posZero) > 0):
		#produce reduced matrix
		matfile = open(matrixFile.replace(".prematrix", ".matrix"), "w")
		for pos, line in enumerate(matLines):
			if (pos == 0):
				matfile.write(numTaxa + " " + str(maximum) + "\n")
			else:
				newline = ""
				for index, char in enumerate(line.split(" ")):
					if (index not in posZero):
						newline += char + " "
				matfile.write(newline.strip() + "\n")

		#read in partition file and produce reduced partition file
		parLines = [line.strip() for line in open(partitionFile) if len(line.strip()) > 0]
		dif = 0
		partitionsOld = []
		parFile = open(partitionFile.replace(".prepartition", ".partition"), "w")
		newSeqLength = 0
		for index, line in enumerate(parLines):
			parRange = line.split("=")[1]
			name = line.split("=")[0].strip()
			
			#partition to keep with old numbering
			start = int(parRange.split("-")[0])
			end = int(parRange.split("-")[1])
			partition = Partition(start, end)
			if (index not in posZero):
				partitionsOld.append(partition)
				parFile.write(name + " = " + partition.sub(dif).toString() + "\n")
				newSeqLength = str(partition.sub(dif).end)
			else:
				dif += partition.length

		#read in phylip file and produce reduced phylip file
		phyLines = [line.strip() for line in open(phyFile) if len(line.strip()) > 0]
		phyFile = open(phyFile.replace(".prephy", ".phy"), "w")
		for index, line in enumerate(phyLines):
			if (index == 0):
				phyFile.write(phyLines[0].split(" ")[0] + " " + newSeqLength + "\n")
			else:
				species = line.split()[0]
				sequence = line.split()[1]
				whitespaces = ""
				newline = species + whitespaces.ljust(line.count(' '))
				for partline in partitionsOld:
					newline += sequence[partline.start:partline.end]
				phyFile.write (newline.strip() + "\n")
		print "Reduction finished!"
	#if no reduction needed only rename files
	else:
		os.rename(matrixFile, matrixFile.replace(".prematrix", ".matrix"))
		os.rename(partitionFile, partitionFile.replace(".prepartition", ".partition"))
		os.rename(phyFile, phyFile.replace(".prephy", ".phy"))
		print "No reduction needed!"


#outsourced functions
def getPositionOfZeros (line):
	result = []
	for index, char in enumerate(line.split()):
		if (char == "0"):
			result.append(index)
	return result

class Partition:
	def __init__(self, start, end):
		self.start = start - 1
		self.end = end
		self.length = self.end - self.start 
	
	def sub (self, dif):
		return Partition(self.start - dif + 1, self.end - dif)

	def toString (self):
		return str(self.start + 1) + "-" + str(self.end) 

def countOnes (line):
	arg = line.split()
	counter = 0
	for index in arg:
		if (index == "1"):
			counter += 1
	return counter

process(len(sys.argv), sys.argv)
