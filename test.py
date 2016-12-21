from Bio import SeqIO
#from numpy import genfromtxt
import csv

def splitChromosomes():	
	line = []
	chromosomeList = []
	#[chromosome number, region start, region end, length, peak position]
	with open('cPgList.csv', 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader: 
			line = list(row[i] for i in [0,1,2,3,4])
			chromosomeList.append(line)
			
	global cOne 
	cOne = []
	global cTwo 
	cTwo = []
	global cThree
	cThree = []
	global cFour
	cFour = [] 
	global cFive
	cFive = [] 
	global cSix
	cSix = []
	global cSeven 
	cSeven = [] 
	global cEight
	cEight = []
	global cNine
	cNine = []
	global cTen 
	cTen = [] 
	global cEleven
	cEleven = []
	global cTwelve 
	cTwelve = [] 
	global cThirteen 
	cThirteen = []
	global cFourteen
	cFourteen = [] 
	global cFifteen
	cFifteen = []
	global cSixteen
	cSixteen = []
	global cSeventeen
	cSeventeen = []
	global cEightteen
	cEightteen = []
	global cNineteen
	cNineteen = [] 
	global cX
	cX = []
	for i in range(4664):
		if (chromosomeList[i][0] == '1'):
			cOne.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '2'):
			cTwo.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '3'):
			cThree.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '4'):
			cFour.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '5'):
			cFive.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '6'):
			cSix.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '7'):
			cSeven.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '8'):
			cEight.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '9'):
			cNine.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '10'):
			cTen.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '11'):
			cEleven.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '12'):
			cTwelve.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '13'):
			cThirteen.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '14'):
			cFourteen.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '15'):
			cFifteen.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '16'):
			cSixteen.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '17'):
			cSeventeen.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '18'):
			cEightteen.append(chromosomeList[i])
		elif (chromosomeList[i][0] == '19'):
			cNineteen.append(chromosomeList[i])
		elif (chromosomeList[i][0] == 'X'):
			cX.append(chromosomeList[i])
#End of SplitChromosomes

#To avoid global variables return each chromosome to a place in an array
#put file name in same spot as chromsome array

#Possible countHistone function here to do loop "For i in range(chromosomeLength)"
#*********************************************************************************


def readChromosome(file, chromosome):	
	handle = open(file, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	#chromosome[chromosome number, region start, region end, length, peak position]
	chromosomeLength = len(chromosome)
	sequence = (records[0].seq)
	chromosomeHistone = []
	for i in range(chromosomeLength): #chromosomeLength
		start = chromosome[i][1]
		end = chromosome[i][2]
		#tail = 0 #Back histone (Last basepair) in the chain ex: -5...5 or -4..2 depending on peakPos
		#head = 0 #Head histone (Front most basepair)
		peakPos = long(chromosome[i][4]) #Highest level of 'cg' middle of 0 histone
		
		left = int(start) - int(peakPos) #Histone length region left of peak pos
		right = int(end) - int(peakPos) #Histone length region left of peak pos
		right = right + int(peakPos)
		left = left + int(peakPos)
		
		hOne = sequence[(int(peakPos) - 100):(int(peakPos) + 100)] 
		#print hOne.count('cg')
		histoneMid = [(hOne.count('cg'))+(hOne.count('CG'))]
		#histoneMid = [(hOne.count('gc')) + (hOne.count('GC'))]
		#histoneMid = [hOne.count('c') + hOne.count('g') + hOne.count('C') +hOne.count('G')]

		
		
		histoneChainRight = []
		i = int(peakPos) + 100
		while (i <= right): 
			histone = sequence[i:i+200] #sequence of 200bp histone length
			numCG = histone.count('cg') + histone.count('CG')
			#numGC = histone.count('gc') + histone.count('GC')
			#numCG = histone.count('c') + histone.count('g') + histone.count('C') +histone.count('G')
			# print histone
			# print "Right side"
			# print "Histone Lenght: " , len(histone)
			# print "Number of cg and CG: " , numCG
			# print "***************************************"

			histoneChainRight.append(numCG)
			i = i+200
			
		#print histoneChain
		
		histoneChainLeft = []
		i = int(peakPos) - 100
		while(i >= left):
			histone = sequence[i-200: i]
			numCG = histone.count('cg') + histone.count('CG')
			#numGC = histone.count('gc') + histone.count('GC')
			#numCG = histone.count('c') + histone.count('g') + histone.count('C') +histone.count('G')
			# print histone
			# print "Left side"
			# print "Histone Lenght: " , len(histone)
			# print "Number of cg and CG: " , numCG
			# print "***************************************"
			histoneChainLeft.append(numCG)
			i = i - 200
		
		#print histoneChain
		histoneChain = [histoneChainLeft,histoneMid,histoneChainRight]
		
		chromosomeHistone.append(histoneChain)
		
	print chromosomeHistone
#End readChromosome


		
#Start main
splitChromosomes()
#print cOne

genome = [cOne,cTwo,cThree,cFour,cFive,cSix,cSeven,cEight,cNine,cTen,cEleven,cTwelve,cThirteen,cFourteen,cFifteen,cSixteen,cSeventeen,cEightteen,cNineteen,cX]

chromosomeFiles = ["mm8Build36/chr1.fa", "mm8Build36/chr2.fa","mm8Build36/chr3.fa","mm8Build36/chr4.fa","mm8Build36/chr5.fa","mm8Build36/chr6.fa","mm8Build36/chr7.fa","mm8Build36/chr8.fa","mm8Build36/chr9.fa","mm8Build36/chr10.fa","mm8Build36/chr11.fa","mm8Build36/chr12.fa","mm8Build36/chr13.fa","mm8Build36/chr14.fa","mm8Build36/chr15.fa","mm8Build36/chr16.fa","mm8Build36/chr17.fa","mm8Build36/chr18.fa","mm8Build36/chr19.fa","mm8Build36/chrX.fa" ]


for i in range(20):
	print "Chromosome: ", (i+1)
	readChromosome(chromosomeFiles[i],genome[i])	
	print "********************************************************************************"

