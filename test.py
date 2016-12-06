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
		elif (chromosomeList[i][0] == '4'):
			cNineteen.append(chromosomeList[i])
		else:
			cX.append(chromosomeList[i])

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
	for i in range(20): #chromosomeLength
		start = chromosome[i][1]
		end = chromosome[i][2]	
		#tail = 0 #Back histone (Last basepair) in the chain ex: -5...5 or -4..2 depending on peakPos
		#head = 0 #Head histone (Front most basepair)
		peakPos = chromosome[i][4] #Highest level of 'cg' middle of 0 histone
		
		left = int(start) - int(peakPos) #Histone length region left of peak pos
		right = int(end) - int(peakPos) #Histone length region left of peak pos
		right = right + int(peakPos)
		left = left + int(peakPos)
		
		hOne = sequence[(int(peakPos) - 100):(int(peakPos) + 100)] 
		#print hOne.count('cg')
		histoneChain=[hOne.count('cg')]
		
		i = int(peakPos) + 100
		while (i <= right): 
			Histone = sequence[i:i+200] #sequence of 200bp
			numCG = Histone.count('cg')
			histoneChain.append(numCG)
			i = i+200
			
		#print histoneChain
		
		i = int(peakPos) - 100
		while(left <= i):
			Histone = sequence[left: left + 200]
			numCG = Histone.count('cg')
			histoneChain.append(numCG)
			left = left + 200
		
		#print histoneChain

		chromosomeHistone.append(histoneChain)
		
	print chromosomeHistone
#End readChromosome
	
	

# start = 3967125
# end = 3968125
# peakPos = 3661437

# firstSeq = sequence[start:end]

# numC = firstSeq.count('c') + firstSeq.count('C')

# numG = firstSeq.count('g') + firstSeq.count('G')


# print firstSeq

# print "Number of C:", numC
# print "Number of G:", numG

#my_data = genfromtxt('cPgList.csv', dilimiter = ',')
#print my_data[1][1]			
			
			
		
		
#Start main
splitChromosomes()
#print cOne

readChromosome("mm8Build36/chr1.fa",cOne)				

