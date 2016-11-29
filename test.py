from Bio import SeqIO
#from numpy import genfromtxt
import csv

# handle = open("mm8Build36/chr1.fa", "rU")
# records = list(SeqIO.parse(handle, "fasta"))
# handle.close()

# sequence = (records[0].seq)

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

line = []
chromosomeList = []
#[chromosome number, region start, region end, length, peak position]
with open('cPgList.csv', 'rb') as csvfile:
	reader = csv.reader(csvfile)
	for row in reader: 
		line = list(row[i] for i in [0,1,2,3,4])
		chromosomeList.append(line)

cOne, cTwo, cThree, cFour, cFive, cSix, cSeven, cEight = []
cNine, cTen, cEleven, cTwelve, cThirteen, cFourteen, cFifteen = []
cSixteen, cSeventeen, cEightteen, cNineteen, cX = []
for i in range(chromosomeList):
	if (chromosomeList[i][0] == '1'):
		cOne.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '2'):
		cTwo.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '3'):
		cThree.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '4'):
		cFour.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '5'):
		cFive.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '6'):
		cSix.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '7'):
		cSeven.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '8'):
		cEight.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '9'):
		cNine.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '10'):
		cTen.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '11'):
		cEleven.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '12'):
		cTwelve.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '13'):
		cThirteen.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '14'):
		cFourteen.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '15'):
		cFifteen.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '16'):
		cSixteen.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '17'):
		cSeventeen.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '18'):
		cEightteen.append(chromosomeList[i])
	if else (chromosomeList[i][0] == '4'):
		cNineteen.append(chromosomeList[i])
	else:
		cX.append(chromosomeList[i])

print cOne
		
		#while(row == 1 in reader):
		#print ', '.join(row)