from numpy import log, floor, abs, bincount, array_split, int64
from numpy import sin, cos, arcsin, radians, sqrt
from numpy import ones, array, zeros, float64, ravel, minimum, maximum, add, concatenate
from math import ceil

import matplotlib.pyplot as plt


### It is a simple classifier for distance based on logorithm magnitude
### example usage>> mybin = Bin(maxbin=36000.0, minbin=2.0, binNum=30)
###	             >> bin.distance(300.0)
class Bin(object):
	def __init__(self, maxbin=36000, minbin=2, binNum=30):
		self.offset = log(minbin)
		self.size = log(maxbin/minbin)
		self.binNum = binNum
		self.counter = zeros(binNum+1, dtype=int64)
	
	# auxiliary function so that to make sure the value is always between 0 and binNum
	def limit(self, value):
		return maximum(minimum(value, self.binNum), 0)
	
	# return the log scale of distances
	def getBin(self, distance):
		scale = ((log(distance)-self.offset)/self.size)*self.binNum
		return self.limit(scale).astype(int)

	def count(self, distances):
		"""
		seperate the distances to different bins and do bin count
		then add the result to bins
		then accumulate the count result to the bin
		"""
		count_result = bincount(self.getBin(distances))
		print len(count_result)
		count_result = concatenate((count_result, zeros(self.binNum+1-len(count_result), dtype=int)))
		self.counter += count_result

	def __add__(self, other):
		"""
		Note that self and other should be of the same max,min and bin
		and it will create a newbin and add the counter together
		"""
		newbin = Bin(self.maxbin, self.minbin, self.binNum)
		newbin.counter = self.counter+other.counter


def distanceParser(input_name):
	"""
	The input file is in following formating for each line:
	<longitute> <latitude> <z>
	<z> value for this task is alway 1 so we ignore it
	it reads the angle form and we should return radian form
	"""
	ras, decs = list(), list()
	with open(input_name) as inputfile:
		for line in inputfile:
			ra, dec, _ = map(float, line.split())
			ras.append(ra)
			decs.append(dec)
	return radians(array(ras, dtype=float64)), radians(array(decs, dtype=float64))

def split_data(ras, decs):
	"""
	It will split the RAs and DECs into smaller chunks which would be better
	for cache coherent
	"""
	size = ceil(len(ras)/2048.0)
	return zip(array_split(ras, size), array_split(decs, size))

def distancePair(group0, group1):
	"""
	The group contains (RAs, DECs) pair, and the formular for calculating
	angle distance is the multi-value version of Professor Brunner's hsAngularDistance
	the source code and explaination is in pcsource.py
	"""
	ra0, dec0 = group0
	ra1, dec1 = group1

	m,n = ra0.size, ra1.size

	ra0.shape = (m, 1)
	dec0.shape = (m, 1)
	repeat0 = ones((1, n))
	cos_dec0 = cos(dec0)
	ra0 = ra0*repeat0
	dec0 = dec0*repeat0

	ra1.shape = (1, n)
	dec1.shape = (1, n)
	repeat1 = ones((m, 1))
	cos_dec1 = cos(dec1)
	ra1 = repeat1*ra1
	dec1 = repeat1*dec1

	deltaRA = 0.5*(ra0-ra1)
	deltaDEC = 0.5*(dec0-dec1)

	angles = 2*(arcsin(minimum(1.0, sqrt(sin(deltaDEC)**2)+(cos_dec0*cos_dec1)*sin(deltaRA)**2)))
	return ravel(abs(angles))


if __name__ == '__main__':
	import sys
	RAs, DECs = distanceParser(sys.argv[1])
	groups = split_data(RAs, DECs)

	mybin = Bin(maxbin=36000, minbin=2, binNum=30)
	i, j = 0,0
	try:
		logfile = open(sys.argv[2], "w")
	except:
		logfile = sys.stdout
	for group0 in groups:
		i += 1
		for group1 in groups:
			j += 1
			distances = distancePair(group0, group1)*36000
			mybin.count(distances)
			print("finished grid %d, %d"%(i,j))
	for i in range(30):
		print("bin:%d number:%d"%(i, mybin.counter[i]))
		logfile.write("%d\n"%mybin.counter[i])