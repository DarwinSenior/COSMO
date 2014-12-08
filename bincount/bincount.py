from numpy import log, floor, abs, squeeze, bincount, vectorize, array_split
from numpy import ones, array, zeros, float64, ravel, minimum, maximum, add, concatenate
from math import ceil


### It is a simple classifier for distance based on logorithm magnitude
### example usage>> mybin = Bin(maxbin=36000.0, minbin=2.0, binNum=30)
###	             >> bin.distance(300.0)
class Bin(object):
	def __init__(self, maxbin=36000, minbin=2, binNum=30):
		self.offset = log(minbin)
		self.size = log(maxbin/minbin)
		self.binNum = binNum
	
	# auxiliary function so that to make sure the value is always between 0 and binNum
	def limit(self, value):
		return maximum(minimum(value, self.binNum), 0)
	
	# return the log scale of distances
	def getBin(self, distance):
		scale = ((log(distance)-self.offset)/self.size)*self.binNum
		return self.limit(scale).astype(int)


def distanceParser(input_name):
	distances = list()
	with open(input_name) as inputfile:
		for line in inputfile:
			distances.append(float(line))
	return array(distances, dtype=float64)

def distancePair(group0, group1):
	m,n = group0.size, group1.size
	group0.shape = (m, 1)
	repeat0 = ones((1, n))
	first = group0*repeat0
	group1.shape = (1, n)
	repeat1 = ones((m, 1))
	second = group1*repeat1
	return ravel(abs(first-second))

def count(distances):
    return bincount(distances)

if __name__ == '__main__':
	import sys
	distances = distanceParser(sys.argv[1])
	distace_list1 = array_split(distances, ceil(len(distances)/2048.0))
	distace_list2 = array_split(distances, ceil(len(distances)/2048.0))

	mybin = Bin(maxbin=36000, minbin=2, binNum=30)
	counter = zeros(30)
	i, j = 0,0
	logfile = open("log_counter", "w")
	for d1 in distace_list1:
		i += 1
		for d2 in distace_list2:
			j += 1
			differences = distancePair(d1, d2)
			bins_result = mybin.getBin(differences)
			count_result = count(bins_result)
			count_result = concatenate((count_result, zeros(30-len(count_result), dtype=int)))
			counter = add(counter, count_result)
			print(str(counter))
			print("finished grid %d, %d"%(i,j))
	for i in range(30):
		logfile.write("bin:%d number:%d"%(i, counter[i]))