from numpy import searchsorted, ravel, array, abs
from numpy import float64, zeros, ones, linspace
def binRange(minval, maxval, levels):
	expos = linspace(0, 1, levels)[1:]
	return (minval*((maxval/minval)**expos))

def addCount(count, value0, value1, bin_range):
	x = searchsorted(bin_range, value0)
	y = searchsorted(bin_range, value1)
	# print x, y, bin_range, value1, value0
	count[x][y] += 1

def distanceParser(input_file):
	"""
	The input file is in following formating for each line:
	<longitute> <latitude> <z>
	<z> value for this task is alway 1 so we ignore it
	it reads the angle form and we should return radian form
	"""
	ras, decs = list(), list()
	for line in input_file:
		ra, dec, _ = map(float, line.split())
		ras.append(ra)
		decs.append(dec)
	return array(ras, dtype=float64), array(decs, dtype=float64)

def distance(group1, group2):
	m, n = group1.size, group2.size
	group1.shape = (m, 1)
	repeat1 = ones((1, n))
	d1 = group1*repeat1

	group2.shape = (1, n)
	repeat2 = ones((m, 1))
	d2 = group2*repeat2

	return ravel(abs(d1-d2))

def main():
	import sys
	RAs, DECs = distanceParser(open(sys.argv[1]))

	ras = distance(RAs, RAs)
	decs = distance(DECs, DECs)

	level = int(sys.argv[2])
	bin_range = binRange(2/3600.0, 1.0, level)
	counter = zeros((bin_range.size, bin_range.size), dtype="int64")
	map(lambda x,y:addCount(counter, x, y, bin_range), ras, decs)

	print counter
	print sum(sum(counter))

if __name__ == '__main__':
	main()


	