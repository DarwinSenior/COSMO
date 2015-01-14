from numpy import array, diff, linspace, searchsorted, reshape, arange
from numpy import ones, array_str, zeros

def decomposeRange(low, high):
	i = 0
	while high!=low:
		if ((low>>i)&1):
			yield(low, low+(1<<i))
			low += 1<<i
		if ((high>>i)&1) and high!=low:
			yield(high-(1<<i), high)
			high -= 1<<i
		i += 1

def sortX(data):
	# assume one-d vector
	data.view("%s,%s"%(data.dtype, data.dtype)).sort()

def splitXY(data):
	# return two vectors of x and y
	newdata = reshape(data, (data.size/2, 2))
	return newdata[:,0], newdata[:,1]

def binRange(minval, maxval, levels):
	expos = linspace(0, 1, levels)
	return minval*((maxval/minval)**expos)

def countLevels(arr, value, bin_range):
	leftcompares = value-bin_range[::-1]
	rightcompares = value+bin_range
	leftresult = searchsorted(arr, leftcompares, side="left")
	rightresult = searchsorted(arr, rightcompares, side="right")
	rightresult[0] = leftresult[-1] = rightresult[1]
	return leftresult, rightresult

def toLevelStores(arr):
	n = (arr.size-1).bit_length()+1
	total = (ones((n, 1))*reshape(arr, (1, arr.size)))
	for i in range(1, n):
		stride = 1<<i
		for x in arange(0, arr.size, stride):
			total[i][x:x+stride].sort()
	return total

def countY(counter_y, low, high, val, bin_range, storesY):
	for l, h in decomposeRange(low, high):
		leftreuslt, rightresult = countLevels(storesY[(h-l).bit_length()-1][l:h], val, bin_range)
		reuslt = diff(leftreuslt)[::-1]+diff(rightresult)
		counter_y += reuslt

def countX(counter, valX, valY, bin_range, storesX, storesY):
	left, right = countLevels(storesX, valX, bin_range)
	for i in range(bin_range.size-1):
		countY(counter[bin_range.size-i-2], left[i], left[i+1], valY, bin_range, storesY)
		countY(counter[i], right[i], right[i+1], valY, bin_range, storesY)

def bincount(data, bin_range):
	sortX(data)
	storesX, storesY = splitXY(data)
	storesY = toLevelStores(storesY)
	counter = zeros((bin_range.size-1, bin_range.size-1), dtype="int32")
	print "sort complete"
	for (xval, yval) in data.view("%s,%s"%(data.dtype, data.dtype)):
		countX(counter, xval, yval, bin_range, storesX, storesY)
	return counter

def distanceParser(input_name):
	"""
	The input file is in following formating for each line:
	<longitute> <latitude> <z>
	<z> value for this task is alway 1 so we ignore it
	it reads the angle form and we should return radian form
	"""
	data = list()
	with open(input_name) as inputfile:
		for line in inputfile:
			ra, dec, _ = map(float, line.split())
			data.append(ra)
			data.append(dec)
	return array(data, dtype="float64")

if __name__ == '__main__':
	import sys

	data = distanceParser(sys.argv[1])
	level = int(sys.argv[2])
	bin_range = binRange(2/3600.0, 1.0, level)
	counter = bincount(data, bin_range)
	print array_str(counter)
	print sum(sum(counter))
