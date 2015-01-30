import json
import numpy
from numpy import searchsorted, ravel, array, abs, ceil
from numpy import float64, zeros, ones, linspace, array_split
import base64

def split_data(ras, decs):
	"""
	It will split the RAs and DECs into smaller chunks which would be better
	for cache coherent
	"""
	size = ceil(len(ras)/256.0)
	return zip(array_split(ras, size), array_split(decs, size))

def addCount(count, value0, value1, bin_range):
	x = searchsorted(bin_range, value0)
	y = searchsorted(bin_range, value1)
	# print x, y, bin_range, value1, value0
	count[x][y] += 1

def binRange(minval, maxval, levels):
	expos = linspace(0, 1, levels)[1:]
	return (minval*((maxval/minval)**expos))

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


def pointsPairsToJson(ra1, dec1, ra2, dec2):
	data = dict()
	data["ra1"] = base64.b64encode(ra1)
	data["ra2"] = base64.b64encode(ra2)
	data["dec1"] = base64.b64encode(dec1)
	data["dec2"] = base64.b64encode(dec2)
	return json.dumps(data)


def mapper(input_file, output_file):
	RAs, DECs = distanceParser(input_file)
	groups = split_data(RAs, DECs)
	for group1 in groups:
		for group2 in groups:
			data = pointsPairsToJson(group1[0], group1[1], group2[0], group2[1])
			output_file.write(data)
			output_file.write("\n")

def reducer(input_file, output_file):
	bin_range = binRange(2/3600.0, 1.0, 5)
	counter = zeros((bin_range.size, bin_range.size), dtype="int64")
	addcount = lambda x,y:addCount(counter, x, y, bin_range)
	for line in input_file:
		data = json.loads(line)
		ra1 = numpy.frombuffer(base64.decodestring(data["ra1"]))
		ra2 = numpy.frombuffer(base64.decodestring(data["ra2"]))
		dec1 = numpy.frombuffer(base64.decodestring(data["dec1"]))
		dec2 = numpy.frombuffer(base64.decodestring(data["dec2"]))
		ras = distance(ra1, ra2)
		decs = distance(dec1, dec2)
		map(addcount, ras, decs)
	counter = ravel(counter)
	output_file.write(base64.b64encode(counter))
	output_file.write("\n")

def combiner(input_file, output_file):
	bin_range = binRange(2/3600.0, 1.0, 5)
	counter = zeros(bin_range.size*bin_range.size, dtype="int64")
	for line in input_file:
		data = numpy.frombuffer(base64.decodestring(line), dtype="int64")
		counter += data
	counter.shape = (bin_range.size, bin_range.size)
	output_file.write(str(counter.tolist()))



