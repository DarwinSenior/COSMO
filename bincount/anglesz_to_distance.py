from COSMO.common.pczsource import pczsource

def parseLineToDistance(line):
	theta, phi, z = map(float, line.split(" "))
	source = pczsource.create(theta, phi, z)
	return source.pcDistance()
"""
coverting the formatting 
theta, phi, z 
to
distance
"""
def anglez_to_distance(input_name, output_name):
	with open(input_name, "r") as input_file:
		with open(output_name, "w") as output_file:
			for line in input_file:
				output_file.write("%f\n"%parseLineToDistance(line))
	input_file.close()
	output_file.close()

# example
if __name__ == '__main__':
	anglezt_to_distance("../theta_phi_z.in", "../distance.out")
	