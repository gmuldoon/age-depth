def cleanDatFile(inputFileName, outputFileName):
	''' Format the input data to be easily readable by pandas'''
	import re
	input = open(inputFileName, "r")
	output = open(outputFileName, "w")

	# read only data lines (skip all header-type lines at top and within data)
	for line in input:
		if line.lstrip().startswith('!'): 
			output.write(line)
		elif line.lstrip().startswith("lm"):
			line = re.sub(' +',' ',line)
			line = line.replace(" ", ",")
			output.write(line)

	input.close()
	output.close()

def trimDatFile4Paraview(filein,fileout):
	import pandas as pd

	# Read data as is
	data=pd.read_csv(filein,sep=',',skiprows=48,names=['name','x','y','z','trace','PST'],index_col=False,usecols=['x','y','z'])

	# Write the data back out with only columns x,y,z and no Landmark header info
	data.to_csv(fileout,sep=',',index=False)




