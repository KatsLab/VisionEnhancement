import sys
import scipy.io as sio
import numpy as np
"""
Constants below
"""
TOLERATION = 0.5 #toleration, a distance range, within which two colors are considered as the same color
META_DIFF = 2.3 #distance in LAB, within which are undistinguishable/potential metamers


minorMatrix = []
retlist = []


"""
This function reads data in a*b*3 form.
It takes the rSize by cSize matrix centered
at x,y and load then into minorMatrix with
each pixel's location info.

fName - file name
x,y - center of matrix
rSize - num row of small taken matrix
cSize - num col of small taken matrix
"""
#TODO: edge checking still to be implemented
def loadMinorMatrix(fName, x, y, rSize, cSize):
	global minorMatrix 
	wholeMatrix = sio.loadmat(fName)
	x_orig = x - int(rSize/2)
	y_orig = y - int(cSize/2)
	temp = wholeMatrix['lab'][x_orig:x_orig + rSize, y_orig:y_orig + cSize]
	for i in range(rSize):
		for j in range(cSize):
			minorMatrix.append(list(np.append(temp[i,j],[x_orig + i,y_orig + j])))

	

"""
This function calculates two pixels' distance in terms of 
LAB distance.
"""
def distance(lab1, lab2):
	distance = (lab1[0]-lab2[0])*(lab1[0]-lab2[0])+\
			   (lab1[1]-lab2[1])*(lab1[1]-lab2[1])+\
			   (lab1[2]-lab2[2])*(lab1[2]-lab2[2])
	return np.sqrt(distance)

	
"""
This function brute forcely examine pixels pairwise
and add all metamers pairs into the retlist, each
pair are in form pixel1_location, pixel2_location, distance in LAB

pixelList - list of pixels in form of [l, a, b, x_loc, y_loc]
toleration - distance within which two colors are the same color
maxDiff - distance within which two colors are potential metamers
"""
def bruteForce(pixelList, toleration, maxDiff):
	global retlist
	for i in range(len(pixelList)):
		for j in range(i+1,len(pixelList)):
			lab1 = pixelList[i][:3]
			lab2 = pixelList[j][:3]
			dis = distance(lab1, lab2)
			if dis > toleration and dis < maxDiff:
				retlist.append([pixelList[i][3:],pixelList[j][3:],dis])
	

"""
This funtion finds all pixels in pixelList, such that the pixel
lies in a sphere, centered at (l,a,b) with diameter of maxDiff

pixelList - see bruteForce
maxDiff - diameter of sphere
l,a,b - center of the sphere
"""
def pixelCenterAt(pixelList, maxDiff, l, a, b):
	inSphere = []
	for i in range(len(pixelList)):
		if distance([l,a,b],pixelList[i][:3]) < maxDiff/2 :
			inSphere.append(pixelList[i])
	return inSphere

"""
This function finds all metamers pairs around a LAB color
"""
def metaCenterAt(pixelList, toleration, maxDiff, l, a, b):
	inSphere = pixelCenterAt(pixelList, maxDiff, l, a, b)
	bruteForce(inSphere, toleration, maxDiff)
	print 'number of pixels lies in Sphere: '+str(len(inSphere))
	
"""
This function roughly catagorize pixelList based on their position in
LAB space into a dictionary. The LAB space is divided into little
cubes, with edge length of toleration. Points in the same cube
are catagorized into the same key.

centers of cubes in LAB - (a*toleration+offset
						   b*toleration+offset
						   c*toleration+offset)
pixelList - list of pixels in form of [l, a, b, x_loc, y_loc]
toleration - edge length of cubes, also a standard for cube centers
offset - a shift parameter of cube centers, in range of [0,toleration)

key format of returned dictionary - (l_center, a_center, b_center)

For example: toleration = 1, offset = 0.75
			centers of cubes should be (a*1+0.75, b*1+0.75, c*1+0.75)
			e.g. (4.75,-3.25, 27.75)
"""
def rawCatagorize(pixelList, toleration, offset):
	if offset>=toleration or offset<0:
		print >> sys.stderr, 'offset should be in range of [0,toleration)'
		sys.exit(1)
	locCatagory = {}
	for i in range(len(pixelList)):
		temp = pixelList[i][:3]
		for j in range(3):
			temp[j] = np.floor(((temp[j]-offset)+ toleration/2.0)/toleration)*toleration
		key = tuple(temp)
		if key not in locCatagory:
			locCatagory[key] = []
		locCatagory[key].append(pixelList[i])
	return locCatagory
	
def main():
	loadMinorMatrix('LAB image_004 data.mat', 256, 256, 512, 512)
	dic =  rawCatagorize(minorMatrix, TOLERATION, 0)
	#we regard 0.5 as toleration for the same color
	
	#bruteForce(minorMatrix, TOLERATION, META_DIFF)
	
	#metaCenterAt(minorMatrix, TOLERATION, META_DIFF, 24, 0, -37) 
		#for good result - 24, 0, -37
	
	
	
	"""
	for key in dic.keys():
		print 'key: ' + str(key)+"  length: "+str(len(dic[key]))
	"""
	"""
	for key in dic[(8.5,1.5,-34.5)]:
		print key
	"""
	print 'total number of pixels: '+str(len(minorMatrix))
	print 'number of catagories: '+ str(len(dic.keys()))
	#print 'number of metamers pairs: '+str(len(retlist))
	
if __name__ == '__main__':
	main()
	