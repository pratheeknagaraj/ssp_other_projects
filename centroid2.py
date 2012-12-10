'''
Pratheek Nagaraj
July 1, 2011
Image Processing and Centriod Project
2. Programming Actvity Part (b)

This program will take an input file for the data array and then calculate the centroid
'''

#Import Packages
from numpy import *
from math import *

def calcCentroid( array ):
    #Total Sum
    'Use Python command to sum all elements in the array'
    sumArray = sum(array);
    
    #X Coordinate
    'Create the sum of the X Weighted Values'
    xSum = 0.0

    'Loop through the array and weight the X Values'
    for element in array:
        xSum = xSum + -1*element[0] + 0*element[1] + 1*element[2] 

    'Find the X Coordinate by dividing by the total sum'
    xCoor = xSum / sumArray
    xCoor = round(xCoor,3)

    #Y Coordinate
    'Create the sum of the Y Weighted Values'
    ySum = 0.0

    'Loop through the array and weight the Y Values'
    weight = -1.0
    for element in array:
        ySum = ySum + weight*sum(element)
        weight = weight + 1

    'Find the Y Coordinate by dividing by the total sum'
    yCoor = ySum / sumArray
    yCoor = round(yCoor,3)

    return xCoor, yCoor

#Input
'Open the containing file'
fileIn = open("input.txt", "r")
fileOut = open("output.txt", "w")

#Create Array
'Create the array with the provided data'
array = zeros((3,3))

#Fill the Array
'Position variables in file'
row = 0
column = 0

'Loop through the input file and write data to the array'
for line in fileIn:
    column = 0
    'Createa an array of elements on each line'
    newLine = line.split()
    'Loop through and insert each element into the array'
    for element in newLine:
        array[ row, column ] = float( element )
        column = column + 1
    row = row + 1

#Find the Centroid
'Call the function defined above and output centroid'
coor = calcCentroid(array)

#Output
print "The Centroid of the data given is", coor
fileOut.write("The Centroid of the data given is " + str(coor))
fileOut.close()

