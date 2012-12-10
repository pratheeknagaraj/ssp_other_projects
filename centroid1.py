'''
Pratheek Nagaraj
July 1, 2011
Image Processing and Centriod Project
2. Programming Actvity Part (a)

This program will determine the centroid of the data given in the table provided
'''

#Import Packages
from numpy import *
from math import *

#Create Array
'Create the array with the provided data'
array = ([[0, 33, 21, 33, 8],
          [0, 56, 51, 53, 26],
          [23, 120, 149, 73, 18],
          [55, 101, 116, 50, 16],
          [11, 78, 26, 2, 10]])

#Total Sum
'Use Python command to sum all elements in the array'
sumArray = sum(array);

#X Coordinate
'Create the sum of the X Weighted Values'
xSum = 0.0

'Loop through the array and weight the X Values'
for element in array:
    xSum = xSum + -2*element[0] + -1*element[1] + 0*element[2] + 1*element[3] + 2*element[4]

'Find the X Coordinate by dividing by the total sum'
xCoor = xSum / sumArray
xCoor = round(xCoor,3)

#Y Coordinate
'Create the sum of the Y Weighted Values'
ySum = 0.0

'Loop through the array and weight the Y Values'
weight = -2.0
for element in array:
    ySum = ySum + weight*sum(element)
    weight = weight + 1

'Find the Y Coordinate by dividing by the total sum'
yCoor = ySum / sumArray
yCoor = round(yCoor,3)

#Output
print "The Centroid of the data given is (", xCoor, ", ", yCoor, ")"
