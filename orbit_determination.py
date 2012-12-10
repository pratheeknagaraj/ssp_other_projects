'''
Pratheek Nagaraj
July 18, 2011
Orbit Determination

Find the orbital elements using a series of observations
'''

##### Import #####

from visual import *
from math import *
from numpy import *

##### Constants #####

k = 0.01720209895
epsilon = 23.43799 #degrees

##### Functions #####

def convertHMStoDegDec( hours, minutes, seconds ):
     if hours < 0:
          return (hours - minutes/60.0 - seconds/3600.0) * 15
     else:
          return (hours + minutes/60.0 + seconds/3600.0) * 15

def convertDMStoDegDec( degrees, minutes, seconds ):
     if degrees < 0:
          return degrees - minutes/60.0 - seconds/3600.0
     else:
          return degrees + minutes/60.0 + seconds/3600.0

def convertDegDectoHMS( degrees, default = "HMS" ):
     degrees = degrees / 15.0
     
     if degrees > 0:
          hours = int(degrees)
     else:
          hours = int(degrees)
          degrees = abs(degrees)
     minutes = int(degrees*60%60)
     seconds = round(degrees*3600%3600%60,4)

     if default == "HMS":
          return hours, minutes, seconds
     elif default == "MS":
          return minutes + 60*hours, seconds
     elif default == "S":
          return seconds + minutes*60 + hours*3600

def convertDegDectoDMS( degreesDec, default = "DMS" ):
     
     if degrees > 0:
          degrees = int(degreesDec)
     else:
          degrees = int(degreesDec)
          degreesDec = abs(degreesDec)
     minutes = int(degreesDec*60%60)
     seconds = round(degreesDec*3600%3600%60,4)

     if default == "DMS":
          return degrees, minutes, seconds
     elif default == "MS":
          return minutes + 60*degrees, seconds
     elif default == "S":
          return seconds + minutes*60 + degrees*3600

def toRadians( deg ):
     return deg * pi / 180

def toDegrees( rad ):
     return rad * 180 / pi

def listToVector( listIn ):
     vectorOut = vector(listIn[0],listIn[1],listIn[2])
     return vectorOut

def readInput():
     # !!!
     #filePath = raw_input("Please enter the file name (without .txt): ")
     #fileIn = open( filePath + ".txt", "r" )
     #print ""
     fileIn = open( "data2.txt", "r" )

     ra = []
     dec = []
     time = []
     obs = 0

     for line in fileIn:
          obs = obs + 1
          column = 0
          newLine = line.split()

          ra.append( convertHMStoDegDec( float(newLine[0]), float(newLine[1]), float(newLine[2]) ) )
          dec.append( convertDMStoDegDec( float(newLine[3]), float(newLine[4]), float(newLine[5]) ) )
          time.append( float(newLine[6]) )

     return ra, dec, time, obs

def getPHats( ra, dec ):
     phat = []

     for alpha, delta in zip(ra, dec):
          phat.append( vector( cos(toRadians(alpha))*cos(toRadians(delta)), \
                               sin(toRadians(alpha))*cos(toRadians(delta)), sin(toRadians(delta)) ) )

     return phat

def getTaus( time ):
     taus = []
     for i in range( 0, len(time) ):
          if i != len(time)/2:
               if i >  len(time)/2:
                    value = k*( time[i] - time[len(time)/2] )
                    taus.append( value )
               else:
                    value = k*( time[i] - time[len(time)/2] )
                    taus.append( value )

     return taus

# !!!!! Not Completed Function !!!!! #
#------------------------------------#
def getBigRVector( time ):
     RVector = [ vector(-0.35795,0.87291,0.37843), vector(-0.46601,0.82851,0.35918), vector(-0.63547,0.72609,0.31478) ]
     return RVector
#------------------------------------#

def tripleScalarProduct( vector1, vector2, vector3 ):
     return dotProduct( crossProduct( vector1, vector2 ), vector3 )

def dotProduct( vector1, vector2 ):
    #Calculate the dot product of a given vector
    if len(vector1) != len(vector2):
        return "Vectors are not same dimensions, please enter the same dimension vectors and try again.\nGoodbye!"
    
    dotValue = 0.0
    for i in range( 0, len(vector1) ):
        dotValue = dotValue + vector1[i]*vector2[i]

    return dotValue

def crossProduct( vector1, vector2 ):
    #Calculate the cross product of a given vector, user choice to display
    if len(vector1) != 3 or len(vector2) != 3:
        return "Vectors must be three dimensional and the same number of dimensions, please try again.\nGoodbye!"

    matrix = zeros(( 3,3 ))
    matrix[0] = (1,1,1)
    matrix[1] = vector1
    matrix[2] = vector2

    vectorNew = ( determinant( minorArray( matrix, 0 ) ), -determinant( minorArray( matrix, 1 ) ), determinant( minorArray( matrix, 2 ) ) )

    return vectorNew

def minorArray( inArray, column ):
    #Get Minor Array
    front = inArray[1:,0:column]
    back = inArray[1:,column+1:]
    newarray = append(front,back,axis = 1)
    return newarray

def determinant( inArray ):
    #Get the Determinant
    det = 0
    for i in range(0,len(inArray)):
        if len(inArray) == 2:
            det = inArray[0,0]*inArray[1,1] - inArray[1,0]*inArray[0,1]
            return det
        else:
            det = det + inArray[0,i]*(-1)**i * determinant( minorArray( inArray, i ) )
    
    return det

def calculateRVectors( RVector, p, phat ):
     rvector = []
     for i in range( 0, len(phat) ):
          rvector.append( p[i]*phat[i] - RVector[i] )
     return rvector

def calculateFandGSmall( rvector, taus ):
     f = []
     g = []
     middle = len(rvector)/2
     for i in range(0, len(taus)):
          fValue = 1. - (1./(2*((mag(rvector[middle]))**3)))*((taus[i])**2)
          gValue = taus[i] - (1./(6*((mag(rvector[middle]))**3)))*((taus[i])**2)
          f.append(fValue)
          g.append(gValue)
     return f, g

def calculateFandGLarge( rvector, rdotvector, taus ):
     f = []
     g = []
     middle = len(rvector)/2
     rMid = rvector[middle]
     rDMid = rdotvector
     for t in taus:
          fValue = (1. - (1./(2*((mag(rMid))**3)))*(t**2) + ((dot( rMid, rDMid ))/(2.0*(mag(rMid)**5)))*t**3 \
                   + (1./24.)*(((3.0)/((mag(rMid))**3))*(dot(rMid,rDMid)/(mag(rMid)**2)-(1./(mag(rMid))**3)) \
                   - ((15.*(dot(rMid,rDMid))**2)/(mag(rMid)**7))+(1./(mag(rMid)**6)))*(t**4))
          gValue = t - (1./(6*((mag(rMid)**3))))*(t**3)+((dot(rMid,rDMid))/(4.*(mag(rMid))**5))*(t**4)
          f.append(fValue)
          g.append(gValue)
     return f, g

def recalculateA( f, g ):
     base = (f[0]*g[1]-f[1]*g[0])*1.0
     a1 = (g[1])/base
     a3 = -(g[0])/base
     return a1, a3

def calculateRDotVector( f, g, rvector ):
     base = g[0]*f[1]-g[1]*f[0]
     rdotvector = ((f[1]*1.0)/(base))*rvector[0] - ((f[0]*1.0)/(base))*rvector[2]
     return rdotvector

def orbitalElements( rvector, rdotvector ):
     rvector = listToVector( rvector[1] )
     rdotvector = listToVector( rdotvector )
     rvector = rotate(rvector, angle = toRadians(-epsilon), axis = (1,0,0) )
     rdotvector = rotate(rdotvector, angle = toRadians(-epsilon), axis = (1,0,0) )
     mu = 1
     t = time[1]
     magr = mag( rvector )
     magrdot = mag( rdotvector )

     a = 1/((2/magr)-((magrdot**2)/mu))
     hvector = cross( rvector, rdotvector )
     
     magh = mag( hvector )
     
     rdotcrossh = listToVector( cross( rdotvector, hvector ) )
     evector = ((rdotcrossh)/mu) - (rvector/magr)
     e = mag( evector )
     
     kvector = vector(0,0,1)
     hvectordotk = dot( kvector, hvector )
     
     i = toDegrees( acos( hvectordotk/magh ) )
     
     nvector = listToVector( cross( kvector, hvector ) )
     magn = mag( nvector )
     ivector = (1,0,0)
     idotn = dot( ivector, nvector )
     omega = toDegrees( acos( idotn/magn ) )
     if nvector.y < 0:
         omega = 360 - omega
     edotn = dot( evector, nvector )
     littleOmega = toDegrees( acos( ( dot( evector, nvector ) )/((e)*(magn)) ) )
     if evector.z < 0:
         littleOmega = 360 - littleOmega
     E = toDegrees(acos((1-(magr/a))/e))
     if dotProduct(rvector, rdotvector) < 0:
         E = 360 - E
     T = t - ((toRadians(E) - e*sin(toRadians(E)))/(k*sqrt((mu/(a**3)))))
     M = (toRadians(E) - e*sin(toRadians(E)))%(2*pi)

     print "Semimajor Axis: ", a
     print "Eccentricity: ", e
     print "Inclination: ", i
     print "Longitude of Ascending Node: ", omega
     print "Argument of Perihelion: ", littleOmega
     print "Time of perihelion passage: ", T
     print "Mean anomaly: ", M

#- Input Data -#

#choice = raw_input("Would you like to import text data (Y/N)?: ")
# !!!
choice = "y"
ra = []
dec = []
time = []
obs = 0
    
if choice == "Y" or choice == "y" or choice == "Yes" or choice == "yes":    
     ra, dec, time, obs = readInput()
else:
     print "Tough Luck"

#- Get p hat, taus, eath sun vector -#

phat = getPHats( ra, dec )
taus = getTaus( time )
RVector = getBigRVector( time )

#- Triple Products -#

R1xp2dotp3 = tripleScalarProduct( RVector[0], phat[1], phat[2] )
R2xp2dotp3 = tripleScalarProduct( RVector[1], phat[1], phat[2] )
R3xp2dotp3 = tripleScalarProduct( RVector[2], phat[1], phat[2] )
p1xR1dotp3 = tripleScalarProduct( phat[0], RVector[0], phat[2] )
p1xR2dotp3 = tripleScalarProduct( phat[0], RVector[1], phat[2] )
p1xR3dotp3 = tripleScalarProduct( phat[0], RVector[2], phat[2] )
p2xR1dotp1 = tripleScalarProduct( phat[1], RVector[0], phat[0] )
p2xR2dotp1 = tripleScalarProduct( phat[1], RVector[1], phat[0] )
p2xR3dotp1 = tripleScalarProduct( phat[1], RVector[2], phat[0] )
p1xp2dotp3 = tripleScalarProduct( phat[0], phat[1], phat[2] )
p2xp3dotp1 = tripleScalarProduct( phat[1], phat[2], phat[0] )

#- Initial Guesses for p -#

a1 = ( time[2] - time[1] )/( time[2] - time[0] )
a3 = -( time[0] - time[1] )/( time[2] - time[0] )

#- Calculate p1, p2, p3 -#

p1 = ( a1*R1xp2dotp3 - R2xp2dotp3 + a3*R3xp2dotp3 )/( a1*p1xp2dotp3 )
p2 = ( a1*p1xR1dotp3 - p1xR2dotp3 + a3*p1xR3dotp3 )/( -1*p1xp2dotp3 )
p3 = ( a1*p2xR1dotp1 - p2xR2dotp1 + a3*p2xR3dotp1 )/( a3*p2xp3dotp1 )
p = [p1,p2,p3]

#- Calculate r -#

rvector = calculateRVectors( RVector, p, phat )

#- F and G series -#

f, g = calculateFandGSmall( rvector, taus )
print "Start f and g"
print f, g

#- Recalculate a1 and a3 -#

a1, a3 = recalculateA( f, g )
print "Recalculated a1 and a3"
print a1, a3

#- Full F and G -#

p1 = ( a1*R1xp2dotp3 - R2xp2dotp3 + a3*R3xp2dotp3 )/( a1*p1xp2dotp3 )
p2 = ( a1*p1xR1dotp3 - p1xR2dotp3 + a3*p1xR3dotp3 )/( -1*p1xp2dotp3 )
p3 = ( a1*p2xR1dotp1 - p2xR2dotp1 + a3*p2xR3dotp1 )/( a3*p2xp3dotp1 )
p = [p1,p2,p3]

rvector = calculateRVectors( RVector, p, phat )
rdotvector = calculateRDotVector( f, g, rvector )
print rvector
print rdotvector

#- Main Iteration Loop -#
cont = True
tolerance = 1e-3
while cont == True:
     pastR = rvector[1]
     f, g = calculateFandGLarge( rvector, rdotvector, taus )
     print "f and g"
     print f, g
     a1, a3 = recalculateA(f, g)
     print "a1 and a3"
     print a1, a3
     p1 = ( a1*R1xp2dotp3 - R2xp2dotp3 + a3*R3xp2dotp3 )/( a1*p1xp2dotp3 )
     p2 = ( a1*p1xR1dotp3 - p1xR2dotp3 + a3*p1xR3dotp3 )/( -1*p1xp2dotp3 )
     p3 = ( a1*p2xR1dotp1 - p2xR2dotp1 + a3*p2xR3dotp1 )/( a3*p2xp3dotp1 )
     p = [p1,p2,p3]
     rvector = calculateRVectors( RVector, p, phat )
     rdotvector = calculateRDotVector( f, g, rvector )
     difference = [fabs(a - b) for a, b in zip(pastR, rvector[1])]
     if difference[0] < tolerance and difference[1] < tolerance and difference[2] < tolerance:
          cont = False
          
#- Convert to orbital elements -#

print rvector
print rdotvector

orbitalElements( rvector, rdotvector )

