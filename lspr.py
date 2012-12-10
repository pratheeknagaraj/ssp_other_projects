'''
Pratheek Nagaraj
July 8, 2011
Least Squares Plate Reduction Project

Text
'''

#Import Packages
from visual import *
from numpy import *
from math import *
import pyfits
import numdisplay

def minorArray( inArray, column ):
    #Get Minor Array
    'Split to make minor array based on column'
    front = inArray[1:,0:column]
    back = inArray[1:,column+1:]
    newarray = append(front,back,axis = 1)
    return newarray

def determinant( inArray ):
    #Get the Determinant
    det = 0
    'Loop through array for determinant'
    for i in range(0,len(inArray)):
        if len(inArray) == 2:
            'Minor is a 2x2'
            det = inArray[0,0]*inArray[1,1] - inArray[1,0]*inArray[0,1]
            return det
        else:
            'Recursively get minor'
            det = det + inArray[0,i]*(-1)**i * determinant( minorArray( inArray, i ) )
    
    return det

def cramerSlice( inArray, column ):
    #Get Minor Array
    'Split to make minor array based on column'
    front = inArray[:,0:column]
    mid = inArray[:,len(inArray):len(inArray)+1]
    back = inArray[:,column+1:]
    'Append arrays onto each other'
    frontmid = append(front,mid,axis = 1)
    newarray = append(frontmid,back,axis = 1)
    finalarray = newarray[0:len(inArray),0:len(inArray)]
    return newarray

def cramer( inArray1 , inArray2, GaussUse, TimeCalc ):

    if TimeCalc:
        start = time.time()
    
    inArray=array([])
    inArray=zeros((len(inArray2),len(inArray2)+1))
    for i in range( 0,len(inArray2) ):
        for j in range( 0, len(inArray2) ):
            inArray[i,j] = inArray1[i][j]
        inArray[i,len(inArray2)]=inArray2[i]
    #Input system of equation to do Cramer's Rule

    if GaussUse or len(inArray2) > 8:
        inArray = Gauss( inArray )
        return "blah"
    else:
        Dbase = 0.0
        'Get D for Cramers'
        Dbase = determinant( cramerSlice( inArray, len(inArray) ) )
        if Dbase == 0:
            return "The coefficient array is singular"
        
        solutionArray = array([])
        solutionArray = zeros(len(inArray))
        'Get D of sub for Cramers'
        for i in range(0,len(inArray)):
            Dsub = determinant( cramerSlice( inArray, i ) )
            solutionArray[i] = ( Dsub/Dbase )

    if TimeCalc:
        end = time.time()
        print "Cramer Evaluation Time: ", end - start
    
    return solutionArray

def Gauss( inArray ):

    for i in range(0,len(inArray)):
        pivotCol = i
        diagElement = 0.0
        diagElement = inArray[pivotCol,pivotCol]
        for j in range(i,len(inArray)+1):
            inArray[pivotCol,j] = inArray[pivotCol,j]/diagElement
        diagElement = inArray[pivotCol,pivotCol]
        for k in range(pivotCol+1,len(inArray)):
            point = 0.0
            point = inArray[k,pivotCol]
            multiple = point/diagElement
            for l in range(pivotCol,len(inArray)+1):
                inArray[k,l] = inArray[k,l] - multiple * inArray[pivotCol,l]
        print inArray
    return inArray

def LSPR():

    choice = raw_input("Would you like to import text data (Y/N)?: ")
    
    if choice == "Y" or choice == "y" or choice == "Yes" or choice == "yes":
        filePath = raw_input("Please enter the file name (without .txt): ")
        fileIn = open( filePath + ".txt", "r" )
        print ""

        ra = []
        dec = []
        xlist = []
        ylist = []
        column = 0
        num = -1

        for line in fileIn:
            num = num + 1
            column = 0
            newLine = line.split()
            
            data = []
            
            for element in newLine:
                data.append( float( element ) )

            if len(data) == 8:
                xlist.append( data[0] )
                ylist.append( data[1] )

                if data[2] < 0:
                    ra.append( data[2] - data[3]/60.0 - data[4]/3600.0 )
                else:
                    ra.append( data[2] + data[3]/60.0 + data[4]/3600.0 )
                    
                if data[5] < 0:
                    dec.append( data[5] - data[6]/60.0 - data[7]/3600.0 )
                else:
                    dec.append( data[5] + data[6]/60.0 + data[7]/3600.0 )
            else:
                xObject = data[0]
                yObject = data[1]

               
    else:
        num = input("Please enter the number of reference stars: ")
        print "\n- Please Input Data -\n"
        ra = []
        dec = []
        xlist = []
        ylist = []

        for i in range(0,num):
            xlist.append( input("Please enter the x coordinate of reference star number " + str(i+1) + ": " ) )
            ylist.append( input("Please enter the y coordinate of reference star number " + str(i+1) + ": " ) )
            
            rah = input("Please enter the RA of reference star number " + str(i+1) + " (hours): ")
            ram = input("Please enter the RA of reference star number " + str(i+1) + " (minutes): ")
            ras = input("Please enter the RA of reference star number " + str(i+1) + " (seconds): ")

            decd = input("Please enter the DEC of reference star number " + str(i+1) + " (degrees): ")
            decm = input("Please enter the DEC of reference star number " + str(i+1) + " (minutes): ")
            decs = input("Please enter the DEC of reference star number " + str(i+1) + " (seconds): ")
            
            if i != num - 1:
                print "\n- Next Reference Star -\n"

            if rah < 0:
                ra.append( rah - ram/60.0 - ras/3600.0 )
            else:
                ra.append( rah + ram/60.0 + ras/3600.0 )
                
            if dec < 0:
                dec.append( decd - decm/60.0 - decs/3600.0 )
            else:
                dec.append( decd + decm/60.0 + decs/3600.0 )
     
        print "\n- Asteroid Inputs -\n"
        xObject = float(input("Please enter the x coordinate of the object: " ))
        yObject = float(input("Please enter the y coordinate of the ojeect: " ))

    xsum = sum(xlist)
    xsqsum = sum([x**2 for x in xlist])
    xysum = sum([x*y for x, y in zip(xlist,ylist)])
    ysum = sum(ylist)
    ysqsum = sum([y**2 for y in ylist])
    asum = sum(ra)
    axsum = sum([a*x for a, x in zip(ra,xlist)])
    aysum = sum([a*y for a, y in zip(ra,ylist)])
    dsum = sum(dec)
    dxsum = sum([d*x for d, x in zip(dec,xlist)])
    dysum = sum([d*y for d, y in zip(dec,ylist)])
 
    coeff1 = [[num, xsum, ysum],[xsum, xsqsum, xysum],[ysum,xysum,ysqsum]]
    const1 = [asum, axsum, aysum]
    coeff2 = [[num, xsum, ysum],[xsum, xsqsum, xysum],[ysum,xysum,ysqsum]]
    const2 = [dsum, dxsum, dysum]

    print const1
    print const2

    sol1 = cramer( coeff1, const1, False, False )
    sol2 = cramer( coeff2, const2, False, False )

    RAObject = sol1[0] + sol1[1]*xObject + sol1[2]*yObject
    DecObject = sol2[0] + sol2[1]*xObject + sol2[2]*yObject

    if RAObject > 0:
        RAObj_h = int(RAObject)
    else:
        RAObj_h = int(RAObject)
        RAObject = abs(RAObject)
    RAObj_m = int(RAObject*60%60)
    RAObj_s = round(RAObject*3600%3600%60,4)

    if DecObject > 0:
        DecObj_d = int(DecObject)
    else:
        DecObj_d = int(DecObject)
        DecObject = abs(DecObject)
    DecObj_m = int(DecObject*60%60)
    DecObj_s = round(DecObject*3600%3600%60,4)


    
    print "RA of Object: ", RAObj_h, RAObj_m, RAObj_s
    print "Dec of Object: ", DecObj_d, DecObj_m, DecObj_s

    RARes = []
    for i in range(num):
        RARes.append( ( ra[i] - ( sol1[0] + sol1[1]*xlist[i] + sol1[2]*ylist[i] ) ) )
        
    DecRes = []
    for i in range(num):
        DecRes.append( ( dec[i] - ( sol2[0] + sol2[1]*xlist[i] + sol2[2]*ylist[i] ) ) )

    print "\n- Star Residuals -\n"
    print "Star #  RA        Dec"
    for i in range(num):
        print "Star " + str(i) + ":%10.6f" % RARes[i] + "%10.6f" %DecRes[i]

    print "\n- Standard Deviation -\n"

    RAsq = sum([a**2 for a in RARes])
    Decsq = sum([d**2 for d in DecRes])

    RAstd = round(pow(RAsq / ( num - 3 ), 0.5), 8)
    Decstd = round(pow(Decsq / ( num - 3 ), 0.5), 8)

    print "RA Standard Deviation: ", RAstd
    print "Dec Standard Deviation: ", Decstd

    covariance = (sum([a*d for a, d in zip(RARes,DecRes)]))/num
    print "\n\nCovariance: ", covariance

    choice2 = raw_input("\nWould you like an output file (Y/N)?: ")
    
    if choice2 == "Y" or choice2 == "y" or choice2 == "Yes" or choice2 == "yes":
        filePath = raw_input("Please enter the output file name (without .txt): ")
        fileOut = open( filePath + ".txt", "w" )

        fileOut.write("LSPR Results")
        fileOut.write("\n------------")
                      
        fileOut.write("\nRA of Object: " + str(RAObj_h) + " " + str(RAObj_m) + " " + str(RAObj_s) )
        fileOut.write("\nDec of Object: "+ str(DecObj_d) + " " + str(DecObj_m) + " " + str(DecObj_s) )

        fileOut.write("\n\n- Star Residuals -\n")
        fileOut.write("\nStar #  RA        Dec ")
        for i in range(num):
            fileOut.write("\nStar " + str(i) + ":%10.6f" % RARes[i] + "%10.6f" %DecRes[i])

        fileOut.write("\n\n- Standard Deviation -\n")
        fileOut.write("\nRA Standard Deviation: " + str(RAstd) )
        fileOut.write("\nDec Standard Deviation: " + str(Decstd) )

        fileOut.write("\n\nCovariance: " + str(covariance))

        fileOut.close()

    print "Done!"

LSPR()
