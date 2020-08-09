import re
import ast
import time
import os

cwd = os.getcwd()
os.chdir(cwd + "/jacobian_arith")
load("c34.sage")
os.chdir(cwd)
load("picard.sage")

data_path = cwd + "/data/hwlpolys_1.txt"


def findLpolys(data_path, C, lowerlimit, limit):
    #Reading file for L polys modulo p
    #Things are stored as a list of lists of the form [p, 1, a1, a2, a3]

    cm_list = []
    limit = prime_pi(limit)-2 #convert into line numbers.

#     import urllib2  # the lib that handles the url stuff
#     file_handle = urllib2.urlopen(path)

    with open(data_path, 'r') as file_handle:

        i = 0
        for line in file_handle:
            if i==0:
                i+=1
                
                # line = "[x^4 + 92*x^2*z^2 + 84*x*z^3 - y^3*z + 10*z^4] 7 4093" for testing with z
                #Read off coefficients of the curve to store in C. The following is to format the defining polynomial.
                
                poly = re.match(r"[^[]*\[([^]]*)\]", "".join(line.split()[:-2])).groups()[0]
                poly = poly.replace(" ", "")[4:] # removing the initial "y^3-"
                
                #Convert to actual polynomial.
                
                tempRing.<x> = QQ[]
                poly = tempRing(poly) 
                list_coef = poly.coefficients(sparse=False)
                
                #The coefficients for C are in the form [1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2]
                
                C = [ZZ(list_coef[0]), ZZ(list_coef[1]), 0, ZZ(list_coef[2]), 0, 1, 0, 0]
            else:
                i += 1
                if(i>1 and i<=limit):
                    line = line[:-1].split(',', 1)
                    arr = ast.literal_eval(line[1])
                    p = int(line[0])
                    cm_list.append([p, matrix(arr)])
                elif(i>1):
                    break
    
    ret = []
    
    print "starting"
    T1 = time.time()
    
    for pcoeffs in polyL:
        p = pcoeffs[0]
        coeffs = pcoeffs[1:]
        if p >= lowerlimit and p%3==2:
            answer=[coeffs,p,liftLpoly(coeffs, p, C)]
            ret.append(answer)
    T2 = time.time()
    
    print "done in "+str(T2-T1)+" seconds"
    print "average time "+str((T2-T1)/len(ret))+" seconds"
    return ret
