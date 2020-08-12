import re
import ast
import time
import os

cwd = os.getcwd()
os.chdir(cwd + "/jacobian_arith")
load("c34.sage")
os.chdir(cwd)
load("picard.sage")

current_data_path = cwd + "/data/hwlpolys_1.txt"

#We make a finite list because these are the only 2 mod 3 Fermat primes < 2^2048.
#Computing the input to the algorithm past this point before the heat death of the universe would be quite hard
    
FermatPrimes = [5,17,257,65537]
            
#For arithmetic in Q(zeta_3).
    
K.<z>=NumberField(x^2+3)
zeta = 1/2 * (-1+z)


def findLpolys(data_path, C, lower_limit, limit):
    """
    Computes the Lpolynomials for primes lowerlimit < p <= limit, given a file of cartier manin matrices at
     data_path for the curve C.
    :param data_path: the file path to the cartier manin matrices formatted as a 3*3 list
    :param C: the curve specified C in the form [1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2]
    :param lower_limit: the lower limit on where to start computing
    :param limit: the upper limit on computing the primes
    :return: a file with Lpolynomials for primes lowerlimit < p <= limit.
    """

    cm_list = []
    limit = prime_pi(limit)-2 #convert into line numbers.

    #import urllib2  # the lib that handles the url stuff
    #file_handle = urllib2.urlopen(path)

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
                
                C = [ZZ(list_coef[0]), ZZ(list_coef[1]), 0, ZZ(list_coef[2]), 0, 0, 0, 0]
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

    for cm_coeffs in cm_list:
        p = cm_coeffs[0]
        CartierManin = cm_coeffs[1]
        if (p >= lower_limit):
            answer = [CartierManin, p, liftLpoly(CartierManin, p, C)]
            ret.append(answer)
    T2 = time.time()

    print("done in " + str(T2 - T1) + " seconds")
    print("average time " + str((T2 - T1) / len(ret)) + " seconds")
    return ret
