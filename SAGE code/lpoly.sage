import re
import ast
import time
from math import sqrt

#The c34.sage file is from https://github.com/emmacneil/c34-curves, and implements Jacobian arithmetic efficiently.

load("c34.sage")
            
#For arithmetic in Q(zeta_3).

K.<z>=NumberField(x^2+3)
zeta = 1/2 * (-1+z)

def findLpolys(cm_path, lower_limit, upper_limit, skipJac):
    """
    Computes the Lpolynomials for primes lowerlimit < p <= limit, given a file of cartier manin matrices at
     cm_path for the curve C.
    :param cm_path: the file path to the cartier manin matrices formatted as a 3*3 list
    :param lower_limit: the lower limit on computing the primes
    :param upper_limit: the upper limit on computing the primes
    :param skipJac: determines whether or not to use the Jacobian when p is 2 mod 3. Skipping (skipJac is True) tends to be slightly faster
    :return: a file with Lpolynomials for primes lowerlimit < p <= limit.
    """
    
    cm_list = []
    line_limit = prime_pi(upper_limit)-2 #convert into line numbers.

    #import urllib2  # the lib that handles the url stuff
    #file_handle = urllib2.urlopen(path)

    with open(cm_path, 'r') as file_handle:

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
                
                C = [ZZ(list_coef[0]), ZZ(list_coef[1]), 0, ZZ(list_coef[2]), 0, 0, 0, 0, 0]
            else:
                i += 1
                if(i>1 and i<=line_limit):
                    line = line[:-1].split(',', 1)
                    arr = ast.literal_eval(line[1])
                    p = int(line[0])
                    cm_list.append([p, matrix(arr)])
                elif(i>1):
                    break
    
    ret = []
    
    print("Starting")

    T1 = time.time()

    for i in range(line_limit-1):
        cm_coeffs=cm_list[i]
        p = cm_coeffs[0]
        CartierManin = cm_coeffs[1]
        if(p>lower_limit):
            poly=liftLpoly(CartierManin, p, C, skipJac)
            ret.append([p,poly])
            
    T2 = time.time()

    print("Done in " + str(T2 - T1) + " seconds")
    print("Average time: " + str((T2 - T1) / len(ret)) + " seconds")
    return ret


#gcds in Z[zeta_3]
CC=ComplexField(30)

#Set up constants here, so that they are not recomputed.
const = RR(1/sqrt(3))
I = CC.0
z3 = CC(-1+I*sqrt(3))/2

def subtract(e1,e2):
    """
    Given e1, e2 in Z[zeta_3], computes e1-e2
    :param e1: an element of Z[zeta_3]
    :param e2: an element of Z[zeta_3]
    :return: e1 - e2 as a pair [x, y] representing x + y*zeta_3
    """
    return [e1[0]-e2[0], e1[1]-e2[1]]

def mult(e1,e2):
    """
    Given e1, e2 in Z[zeta_3], computes e1*e2
    :param e1: an element of Z[zeta_3]
    :param e2: an element of Z[zeta_3]
    :return: e1 * e2 as a pair [x, y] representing x + y*zeta_3
    """
    return [e1[0]*e2[0]-e1[1]*e2[1], e1[1]*e2[0]+e1[0]*e2[1]-e1[1]*e2[1]]
    
def find_rem(e1,e2):
    """
    Given e1, e2 in Z[zeta_3], computes e1 mod e2
    :param e1: an element of Z[zeta_3]
    :param e2: an element of Z[zeta_3]
    :return: e1 mod e2 as a pair [x, y] representing x + y*zeta_3
    """
    #Computes the remainder of two elements of Z[zeta_3].
    
    e1complex = e1[0]+z3*e1[1]
    e2complex = e2[0]+z3*e2[1]
    
    quot = e1complex/e2complex
    
    a = quot.real()
    b = quot.imag()
    
    kappa = [(round(a+const*b)),(round(2*const*b))]
    remainder = subtract(e1, mult(kappa,e2))
    return remainder

def Euclidean_GCD(e1,e2):
    """
    Given e1, e2 in Z[zeta_3], computes gcd(e1, e2)
    :param e1: an element of Z[zeta_3]
    :param e2: an element of Z[zeta_3]
    :return: gcd(e1, e2) as a pair [x, y] representing x + y*zeta_3
    """
    #Store e_i = [a,b] for a+b*zeta_3
    #Recursively use Euclidean algorithm.
    
    if(e2==[0,0]):
        return e1
    else:
        remainder = find_rem(e1,e2)
        return Euclidean_GCD(e2, remainder)

def spanorder(P1,P2,p,E,C):
    #g,h are elements we are finding the span of
    #Order of G is p^E.
    
    e1 = findorderexp(P1,E,p,C)
    e2 = findorderexp(P2,E,p,C)
    
    print(e1,e2)
    
    if(e1==0 or e2==0): return [P1,P2,e1,e2]
    
    if(e1<=e2): 
        (x,y) = DL(P2,P1,e2,p,C)
        if(-x/y>=0): Q1 = P1+P2*(-x/y)
        elif(-x/y<0): Q1 = P1+P2*(int(p^E-x/y))
        Q2 = P2
        n1 = findorderexp(Q1,E,p,C)
        n2 = findorderexp(Q2,E,p,C)
        return([Q1,Q2,n1,n2])
    else: 
        (x,y) = DL(P1,P2,e1,p,C)
        if(-x/y>=0): Q1 = P2+P1*(-x/y)
        elif(-x/y<0): Q1 = P2+P1*(int(p^E-x/y))
        Q2 = P1
        n1 = findorderexp(Q1,E,p,C)
        n2 = findorderexp(Q2,E,p,C)
        return([Q1,Q2,n1,n2])
    
def DL(P1,P2,e,p,C):
    exp=1
    while(True):
        x = pohlighellman(P1,P2,e,p,C)
        if(x != None): return(x,exp)
        exp *= p
        P2=P2*p

def findorderexp(D,E,p,C):
    # Order of G is p^E.
    for i in range(E+1):
        if(p^(i)*D == C.zero_divisor()): return i

def pohlighellman(g,h,e,p,C):
    x_k = 0
    gamma = g*(p^(e-1))
    
    for k in range(e):
        h_k = (g*(-x_k)+h)*(p^(e-1-k))
        #Find d_k so gamma^d_k = h_k using BSGS
        d_k = logp(gamma,h_k,p,C)
        
        if(d_k == None): return None
        
        x_k = x_k+p^k * d_k
        
    return x_k

def logp(alpha,beta,p,C):
    #Find x so alpha^x=beta
    #We use this for p=3 so this is trivial
    #In general, should use BSGS
    
    if(beta == C.zero_divisor()): return 0
    if(alpha==beta): return 1
    if(2*alpha==beta): return 2

    # Solution not found
    return None

def threesylowislarge(N3free,E,C,bound):
    
    P1 = C.random_divisor()
    P2 = C.random_divisor()
    
    P1 = P1*(N3free)
    P2 = P2*(N3free)
    
    L = spanorder(P1,P1,3,E,C)
    
    if(L[2]+L[3]>bound): return True
    
    Q1 = L[0]
    Q2 = L[1]
    n1 = L[2]
    n2 = L[3]

    P3 = C.random_divisor()
    P3 = P3*(N3free)

    InSpan = False

    for i in range(3^n1):
        for j in range(3^n2):
            if(P3 == Q1*i+Q2*j): InSpan = True

    if(n1+n2 >= bound and InSpan == False): return True
    
    return False

def liftLpoly(CartierManin, p, C, skipJac):
    """
    lifts a cartier manin matrix for the fixed prime p to the almost surely unique Lpolynomial.
    :param CartierManin: formatted as a matrix [[*,*,*], [*,*,*], [*,*,*]]
    :param p: the prime p
    :param C: the curve specified C in the form [1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2]
    :param skipJac: boolean for whether or not the Jacobian is used or not in the 2 mod 3 case
    :return: the lifted L polynomials Lp with its coefficients [1, x, x^2, x^3, x^4, x^5, x^6]
    """
        
    if(p%3==1):
        #We are now in the split case. We require the prime to be ordinary, so check this first.
        
        if(CartierManin.rank()<3): return "Non-ordinary prime."
        
        #First step is to find a cube root mod p and factor p in Z[zeta_3].
        
        #Finding cube root (using traces on elliptic curve with CM by Q(sqrt(-3)) and disc = -12):
        
        k = GF(p)
        
        sqrt_3 = (k(-3)).sqrt()
        zeta3 = k(-1+sqrt_3)/(k(2))

        zeta3bar = zeta3^2
        
        #Now we solve the linear system (3.3) in the paper.
        
        M = matrix([[k(1),zeta3],[k(1),zeta3bar]])
        Minv = M.inverse()
        
        #Find r_p, r_pbar from CartierManin. These are coefficients in the characteristic polynomials of the blocks
        #which correspond to the trace in each case, which is why they are just adding diagonal entries.
        
        r_pbar = k(CartierManin[2][2])
        r_p = k(CartierManin[0][0]+CartierManin[1][1])
        
        #This vector has the x and y coordinates of a_p (Z, Z*zeta_3 parts in Z[zeta_3])
        
        a_p = Minv * vector([r_p,r_pbar])
        
        #Determine a_p from its reduction modulo p.
        
        a_p = [a_p[0],a_p[1]]
        
        for i in range(2):
            entry = a_p[i]
            if(ZZ(entry)>p/2):
                a_p[i]=ZZ(entry)-p
                
        #Determine pi via the Extended Euclidean algorithm in Z[zeta_3].
        #The objects K and zeta are defined once outside of this function because they are expensive to construct.
        
        pi_ideal = K.ideal(-ZZ(zeta3)+zeta,p)
        
        #Here, we use [a,b] to denote a+b*zeta (zeta is the cube root of unity in C)
        
        pi = Euclidean_GCD([-ZZ(zeta3),1],[p,0])
        pi = (pi[0])+(pi[1])*zeta

        
        #(Optional sage implementation of factoring - speed is *very* slightly slower, small difference)
        #pi = pi_ideal.gens_reduced()[0]
        
        pibar = pi[0]-z*pi[1]
        
        #Sixth root of unity.
        
        zeta6 = zeta+1
        
        #Apply the map sigma to get elements of F_p (the residue field)
        
        Kp = pi_ideal.residue_field()
        
        sigma_zeta6 = (Kp(zeta6))
        sigma_pibar = (Kp(pibar))
        
        #Determine the sextic root of unity for finding the value of c_p. 
        #This amounts to finding the correct exponent.
        
        s_p = k(CartierManin[0][0]*CartierManin[1][1]-CartierManin[0][1]*CartierManin[1][0])
        sextic_val = s_p/(r_pbar*sigma_pibar)
        
        sextic_power = 0
        
        #Determine the correct power of zeta6.
        
        for i in range(6):
            if(sigma_zeta6^i == sextic_val): 
                sextic_power = i
        
        # Now we have a_p, b_p, c_p. We may return L_p(C,T).
        
        a_p = ZZ(a_p[0])+ZZ(a_p[1])*zeta
        a_pbar = a_p[0]-z*a_p[1]
        c_p = zeta6^(sextic_power)*pibar*p
        c_pbar = c_p[0]-z*c_p[1]
        b_p = c_p * a_pbar/p
        b_pbar = b_p[0]-z*b_p[1]
        
        #Coefficients of L-polynomial.
        
        A1 = -(a_p+a_pbar)
        A2 = (a_p*a_pbar)+(b_p+b_pbar)
        A3 = -(c_p+c_pbar)-(a_p*b_pbar + a_pbar*b_p)
        return [1,A1, A2, A3, p*A2, p^2*A1, p^3]

    elif (p%3==2 and p>= 877 and skipJac==False):
        
        #Start by computing coefficients A_1, A_2, A_3 modulo p for the L-polynomial.

        A_1 = 0
        A_3 = 0

        #A_1,A_3 are now correct. We read off A_2 from the Cartier-Manin matrix.
        
        #L_p(t) = (pt^2+1)(p^2 t^4+At^2-pt^2+1). So A_1=0, A_3=0. A_2 = A (mod p).
        
        k = GF(p)
        A_2modp = ZZ(k(-CartierManin[2][0]*CartierManin[0][2]-CartierManin[2][1]*CartierManin[1][2]))
        
        #List to store candidate values for the order of the Jacobian.
        
        jacvals=[]
        
        curve = C34Curve(k, C)

        #We construct this polynomial to check irreducibility and use Lemma 3.9.
        
        PolyRing.<t> = PolynomialRing(k)
        
        f = t^4 + k(C[3])*t^2 + k(C[1])*t + k(C[0])
        
        mod3res = 0
        
        if((f).is_irreducible()): mod3res = 1
        
        #We have now computed the residue of the central coefficient of the quartic factor modulo 3.
        #Now, search for possible candidates within the bounds of Lemma 3.8.
        #The bounds say A_2 is between -p and 3p, as the central coefficient has |A|<=2p.

        for i in range(-2,3):
            A_2 = A_2modp + i*p
            if(A_2%3==mod3res):
                JC = 1+A_2+A_2*p+p^3
                jacvals.append([JC, [1,0,A_2,0,p*A_2,0,p^3]])
       
        numruns = 0
        
        while(len(jacvals)>1):
            
            D = curve.random_divisor()
            numruns+=1

            if(numruns >= 5 and (p in FermatPrimes)): 
                return "Use O(p^1/4) algorithm instead."
                #This is effectively finitely many cases (which don't occur in practice), so not worth implementing.
            
            if(numruns >= 5 and jacvals[1][0]==(p+1)^3):
                
                #We are in the exceptional case where s_p=p. The possibilities are either Jac = (p+1)^3, or (p+1)(p^2-p+1).
                #The following tests for (p+1)^3 using the 3-sylow.
                
                if(numrums==5):
                    N = p+1
                    E = 0
                    while(N%3==0): 
                        N=int(N/3)
                        E+=1

                    bound = E+1

                    #Get 3-free part N + exponent E of 3-sylow for the candidate with Jacobian order (p+1)^3
                    E = 3*E
                    N = N^3

                #If #Jac = (p+1)^3, then with O(1) attempts we will prove the 3-sylow is large enough
                #Can eliminate the other candidate (whose 3-sylow is of order 3^(bound))
                
                if(threesylowislarge(N,E,C,bound)): return jacvals[1][1]
                
                #If this method cannot find a large enough 3-sylow, it is found by the rest of the loop: 
                #There exists a prime dividing the group exponent lambda but not p+1.
            
            Jac = jacvals[0]
                
            N_1D = Jac[0]*D

            if N_1D != curve.zero_divisor():
                jacvals.remove(Jac)
            
            else:
                #Compute D2 so that we use less additions.
                D2 = (3*(p+1)*p)*D
                
                if(N_1D + D2 != curve.zero_divisor()):
                    jacvals.remove(jacvals[1])

        #Return the only candidate L-polynomial left.
        return jacvals[0][1]
    
    elif(p%3==2 and skipJac):
        
        k=GF(p)

        #The following formula works because all other entries besides those used below are 0 in the 3x3 Cartier-Manin matrix.

        A2modp = ZZ((-CartierManin[2,0]*CartierManin[0,2]-CartierManin[2,1]*CartierManin[1,2])%p)
        
        R.<x> = PolynomialRing(k, 'x')

        f = x^4 + k(C[3])*x^2 + k(C[1])*x + k(C[0])

        #The factorization of f determines the characteristic polynomial of Frob_p on the geometric 3-torsion. The end result gives a residue agreeing with the output 0/1 of polisirreducible.

        mod3res = 0
        
        if((f).is_irreducible()): mod3res = 1

        candidates = 0;
        jacvals = []

        for i in range(-1,3):
            A2 = A2modp + i*p
            if(Mod(A2,3)==mod3res):
                candidates += 1
                jacvals.append([1,0,A2,0,p*A2,0,p^3])

        if(candidates == 1):
            return jacvals[0]

        f0 = k(C[0])
        f1 = k(C[1])
        f2 = k(C[3])

        # Magic polynomial to compute F_p(Jac[2]) via spl(psi_f) (in this case). We only care about its roots over F_p to compute L modulo 2 as they determine the characteristic polynomial of Frob_p.

        psi_f = R(x^9 + 24*f2*x^7-168*f1*x^6+ (1080*f0-78*f2^2)*x^5+ 336*f1*f2*x^4+ (1728*f0*f2-636*f1^2+ 80*f2^3)*x^3+ (-864*f0*f1-168*f1*f2^2)*x^2+ (-432*f0^2 + 216*f0*f2^2-120*f1^2*f2-27*f2^4)*x-8*f1^3)

        binp = bin(p)
        
        S.<y>=R.quotient(psi_f)

        xp_bar = 1
        res = S(y)

        #Use repeated squaring to efficiently compute x^p mod psi_f.

        for i in range(len(binp)-2):
            if(binp[len(binp)-i-1]=='1'):
                xp_bar = xp_bar*res
            res = res*res

        #This computation is now O(1) because we have reduced the degree. This encodes whether psi_f has roots or not.

        psiGCDdeg = ((lift(xp_bar)-x).gcd(psi_f)).degree()

        #The parity of A2.

        parity = (psiGCDdeg > 0)

        for i in range(2):
            if(jacvals[i][2]%2 == +(parity)): 
                return(jacvals[i])