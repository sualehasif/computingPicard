
#gcds in Z[zeta_3] (as it turns out, SAGE's method is faster. Maybe this can be optimized)
CC=ComplexField(30)

#Set up constants here, so that they are not recomputed.
const = RR(1/sqrt(3))
I = CC.0
z3 = CC(-1+I*sqrt(3))/2

def subtract(e1,e2):
    return [e1[0]-e2[0], e1[1]-e2[1]]

def mult(e1,e2):
    return [e1[0]*e2[0]-e1[1]*e2[1], e1[1]*e2[0]+e1[0]*e2[1]-e1[1]*e2[1]]
    
def find_rem(e1,e2):
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
    #Store e_i = [a,b] for a+b*zeta_3
    #Recursively use Euclidean algorithm.
    
    if(e2==[0,0]):
        return e1
    else:
        remainder = find_rem(e1,e2)
        return Euclidean_GCD(e2, remainder)

def Dpowers(l,D):
    Dpows = [D]
    
    for i in range(l-1):
        Dtemp = Dpows[len(Dpows)-1]
        D = 2*D
        Dpows.append(D)
    #Dpows now has all needed powers of 2^k D.
    return Dpows

def multbyN(N,D,Dpows):
    #f0 = C[0]
    #f1 = C[1]
    #f2 = C[2]
    #h0 = C[3]
    #h1 = C[4]
    #h2 = C[5]
    #h3 = C[6]
    #Fast conversion to binary.
    binaryN = bin(N)[2:]
    #Store log_2 N
    l = len(binaryN)
    #Store powers.
    Dtemp = []
    
    for i in range(l):
        D = Dpows[l-1-i]
        if(binaryN[i]=='1'):
            if(Dtemp == []): Dtemp = D
            else: Dtemp = D+Dtemp
                
    return Dtemp

