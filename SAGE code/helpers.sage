
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