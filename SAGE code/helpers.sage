
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
