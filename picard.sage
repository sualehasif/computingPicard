from math import sqrt
load("helpers.sage")

def liftLpoly(coeffs, p, C):
    if p<=151:
        return Lpoly(coeffs,p,C)
        #For small cases. Can probably be dramatically improved; was somewhat ignored for this project.
        
    elif(p%3==1 and p>11):
        #We are now in the split case. We require the prime to be ordinary, so check this first.
        
        if(CartierManin.rank()<3): return "Non-ordinary prime."
        
        #First step is to find a cube root mod p and factor p in Z[zeta_3].
        
        #Finding cube root (using Tonelli-Shanks for speed):
        
        k = GF(p)
        zeta3 = k(-1+square_root_mod_prime(mod(-3,p)))/(k(2))
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
        return [1,A1, A2, A3, p*A1, p^2*A2, p^3*A3]
    elif (p%3==2 and p>= 877):
        
        
        A_1 = 0
        A_3 = 0
        #A_1,A_3 are now correct
        
        #L_p(t) = (pt^2+1)(p^2 t^4+(a/2)t^2-pt^2+1). So A_1=0, A_3=0. A_2=a/2.
        
        k = GF(p)
        A_2modp = ZZ(k(-CartierManin[2][0]*CartierManin[0][2]-CartierManin[2][1]*CartierManin[1][2]))
        
        jacvals=[]
        
        curve = C34Curve(k, C)
        
        PolyRing.<t> = PolynomialRing(k)
        
        f = t^4 + k(C[3])*t^2 + k(C[1])*t + k(C[0])
        
        mod3res = 0
        
        if((f).is_irreducible()): mod3res = 1
        
        for i in range(-2,3):
            A_2 = A_2modp + i*p
            if(A_2%3==mod3res):
                JC = 1+A_2+A_2*p+p^3
                jacvals.append([JC, [1,0,A_2,0,p*A_2,0,p^3]])
                    
        l = len(bin(p+1))*3+1
        numruns = 0
        
        while(len(jacvals)>1):
            
            D = curve.random_divisor()
            Dpows = Dpowers(l,D,C)
            numruns+=1
            
            #Use Monte-Carlo algorithm with M=10.
            
            if(numruns > 10 and (p in FermatPrimes)): 
                return "Use O(p^1/4) algorithm instead."
            
            if(numruns > 10 and jacvals[1][0]==(p+1)^3):
                jacvals.remove(jacvals[0])
            
            Jac = jacvals[0]
                
            N_1D = multbyN(Jac[0]-1,D,Dpows)

            if N_1D != 0:
                jacvals.remove(Jac)
            
            else:
                D2 = multbyN(3*(p+1)*p, D, Dpows)
                
                if(N_1D + D2 != 0):
                    jacvals.remove(jacvals[1])

        return jacvals[0][1]
