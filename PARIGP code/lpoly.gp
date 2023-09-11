lift2mod3(p,cartierManin,C)=
{
    cartierManin = matrix(3,3,a,b,cartierManin[a][b]);

    \\ The following formula works because all other entries besides those used below are 0 in the 3x3 Cartier-Manin matrix.

    A2modp = lift(Mod(-cartierManin[3,1]*cartierManin[1,3]-cartierManin[3,2]*cartierManin[2,3],p));

    f = x^4 + Mod(C[3],p)*x^2 + Mod(C[2],p)*x + Mod(C[1],p);

    \\ The factorization of f determines the characteristic polynomial of Frob_p on the geometric 3-torsion. The end result gives a residue agreeing with the output 0/1 of polisirreducible.
    
    mod3res = polisirreducible(f);

    candidates = 0;
    jacvals = List();

    for(i=-1,2, A2 = A2modp + i*p; if(Mod(A2,3)==mod3res, candidates += 1; listput(jacvals,[1,0,A2,0,p*A2,0,p^3]);););

    if(candidates == 1, ret = jacvals[1]; return(concat(["1 + (",Str(ret[2]),")T + (",Str(ret[3]),")T^2 + (",Str(ret[4]),")T^3 + (",Str(ret[5]),")T^4 + (",Str(ret[6]),")T^5 + (",Str(ret[7]),")T^6"])););

    f0 = Mod(C[1],p);
    f1 = Mod(C[2],p);
    f2 = Mod(C[3],p);

    \\ Magic polynomial to compute F_p(Jac[2]) via spl(psi_f) (in this case). We only care about its roots over F_p to compute L modulo 2 as they determine the characteristic polynomial of Frob_p.

    psi_f = x^9 + 24*f2*x^7-168*f1*x^6+ (1080*f0-78*f2^2)*x^5+ 336*f1*f2*x^4+ (1728*f0*f2-636*f1^2+ 80*f2^3)*x^3+ (-864*f0*f1-168*f1*f2^2)*x^2+ (-432*f0^2 + 216*f0*f2^2-120*f1^2*f2-27*f2^4)*x-8*f1^3;

    binp = binary(p);
    
    xp_bar = 1;
    res = Mod(x,psi_f);

    \\ Use repeated squaring to efficiently compute x^p mod psi_f.
    
    for(i=1, length(binp), if(binp[length(binp)+1-i]==1, xp_bar *= res;); res = res*res;);

    \\ This computation is now O(1) because we have reduced the degree. This encodes whether psi_f has roots or not.

    psiGCDdeg = poldegree(gcd(lift(xp_bar)-x,psi_f));

    \\The parity of A2.

    parity = (psiGCDdeg > 0);

    for(i=1,2, if(Mod(jacvals[i][3],2) == parity, ret=jacvals[i]; return(concat(["1 + (",Str(ret[2]),")T + (",Str(ret[3]),")T^2 + (",Str(ret[4]),")T^3 + (",Str(ret[5]),")T^4 + (",Str(ret[6]),")T^5 + (",Str(ret[7]),")T^6"]));););
}

nf = bnfinit(z^2+3);
zzeta = nfeltdiv(nf, nfeltadd(nf, -1, z), 2);

lift1mod3(p,cartierManin,t,useTr)=
{
    if(p<53, return("Use naive point counting."););

    cartierManin = matrix(3,3,a,b,cartierManin[a][b]);

    \\ Compute zeta_3 from trace of Frobenius, or using a randomized algorithm otherwise.

    if(useTr==True, v = sqrtint((4*p-t^2)/12); sqrt_3 = Mod(t/(2*v),p); zeta3 = (Mod(-1,p)+ sqrt_3)/(Mod(2,p)););
    if(useTr==False, zeta3 = (Mod(-1,p)+sqrt(Mod(-3,p)))/Mod(2,p););

    zeta3bar = zeta3^2;

    M = [Mod(1,p), zeta3; Mod(1,p), zeta3bar];
    r_p = Mod(cartierManin[1,1]+cartierManin[2,2], p);
    r_pbar = Mod(cartierManin[3,3], p);

    \\ Solve the linear system for a_p. centerlift makes sure we pick the right lift of the residue so |a_p| is small.

    a_p = matsolve(M, Col([r_p, r_pbar]));
    a_p = [centerlift(a_p[1]), centerlift(a_p[2])];

    \\ To deduce the L-polynomial, it remains to find c_p. By Lemma 3.1, we just need to compute the corresponding 6th root of unity assuming p is ordinary.
    \\ Create 6th root of unity in Q(zeta_3) (this is "nf"). Then compute the ideal (pi), where (p)=(pi)(pibar).

    zeta6 = nfeltadd(nf, zzeta, 1);

    pi_ideal_element = nfeltadd(nf, -1*lift(zeta3), zzeta);
    pi_ideal_1 = idealhnf(nf, pi_ideal_element);
    pi_ideal_2 = idealhnf(nf, p);
    pi_ideal = idealadd(nf, pi_ideal_1, pi_ideal_2);

    \\ Compute actual generator of (pi), as well as (pibar).

    pi_00 = bnfisprincipal(nf, pi_ideal)[2];
    pibar_00 = nfeltadd(nf, pi_00[1]-pi_00[2], -1*pi_00[2]*zzeta);

    \\ Create the residue field to put into nfmodpr.

    Kp = [p, pi_00,1,1,[pibar_00[1],-pibar_00[2]; pi_00[2],-pi_00[1]]];

    sigma_zeta6 = nfmodpr(nf, zeta6, Kp);
    sigma_pibar = nfmodpr(nf, pibar_00, Kp);

    \\ Compute s as a determinant of the upper 2x2 block in the CM matrix (g_p^sigma is a reverse characteristic polynomial for the block) then pass to the residue field to get s_p.

    s = cartierManin[1,1] *  cartierManin[2,2] - cartierManin[1,2]*cartierManin[2,1];
    s_p = nfmodpr(nf, s, Kp);

    \\Use the image sextic_val = sigma(zeta6) to determine the power sextic_power=i of zeta6 so c_p = zeta6^i * pibar * p.

    r_pbar_ffelt = nfmodpr(nf, cartierManin[3,3], Kp);
    sextic_val = s_p / (r_pbar_ffelt * sigma_pibar );

    my(sextic_power = 0);
    sextic_power = 0;

    for(i = 0, 5, if(sigma_zeta6^i == sextic_val, sextic_power = i));

    \\ Find a_p, c_p as field elements in Q(zeta3)=nf. We may compute b_p from these, and thus both of L_p^sigma, L_p^\sigmabar.
    
    if(a_p[2]!=0,
       a_pbar = nfeltadd(nf, lift(a_p[1]-a_p[2]), -zzeta*lift(a_p[2]));
       a_p = nfeltadd(nf, lift(a_p[1]), lift(a_p[2])*zzeta);
    );
    if(a_p[2]==0,
       a_pbar = Col([a_p[1],0]);
       a_p = a_pbar;
    );
    c_p = nfeltmul(nf, nfeltmul(nf, nfeltpow(nf, zeta6, sextic_power), pibar_00), p);
    c_pbar = nfeltadd(nf, c_p[1]-c_p[2], -1*c_p[2]*zzeta);

    b_p = nfeltdiv(nf, nfeltmul(nf, c_p, a_pbar),p);
    b_pbar = nfeltadd(nf, b_p[1]-b_p[2], -1*b_p[2]*zzeta);

    \\ Compute the L polynomial explicitly.
    
    A1 = -2*(a_p[1])+a_p[2];
    A2 = nfeltadd(nf, nfeltmul(nf, a_p, a_pbar), nfeltadd(nf, b_p, b_pbar));
    A3 = -nfeltadd(nf, nfeltadd(nf, c_p, c_pbar), nfeltadd(nf, nfeltmul(nf, a_p, b_pbar), nfeltmul(nf, a_pbar, b_p)));
    return(concat(["1 + (",Str(A1),")T + (",Str(A2),")T^2 + (",Str(A3),")T^3 + (",Str(p*A2),")T^4 + (",Str(p^2*A1),")T^5 + (",Str(p^3),")T^6"]));
}

\\ To use this function: pre-compute Cartier-Manin matrices for C, and compute a list of Frobenius traces of some elliptic curve
\\ E/Q with CM by Q(sqrt(-3)) and CM discriminant -12, such as y^2 = x^3 - 15x + 22.
\\ The Frobenius traces begin at the prime 2; depending on primes skipped in CM matrices, you'll want to shift the index.
\\ You can remove TrList and compute zeta_3 directly from built-in PARI/GP functions via liftLpolysAlt if you don't want to deal with this, the time is not too much slower.
\\ The Cartier-Manin matrices are put in CMList in the format [p,array for CM matrix].
\\ We have provided sample files to use for CMList, TrList in the repository.
\\ To use these with this code, read in CMList via readstr("data_file") - it is processed correctly in the function.
\\ TrList can be directly inputted as a list.
\\ The list C = [C[1],C[2],C[3]] corresponds to the Picard curve y^3 = x^4+C[3]x^2+C[2]x+C[1].

liftLpolys(CMList,TrList,C,i)=
{
    \\Process string into usable form.

    string = CMList[i+1];
    l = length(string);
    procstring = Vec(string);
    index = 0;
    for(i=1, l, if(procstring[i]==",", index=i; i=l+1;););

    p = "";
    for(i=1, index-1, p=concat(p,procstring[i]););
    p = eval(p);

    CM = "";
    for(i=index+1, l, CM=concat(CM,procstring[i]););
    CM = eval(CM);

    t = TrList[i];
    if(Mod(p,3)==1, return(lift1mod3(p,CM,t,True)););
    if(Mod(p,3)==2, return(lift2mod3(p,CM,C)););
}

liftLpolysAlt(CMList,C,i)=
{
    \\Process string into usable form.

    string = CMList[i+1];
    l = length(string);
    procstring = Vec(string);
    index = 0;
    for(i=1, l, if(procstring[i]==",", index=i; i=l+1;););

    p = "";
    for(i=1, index-1, p=concat(p,procstring[i]););
    p = eval(p);

    CM = "";
    for(i=index+1, l, CM=concat(CM,procstring[i]););
    CM = eval(CM);

    if(Mod(p,3)==1, return(lift1mod3(p,CM,1,False)););
    if(Mod(p,3)==2, return(lift2mod3(p,CM,C)););
}