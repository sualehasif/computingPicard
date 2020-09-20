
mulconst(a, c)=
{
    \\       \\ returns c*a, where a is supposed to be a field element and c a small
    \\       \\ constant. Put into a separate so that each occurrence of
    \\       \\ "*" in the formulae corresponds to a multiplication of two field
    \\       \\ elements. It should normally be computed by repeated additions to
    \\       \\ simplify, we use a multiplication here.
   return(c*a);
}

divconst(a, c)=
{
    \\       \\ returns a/c, where a is supposed to be a field element and c a small
    \\       \\ constant. It should normally be computed in a special way, see our
    \\       \\ article to simplify, we use a normal division here.
   return(a/c);
}


invert(a)=
{
if(a!=0, 
  return(1/a), 
  return(0));
}


    \\     \\ ^^^*********************************************************************



interpolate(X0, X1, Xm1, Xinf)=
{
    \\      \\ Cost:    for free
    \\      \\ interpolation on 0, 1, -1  and infinity
    \\  
    \\      \\ Interpolate

    W3 = Xinf;
    W0 = X0;
    W2 = divconst(X1 + Xm1, 2) - W0;
    W1 = X1 - ( W0 + W2 + W3 );
    
    return([W0,W1,W2,W3]);
}

    \\     \\ ***************************************************************************



toom_cook_1_1(X0,X1,Y0,Y1)=
{
    \\      \\ Let W = X * Y first multiply without the leading coefficient by
    \\      \\ interpolation on 0, 1 and infinity

    Xval1  = X1 + X0;
    Yval1  = Y1 + Y0;

    \\      \\ Evaluate

    Wval0   = X0 * Y0;
    Wval1   = Xval1 * Yval1;
    Wvalinf = X1*Y1;

    \\      \\ Interpolate

    W2 = Wvalinf;
    W0 = Wval0;
    W1 = Wval1 - ( W0 + W2);

    return([W0,W1,W2]);
}

    \\     \\ ***************************************************************************


toom_cook_1_2(X0,X1,Y0,Y1,Y2)=
{
    \\      \\ Let W = X * Y first multiply without the leading coefficient by
    \\      \\ Cost: 4M
    \\      \\ interpolation on 0, 1, -1  and infinity

    Xval1  = X1 + X0;
    Yval1  = Y2 + Y1 + Y0;
    Xvalm1 = - X1 + X0;
    Yvalm1 = Y2 - Y1 + Y0;

    \\     \\ Evaluate

    Wval0   = X0 * Y0;
    Wval1   = Xval1 * Yval1;
    Wvalm1  = Xvalm1 * Yvalm1;
    Wvalinf = X1*Y2;

    \\     \\ Interpolate

    W3 = Wvalinf;
    W0 = Wval0;
    W2 = divconst(Wval1+Wvalm1,2) - W0;
    W1 = Wval1 - ( W0 + W2 + W3 );

    return([W0,W1,W2,W3]);
}


    \\     \\ ***************************************************************************


toom_cook_1_3(X0,X1,Y0,Y1,Y2,Y3)=
{
    \\     \\ Let W = X * Y first multiply without the leading coefficient by
    \\     \\ Cost: 5M
    \\     \\ interpolation on 0, 1, -1, 2  and infinity

    Xval1  = X1 + X0;
    Yval1  = Y3 + Y2 + Y1 + Y0;
    Xvalm1 = - X1 + X0;
    Yvalm1 = - Y3 + Y2 - Y1 + Y0;
    Xval2  = mulconst(X1,2) + X0;
    Yval2  = mulconst (Y3, 8) + mulconst (Y2, 4) + mulconst (Y1, 2) + Y0;
    \\     \\ Yval2  = mulconst ( mulconst ( mulconst (Y3, 2) + Y2, 2) + Y1,2) + Y0

    \\     \\ Evaluate

    Wval0   = X0 * Y0;
    Wval1   = Xval1 * Yval1;
    Wvalm1  = Xvalm1 * Yvalm1;
    Wval2   = Xval2 * Yval2;
    Wvalinf = X1*Y3;

    \\     \\ Interpolate

    W4 = Wvalinf;
    W0 = Wval0;
    W2 = -(Wval0 + Wvalinf) + divconst(Wval1 + Wvalm1,2);
    W3 = - mulconst(Wvalinf,2) + divconst(Wval0 - Wval1 + divconst(Wval2 - Wvalm1,3),2);
    W1 = Wval1 - ( W0 + W2 + W3 + W4);

    return([W0,W1,W2,W3,W4]);
}

    \\     \\ ***************************************************************************


toom_cook_2_2(X0,X1,X2,Y0,Y1,Y2)=
{
    \\      \\ Let W = X * Y first multiply without the leading coefficient by
    \\      \\ Cost: 5M
    \\      \\ interpolation on 0, 1, -1, 2 and infinity

    Xval1  = X2 + X1 + X0;
    Yval1  = Y2 + Y1 + Y0;
    Xvalm1 = X2 - X1 + X0;
    Yvalm1 = Y2 - Y1 + Y0;
    Xval2  = mulconst (X2, 4) + mulconst (X1, 2) + X0;
    Yval2 =  mulconst (Y2, 4) + mulconst (Y1, 2) + Y0;
    \\      \\ Xval2  = mulconst (mulconst (X2, 2) + X1, 2) + X0
\\  \\ Yval2  = mulconst (mulconst (Y2, 2) + Y1, 2) + Y0

\\ \\ Evaluate

    Wval0   = X0 * Y0;
    Wval1   = Xval1 * Yval1;
    Wvalm1  = Xvalm1 * Yvalm1;
    Wval2   = Xval2 * Yval2;
    Wvalinf = X2*Y2;

\\ \\ Interpolate

    W4 = Wvalinf;
    W0 = Wval0;
    W2 = divconst (Wval1 + Wvalm1, 2) - (W0 + W4);
    W3 = divconst (divconst (Wval2 - Wval1 + Wvalm1 - W0, 2) - mulconst (mulconst (W4, 4) + W2, 2), 3);
    W1 = Wval1 - ( W0 + W2 + W3 + W4 );

    return([W0,W1,W2,W3,W4]);
}
\\ \\ ***************************************************************************


toom_cook_square_2(X0,X1,X2)={
\\ Let W = X^2  first multiply without the leading coefficient by
\\   5SQ

\\ interpolation on 0, 1, -1, 2 and infinity

    Xval1  = X2 + X1 + X0;
    Xvalm1 = X2 - X1 + X0;
    Xval2  = mulconst (X2, 4) + mulconst (X1, 2) + X0;
\\ \\ Xval2  = mulconst (mulconst (X2, 2) + X1, 2) + X0

\\ Evaluate

    Wval0   = X0^2;
    Wval1   = Xval1^2;
    Wvalm1  = Xvalm1^2;
    Wval2   = Xval2^2;
    Wvalinf = X2^2;

\\ Interpolate

    W4 = Wvalinf;
    W0 = Wval0;
    W2 = divconst (Wval1 + Wvalm1, 2) - (W0 + W4);
    W3 = divconst (divconst (Wval2 - Wval1 + Wvalm1 - W0, 2) - mulconst (mulconst (W4, 4) + W2, 2), 3);
    W1 = Wval1 - ( W0 + W2 + W3 + W4 );

    return([W0,W1,W2,W3,W4]);
}

toom_cook_2_3 (X0,X1,X2,Y0,Y1,Y2,Y3)={
\\ Let W = X * Y first multiply without the leading coefficient by
\\ interpolation on 0, 1, -1, 2 , -2 and infinity

    Xval1  =  X2 + X1 + X0;
    Yval1  =  Y3 + Y2 + Y1 + Y0;
    Xvalm1 =  X2 - X1 + X0;
    Yvalm1 = -Y3 + Y2 - Y1 + Y0;
    Xval2  =  mulconst (X2, 4) + mulconst (X1, 2) + X0;
    Yval2  =  mulconst (Y3, 8) + mulconst (Y2, 4) + mulconst (Y1, 2) + Y0;
    Xvalm2 =  mulconst (X2, 4) - mulconst (X1, 2) + X0;
    Yvalm2 = -mulconst (Y3, 8) + mulconst (Y2, 4) - mulconst (Y1, 2) + Y0;

\\ Evaluate

    Wval0   = X0*Y0;
    Wval1   = Xval1*Yval1;
    Wvalm1  = Xvalm1*Yvalm1;
    Wval2   = Xval2*Yval2;
    Wvalm2  = Xvalm2*Yvalm2;
    Wvalinf = X2*Y3;

\\ interpolate

    W0 = Wval0;
    W5 = Wvalinf;
    W3 = -mulconst(Wvalinf,5) + divconst(divconst(Wvalm1 - Wval1 + divconst(Wval2 - Wvalm2,2),2),3);
    W4 = divconst(divconst(-(Wval1 + Wvalm1)+ divconst( mulconst(Wval0, 3) + divconst(Wval2 + Wvalm2,2) ,2) ,2),3);
    W2 = divconst( Wval1 + Wvalm1, 2) -( W0 + W4 );
    W1 = Wval1 - ( W0 + W3 + W2 + W4 + W5 );

    return([W0,W1,W2,W3,W4,W5]);
}

\\ ***************************************************************************

toom_cook_3_3 (X0,X1,X2,X3,Y0,Y1,Y2,Y3)={
\\ Let W = X * Y first multiply without the leading coefficient by
\\ interpolation on 0, 1, -1, 2 , -2, 3 and infinity
\\ 7M

    Xval1  =  X3 + X2 + X1 + X0;
    Yval1  =  Y3 + Y2 + Y1 + Y0;
    Xvalm1 = -X3 + X2 - X1 + X0;
    Yvalm1 = -Y3 + Y2 - Y1 + Y0;
    Xval2  =  mulconst (X3, 8) + mulconst (X2, 4) + mulconst (X1, 2) + X0;
    Yval2  =  mulconst (Y3, 8) + mulconst (Y2, 4) + mulconst (Y1, 2) + Y0;
    Xvalm2 = -mulconst (X3, 8) + mulconst (X2, 4) - mulconst (X1, 2) + X0;
    Yvalm2 = -mulconst (Y3, 8) + mulconst (Y2, 4) - mulconst (Y1, 2) + Y0;
    Xval3  =  mulconst (X3, 27)+ mulconst (X2, 9) + mulconst (X1, 3) + X0;
    Yval3  =  mulconst (Y3, 27)+ mulconst (Y2, 9) + mulconst (Y1, 3) + Y0;

\\ Evaluate

    Wval0   = X0*Y0;
    Wval1   = Xval1*Yval1;
    Wvalm1  = Xvalm1*Yvalm1;
    Wval2   = Xval2*Yval2;
    Wvalm2  = Xvalm2*Yvalm2;
    Wval3   = Xval3*Yval3;
    Wvalinf = X3*Y3;

\\ interpolate

    W0 = Wval0;
    W6 = Wvalinf;
    W4 = -mulconst(Wvalinf,5) + divconst(divconst(-(Wval1+Wvalm1) + divconst(mulconst(Wval0,3) + divconst(Wvalm2+Wval2,2),2),2),3);
    W3 = mulconst(mulconst(Wvalinf,3),5) + divconst(divconst(divconst( mulconst(Wval0,5)  - mulconst(Wval1,7) +divconst(- Wval3 + mulconst(Wval2,7) - Wvalm2 - Wvalm1,2),2),2),3);
    W5 = -mulconst(Wvalinf,3) + divconst(divconst(divconst(Wval1-Wval0 + divconst(Wvalm1-Wval2+ divconst(Wval3-Wvalm2,5),2),2),2),3);
    W2 = divconst( Wval1 + Wvalm1, 2) -( W0 + W4 + W6);
    W1 = Wval1 - ( W0 + W2 + W3 + W4 + W5 + W6 );

    return([W0,W1,W2,W3,W4,W5,W6]);

}

\\ ***************************************************************************


\\ ***************************************************************************
\\ ***************************************************************************
\\ ***************************************************************************

modular_inverse_3n_2 (u0, u1, u2, v0, v1, v2)={
\\ For u= x^3 + u2*x^2 + u1*x + u0 and v= v2*x^2 + v1*x + v0
\\ Compute d = gcd (u, v) and s
\\ such that d = (something)*u + s*v
\\ note that d is not i.g. the resultant(u,v).

\\ Cost  13M + 2SQ

\\ First pseudodivision quotient into q interpolation on
\\ 0, 1 for remainder r

    uval1 = u2 + u1 + u0;
    vval1 = v2 + v1 + v0;

    q1 = v2;
    q0 = v2 * u2 - v1;
    lc2 = v2^2;
    r0 = lc2 * u0 - q0 * v0;
    r1 = lc2 * (uval1 + 1) - (q1 + q0) * (vval1) - r0;
\\ Second pseudodivision quotient into Q, remainder into d
    Q1 = r1 * v2;
    Q0 = r1 * v1 - r0 * v2;
    lc2 = r1^2;
    d = lc2 * v0 - Q0 * r0;
\\ s = lc2 + Q q
    s0 = q0 * Q0;
    s2 = q1 * Q1;
    s1 = (q0 + q1) * (Q0 + Q1) - s0 - s2;
    s0 = s0 + lc2;

    return([s0,s1,s2,d]);
}

modular_inverse_3n_4 (u0, u1, u2, c0, c1, c2, c3, c4)={
\\ For u= x^3 + u2*x^2 + u1*x + u0 and c= c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
\\ Compute d = gcd (u, v) and s
\\ such that d = (something)*u + s*v
\\ note that d is not i.g. the resultant(u,v).
\\ 17M +2SQ

\\ Euclidian algorithm between c and u

\\ q1 = quotient of c by u
    q11 = c4;
    q10 = c3 - c4 * u2;
\\ r1 = c - q1 u by interpolation on 0, 1, -1;
    r10 = c0 - q10 * u0;
    uval1 = 1 + u2 + u1 + u0;
    uvalm1 = -1 + u2 - u1 + u0;
    val1 = c4 + c3 + c2 + c1 + c0 - (q11 + q10) * uval1;
    valm1 = c4 - c3 + c2 - c1 + c0 - (q10 - q11) * uvalm1;
    r11 = divconst (val1 - valm1, 2);
    r12 = val1 - r11 - r10;
\\ q2 = quotient of r12^2 u by r1;
    q21 = r12;
    q20 = r12 * u2 - r11;
\\ r2 = r12^2 u - q2 r1 by interpolation on 0 and 1;
    lc2 = r12^2;
    r20 = lc2 * u0 - q20 * r10;
    r21 = lc2 * uval1 - (q21 + q20) * (r12 + r11 + r10) - r20;
\\ q3 = quotient of r21^2 r1 by r2;
    q31 = r21 * r12;
    q30 = r21 * r11 - r20 * r12;
    lc2 = r21^2;
    d = lc2 * r10 - q30 * r20;
\\ s = r21^2 + q2 q3;
    s0 = q20 * q30;
    s2 = q21 * q31;
    s1 = (q21 + q20) * (q31 + q30) - s0 - s2;
    s0 = s0 + lc2;

    return([s0,s1,s2,d]);
}

modular_inverse_3n_3 (u0, u1, u2, c0, c1, c2, c3)={
\\ For u= x^3 + u2*x^2 + u1*x + u0 and c=  c3*x^3 + c2*x^2 + c1*x + c0
\\ Compute d = gcd (u, v) and s
\\ such that d = (something)*u + s*v
\\ note that d is not i.g. the resultant(u,v).

\\ 16M + 2SQ

\\ Euclidian algorithm between c and u

\\ q1 = quotient of c by u

    q10 = c3;
\\ r1 = c - q1 u by interpolation on 0, 1, -1;
    r10 = c0 - q10 * u0;
    uval1 = 1 + u2 + u1 + u0;
    uvalm1 = -1 + u2 - u1 + u0;
    val1 = c3 + c2 + c1 + c0 - q10 * uval1;
    valm1 = - c3 + c2 - c1 + c0 - q10 * uvalm1;
    r11 = divconst (val1 - valm1, 2);
    r12 = val1 - r11 - r10;
\\ q2 = quotient of r12^2 u by r1;
    q21 = r12;
    q20 = r12 * u2 - r11;
\\ r2 = r12^2 u - q2 r1 by interpolation on 0 and 1;
    lc2 = r12^2;
    r20 = lc2 * u0 - q20 * r10;
    r21 = lc2 * uval1 - (q21 + q20) * (r12 + r11 + r10) - r20;
\\ q3 = quotient of r21^2 r1 by r2;
    q31 = r21 * r12;
    q30 = r21 * r11 - r20 * r12;
    lc2 = r21^2;
    d = lc2 * r10 - q30 * r20;
\\ s = r21^2 + q2 q3;
    s0 = q20 * q30;
    s2 = q21 * q31;
    s1 = (q21 + q20) * (q31 + q30) - s0 - s2;
    s0 = s0 + lc2;

    return([s0,s1,s2,d]);
}

\\ ^^***********************************************************************
\\ ***************************************************************************
\\ ***************************************************************************

remainder_4_3n (t0, t1, t2, t3, t4, u0, u1, u2)={
\\ For u= x**3 + u2*x**2 + u1*x + u0 and t= t4*x**4 + t3*x**3 + t2*x**2 + t1*x + t0
\\ Compute Remainder(t,u)
\\ Cost 4M

    uval1  = u2 + u1 + u0;
    uvalm1 = u2 - u1 + u0;

\\ q = quotient of t by u
    q1 = t4;
    q0 = t3 - t4 * u2;
    \\ t = t - q u;
    T0 = t0 - q0 * u0;
    tval1  = (t4 + t3 + t2 + t1 + t0) - (q0 + q1) * (uval1 + 1);
    tvalm1 = (t4 - t3 + t2 - t1 + t0) - (q0 - q1) * (uvalm1 - 1);
    T1 = divconst (tval1 - tvalm1, 2);
    T2 = tval1 - T0 - T1;

    return([T0,T1,T2]);
}

\\ ***************************************************************************
\\ ***************************************************************************

quotient_6n_3n (u3, u4, u5, v0, v1, v2)={
\\ For u= x^6 + u5*x^5 + u4*x^4 + u3 + ... and v= x^3 + v2*x^2 + v1*x + v0
\\ or for u= x^9 + u5*x^8 + u4*x^7 + u3 * x^6 + ... and v= x^6 + v2*x^5 + v1*x^4 + v0*x^3
\\ Compute Quotient(u,v)
\\ Cost 3M

    q2 = u5 - v2;
    q1 = u4 - v1 - v2 * q2;
    r1 = u3 - v0 - v1 * (u5 - v2);
    q0 = - v2 * q1 + r1;

    return([q0,q1,q2]);
}

\\ Addition of divisors.

add(D1, D2, C)={

    u10 = D1[1];
    u11 = D1[2];
    u12 = D1[3];

    v10 = D1[4];
    v11 = D1[5];
    v12 = D1[6];


    u20 = D2[1];
    u21 = D2[2];
    u22 = D2[3];

    v20 = D2[4];
    v21 = D2[5];
    v22 = D2[6];

    f0 = C[1];
    f1 = C[2];
    f2 = C[3];
    h0 = C[4];
    h1 = C[5];
    h2 = C[6];
    h3 = C[7];


\\ step 1

\\  compute the inverse t1 of v1-v2 modulo u2
\\  13M + 2SQ

    u20u21_modular_inverse = modular_inverse_3n_2 (u20, u21, u22, v10-v20, v11-v21, v12-v22);

    t10 = u20u21_modular_inverse[1];
    t11 = u20u21_modular_inverse[2];
    t12 = u20u21_modular_inverse[3];
    res1  = u20u21_modular_inverse[4];

\\  compute c=(u1-u2)t1
\\  5M

    u10u20toom_cook_2_2 = toom_cook_2_2 (u10-u20,u11-u21,u12-u22,t10,t11,t12);

    c0 = u10u20toom_cook_2_2[1];
    c1 = u10u20toom_cook_2_2[2];
    c2 = u10u20toom_cook_2_2[3];
    c3 = u10u20toom_cook_2_2[4];
    c4 = u10u20toom_cook_2_2[5];

\\ compute the remainder r of c by u2
\\  4M
\\ Vorsicht={ mit den Zeichen r0,r1,r2={ die habe ich aus den anderen Implem. so uebernommen ...

    c0c1remainder_4_3 = remainder_4_3n(c0, c1, c2, c3, c4, u20, u21, u22);

    r2 = c0c1remainder_4_3[1];
    r1 = c0c1remainder_4_3[2];
    r0 = c0c1remainder_4_3[3];

\\ Solve the linear equation

\\  5M + 1SQ

    tmp1 = res1 * (v12+v22);
    tmp2 = v12^2;
    tmp3 = r0 * tmp2;
    tmp4 = r1 * tmp2;
    tmp5 = tmp3 * u22;
    gamma1_p =tmp1-(tmp4 - tmp5);
    Gamma1 = gamma1_p * res1;

\\  compute the part s of w= res1 r0 y^2 + s y + t with interpolation by 0 and 1:
\\  8M

    tmp6  = res1 * r0;
    tmp7  = r0 * tmp3;
    tmp8  = r2 * gamma1_p;
    tmp9  = tmp7 * u20;
    tmp10 = tmp6 * (v10 + v20);
    sval0 = tmp8 - tmp9 - tmp10;
    rval1  = r0 + r1 + r2;
    u2val1 = 1 + u22 + u21 + u20;
    v1val1 = v12 + v11 + v10;
    v2val1 = v22 + v21 + v20;
    tmp11 = rval1 * (tmp3 + gamma1_p);
    tmp12 = tmp7 * u2val1;
    tmp13 = tmp6 * (v1val1 + v2val1);

    sval1  = tmp11 - tmp12 - tmp13;

\\ s = s0 + s1 x
    s0 = sval0;
    s1 = sval1 - sval0;

\\ Compute Coeff(w,x,3) 
\\ 4M

    tmp14 = mulconst(tmp6 * v11, 2);
    tmp15 = v12 * (s1 + tmp14);
    tmp16 = tmp6 * tmp2;
    tmp17 = u12 * tmp16;
    coeff_w_x_3 = Gamma1 + tmp17 - tmp15;


\\ compute simultaneous inverse, such that the computing of the reultant should be easy and monic ...
\\ Cost :  if h3 == 0 ===> 2M + 1SQ + 1I
\\         if h3 ne 0 ===> 8M + 2SQ + 1I

    tmp18    = tmp6 * coeff_w_x_3;
    tmp19    = tmp6^2;

    if(h3 == 0,
        inv = invert(tmp18);
        coeffwx3_res1r0_m1 = inv;
        tmp20    = inv * coeff_w_x_3;
        res1r0m1 = tmp20,
    \\ else
        tmp20    = h3 * tmp19;
        tmp21    = h3 * ( -mulconst(tmp18,2) + tmp20);
        tmp22    = coeff_w_x_3^2;
        B        = tmp22 + tmp21;
        tmp23    = tmp18 * B;
        inv      = invert(tmp23);
        tmp24    = B * inv;
        coeffwx3_res1r0_m1 = tmp24;
        tmp25    = coeff_w_x_3 * coeffwx3_res1r0_m1;
        res1r0m1 = tmp25 ;
        tmp26    = inv * tmp18;
        Bm1      = tmp26;
        tmp27    = tmp22 * Bm1;
        Am1      = tmp27
    );
   \\ we use Am1 later, as we need to make the resltant monic.

\\  Compute sy = sy1 x + sy 0 ( = s /( res1 r0))
\\  2M

    tmp28 = s0 * res1r0m1;
    tmp29 = s1 * res1r0m1;
    sy0 = tmp28;
    sy1 = tmp29;

\\ Compute ty = u1 ( v12^2 x + gamma1) -v1 (v1 +sy) and tx = ty (res1 r0) ^2 ( res1 r0 coeff(w,x,3) )^(-1)

\\ 12M
    tmp30 = Gamma1 * res1r0m1;
    gamma1 =  tmp30;
    tmp31 = u10 * gamma1;
    tmp32 = v10 * (v10 + sy0);
    tmp33 = (1+ u12 + u11 + u10) * (tmp2 + gamma1);
    tmp34 = (v12 + v11 + v10) * (v12 + v11 + v10 + sy1 + sy0);
    tmp35 = (-1 + u12 - u11 + u10) * (- tmp2 + gamma1);
    tmp36 = (v12 - v11 + v10) * (v12 - v11 + v10 - sy1 + sy0);
    tmp37 = coeff_w_x_3 * res1r0m1;
    tyval0   = tmp31 - tmp32;
    tyval1   = tmp33 - tmp34;
    tyvalm1  = tmp35 - tmp36;
    tyvalinf = tmp37;

\\ for the following operation, no multiplication are necessary...
    ty0ty1ty2ty3 = interpolate(tyval0, tyval1, tyvalm1, tyvalinf);
    ty0 = ty0ty1ty2ty3[1];
    ty1 = ty0ty1ty2ty3[2];
    ty2 = ty0ty1ty2ty3[3];
    ty3 = ty0ty1ty2ty3[4];


    tmp38 = tmp19 * coeffwx3_res1r0_m1;
    tmp39 = ty0 * tmp38;
    tmp40 = ty1 * tmp38;
    tmp41 = ty2 * tmp38;
    tx0 = tmp39;
    tx1 = tmp40;
    tx2 = tmp41;

\\ Step 2:

\\ first coeff von tx^3 = x^9 + tx3_8 * x^8 + tx3_7 * x^7 + tx3_6 * x^6
\\ 1SQ + 1M

    tmp42 = tx2^2;
    tmp43 = tx2 * ( mulconst ( mulconst ( tx1,2 ) , 3 ) + tmp42 );
    tx3_8 = mulconst (tx2,3);
    tx3_7 = mulconst ( tx1 + tmp42, 3 );
    tx3_6 = tmp43 + mulconst (tx0,3);

\\ Berechnung von sy^2 = sy2_0 + sy2_1 * x + sy2_2 *x^2
\\  3SQ

    tmp44 = sy1^2;
    tmp45 = (sy1 + sy0)^2;
    tmp46 = sy0^2;


    sy2_0 = tmp46;
    sy2_1 = (tmp45 - tmp44 - tmp46);
    sy2_2 = tmp44;

\\ Computing of sy * ty = syty4 * x^4 + syty3 * x^3 + syty2 * x^2 + syty1 * x + syty0
\\ 5M
    sy0sy1toom_cook_1_3 = toom_cook_1_3 (sy0,sy1,ty0,ty1,ty2,ty3);

    syty0 = sy0sy1toom_cook_1_3[1];
    syty1 = sy0sy1toom_cook_1_3[2];
    syty2 = sy0sy1toom_cook_1_3[3];
    syty3 = sy0sy1toom_cook_1_3[4];
    syty4 = sy0sy1toom_cook_1_3[5];

\\ first coeff von sy^3 = sy3_3*x^3 + sy3_2 * x^2
\\ 2M

    tmp47 = sy1 * sy2_2;
    tmp48 = mulconst( sy2_2 * sy0, 3);
    sy3_3 = tmp47;
    sy3_2 = tmp48;


\\  Computing of H1= H18 * x^8  + H17 * x^7 + H16 * x^6 
\\  1M

    tmp49 = (1 - mulconst( syty4 ,3 )) * f2;
    H18 = 1 - mulconst( syty4 ,3 ) ;
    H17 = sy3_3 - mulconst( syty3 ,3 );
    H16 = f2 + sy3_2 - mulconst( syty2 ,3 ) + tmp49;

\\ Computing H2:

    if ( h0 == 0 && h1 == 0 && h2 == 0 && h3 == 0, 
\\ Picard curves case:
        H23=0;
        H24=0;
        H25=0;
        H26=0;
        H27=0;
        H28=0;
        H29=0,
        
        if ( h1 == 0 && h2 == 0 && h3 == 0,
            \\ h is a constant (h3=h2=h1=0 && h0 ne 9)
            \\ 1M + 1SQ
            tmp50 = - mulconst( ty3^2 , 2 );
            tmp51 = h0 * tmp50;
            H23=0;
            H24=0;
            H25=0;
            H26=tmp51;
            H27=0;
            H28=0;
            H29=0,

            if( h2 == 0 && h3 == 0,
                \\ h is linear
                \\ 6M

                ty2ty3toom_cook_1_1 = toom_cook_1_1 (ty2, ty3, sy2_2 - mulconst( ty2, 2 ), - mulconst( ty3, 2 ));
                P0 = ty2ty3toom_cook_1_1[1];
                P1 = ty2ty3toom_cook_1_1[2];
                P2 = ty2ty3toom_cook_1_1[3];
                
                h25h26h27toom_cook_1_1 = toom_cook_1_1 (h0, h1, P1 + sy1, P2);
                H25 = h25h26h27toom_cook_1_1[1];
                H26 = h25h26h27toom_cook_1_1[2];
                H27 = h25h26h27toom_cook_1_1[3];
                
                H23=0;
                H24=0;
                H28=0;
                H29=0,
                
                if(h3 == 0,
                    \\ h is quadratic
                    \\ 10M

                    P0P1P2toom_cook_2_2 = toom_cook_2_2 (ty1, ty2, ty3, sy2_1 - mulconst( ty1, 2 ) + h1, sy2_2 - mulconst( ty2, 2 ) + h2, - mulconst( ty3, 2 ));
                    
                    P0 = P0P1P2toom_cook_2_2[1];
                    P1 = P0P1P2toom_cook_2_2[2];
                    P2 = P0P1P2toom_cook_2_2[3];
                    P3 = P0P1P2toom_cook_2_2[4];
                    P4 = P0P1P2toom_cook_2_2[5];

                    H24H25H26toom_cook_2_2 = toom_cook_2_2(h0, h1, h2,  P2 +sy0,  P3 +sy1, P4);

                    H24 = H24H25H26toom_cook_2_2[1];
                    H25 = H24H25H26toom_cook_2_2[2];
                    H26 = H24H25H26toom_cook_2_2[3];
                    H27 = H24H25H26toom_cook_2_2[4];
                    H28 = H24H25H26toom_cook_2_2[5];
                    
                    H29=0,

                    \\ else
                        \\ h is cubic
                        \\ 15M

                        P0P1P2toom_cook_3_3 = toom_cook_3_3(ty0, ty1, ty2, ty3, sy2_0 - 2* ty0 + h0, sy2_1 - 2* ty1 + h1, sy2_2 - 2* ty2 + h2, - 2 * ty3 + h3);
                        P0 = P0P1P2toom_cook_3_3[1];
                        P1 = P0P1P2toom_cook_3_3[2];
                        P2 = P0P1P2toom_cook_3_3[3];
                        P3 = P0P1P2toom_cook_3_3[4];
                        P4 = P0P1P2toom_cook_3_3[5];
                        P5 = P0P1P2toom_cook_3_3[6];
                        P6 = P0P1P2toom_cook_3_3[7];
                        
                        tmp52 = f2 * sy1;
                        H23H24H25toom_cook_3_3 = toom_cook_3_3 (h0, h1, h2, h3, P3 + tmp52 , P4 + sy0,  P5 + sy1, P6);

                        H23 = H23H24H25toom_cook_3_3[1];
                        H24 = H23H24H25toom_cook_3_3[2];
                        H25 = H23H24H25toom_cook_3_3[3];
                        H26 = H23H24H25toom_cook_3_3[4];
                        H27 = H23H24H25toom_cook_3_3[5];
                        H28 = H23H24H25toom_cook_3_3[6];
                        H29 = H23H24H25toom_cook_3_3[7]
                )
            )
        )
    );


\\ Berechnung der Resultante:
\\ resWC= (k1^3 *( H1 + H2 ) + tx^3) * Am1 , notice that Am1 is equal to 1 if deg(h) < 3!!!


\\   1SQ  + 4M
    k1    = tmp38;
    tmp53 = k1^2;
    tmp54 = k1 * tmp53;
    tmp55 = tmp54 *( H16 + H26 );
    tmp56 = tmp54 *( H17 + H27 );
    tmp57 = tmp54 *( H18 + H28 );

    if(h3 == 0,
        ReswC6 = tmp55 + tx3_6;
        ReswC7 = tmp56 + tx3_7;
        ReswC8 = tmp57 + tx3_8;
        Am1=1,

            tmp58 = (tmp55 + tx3_6) * Am1;
            tmp59 = (tmp56 + tx3_7) * Am1;
            tmp60 = (tmp57 + tx3_8) * Am1;
            ReswC6 = tmp58;
            ReswC7 = tmp59;
            ReswC8 = tmp60
    );


\\ Computing the first term of U=u1*u2 = x^6 + U5*x^5 + U4*x^4 + U3 * x^3 + ....
\\ 3M
    tmp61 = u12 * u22;
    tmp62 = u12 * u21;
    tmp63 = u11 * u22;

    U5 = u12 + u22;
    U4 = u21 + u11 + tmp61;
    U3 = u20 + u10 + tmp62 + tmp63;

\\ 3M
    um0um1um2quotient_6n_3n = quotient_6n_3n (ReswC6, ReswC7, ReswC8, U3, U4, U5);

    Um0 = um0um1um2quotient_6n_3n[1];
    Um1 = um0um1um2quotient_6n_3n[2];
    Um2 = um0um1um2quotient_6n_3n[3];

\\  16M + 2SQ
    a10a11a12modular_inverse_3n_3 = modular_inverse_3n_3 (Um0, Um1, Um2, ty0 - sy2_0 - h0, ty1 - sy2_1 - h1, ty2 - sy2_2 - h2, ty3 - h3);

    alpha10 = a10a11a12modular_inverse_3n_3[1];
    alpha11 = a10a11a12modular_inverse_3n_3[2];
    alpha12 = a10a11a12modular_inverse_3n_3[3];
    res2 = a10a11a12modular_inverse_3n_3[4];

\\ 4M
    Rem10Rem11Rem12 = remainder_4_3n (syty0 - f0 , syty1 -f1, syty2 - f2, syty3, syty4 - 1 , Um0, Um1, Um2);
    
    Rem10 = Rem10Rem11Rem12[1];
    Rem11 = Rem10Rem11Rem12[2];
    Rem12 = Rem10Rem11Rem12[3];

\\  5M
    Rem20Rem21Rem22Rem23Rem24 = toom_cook_2_2 (Rem10,Rem11,Rem12,alpha10,alpha11,alpha12);
    Rem20 = Rem20Rem21Rem22Rem23Rem24[1];
    Rem21 = Rem20Rem21Rem22Rem23Rem24[2];
    Rem22 = Rem20Rem21Rem22Rem23Rem24[3];
    Rem23 = Rem20Rem21Rem22Rem23Rem24[4];
    Rem24 = Rem20Rem21Rem22Rem23Rem24[5];

\\ 4M
    fauxV12_0fauxV12_1fauxV12_2 = remainder_4_3n (Rem20 , Rem21, Rem22, Rem23, Rem24 , Um0, Um1, Um2);
    fauxV12_0 = fauxV12_0fauxV12_1fauxV12_2[1];
    fauxV12_1 = fauxV12_0fauxV12_1fauxV12_2[2];
    fauxV12_2 = fauxV12_0fauxV12_1fauxV12_2[3];


\\ 5M + 1I
    tmp64   = res2 *fauxV12_2;
    inv2    = invert(tmp64);
    tmp65   = inv2*fauxV12_2;
    invres2 = tmp65;


    vD1D2_0 = invres2 * fauxV12_0;
    vD1D2_1 = invres2 * fauxV12_1;
    vD1D2_2 = invres2 * fauxV12_2;

\\ uD1D2=((inv2*res2^2*vD1D2)^3+(inv2*res2^2)^3*UnivariatePolynomial(f4)-(inv2*res2^2)^3*h1*vD1D2^2+(inv2*res2^2)^3*h2*vD1D2) div u_D1D2 
\\  5M + 3SQ

    tmp66 = res2^2;
    tmp67 = inv2 * tmp66;
    tmp68 = tmp67^2;
    tmp69 = tmp68 * tmp67;
    tmp70 = tmp67 * vD1D2_1;
    tmp71 = tmp70^2;
    tmp72 = tmp67 * vD1D2_0;
    tmp73 = tmp70 * (tmp71 + mulconst ( mulconst ( tmp72 , 2), 3));

    if ( h0 == 0 && h1 == 0 && h2 == 0 && h3 == 0,
    \\ Picard curves case:

        tot5 = mulconst( tmp70, 3);
        tot4 = mulconst(tmp72 + tmp71, 3) - tmp69;
        tot3 = tmp73,


        if (h1 == 0 && h2 == 0 && h3 == 0,
            \\ h is a constant (h3=h2=h1=0 && h0 ne 9)

            tot5 = mulconst(tmp70, 3);
            tot4 = mulconst(tmp72+ tmp71, 3) - tmp69;
            tot3 = tmp73,

            if ( h2 == 0 && h3 == 0,
                \\ h is linear
                \\   2M
                tmp74 = h1 * vD1D2_2; 
                tmp75 = tmp69 * tmp74;  \\ coeff bei x^3

                tot5 = mulconst( tmp70, 3);
                tot4 = mulconst(tmp72 + tmp71, 3) - tmp69;
                tot3 = tmp73 + tmp75,

                if ( h3 == 0,
                    \\ h is quadratic
                    \\ 5M
                    tmp74tmp75tmp76 =toom_cook_1_1 (vD1D2_1,vD1D2_2,h1,h2);
                    tmp74 = tmp74tmp75tmp76[1];
                    tmp75 = tmp74tmp75tmp76[2];
                    tmp76 = tmp74tmp75tmp76[3];
                    
                    tmp77 = tmp69 * tmp75;  \\ coeff bei x^3
                    tmp78 = tmp69 * tmp76;  \\ coeff bei x^4

                    tot5 = mulconst( tmp70, 3);
                    tot4 = mulconst( tmp72 + tmp71, 3) - tmp69 + tmp78;
                    tot3 = tmp73 + tmp77,

                    \\ h is cubic
                    \\ 8M 
                        tmp74tmp75tmp76tmp77tmp78 = toom_cook_2_2 (vD1D2_0, vD1D2_1, vD1D2_2 , h1,h2, h3 );
                        tmp74 = tmp74tmp75tmp76tmp77tmp78[1];
                        tmp75 = tmp74tmp75tmp76tmp77tmp78[2];
                        tmp76 = tmp74tmp75tmp76tmp77tmp78[3];
                        tmp77 = tmp74tmp75tmp76tmp77tmp78[4];
                        tmp78 = tmp74tmp75tmp76tmp77tmp78[5];

                        tmp79 = tmp69 * tmp76;  \\ coeff bei x^3
                        tmp80 = tmp69 * tmp77;  \\ coeff bei x^4
                        tmp81 = tmp69 * tmp78;  \\ coeff bei x^5

                        tot5 = mulconst( tmp70, 3) + tmp81; 
                        tot4 = mulconst( tmp72+ tmp71, 3) - tmp69 +tmp80;
                        tot3 = tmp73 + tmp79
                )
            )
        )
    );
        
    

\\ 3M
    UD1D2_0UD1D2_1UD1D2_2 = quotient_6n_3n (tot3, tot4, tot5, Um0, Um1, Um2);
    UD1D2_0 = UD1D2_0UD1D2_1UD1D2_2[1];
    UD1D2_1 = UD1D2_0UD1D2_1UD1D2_2[2];
    UD1D2_2 = UD1D2_0UD1D2_1UD1D2_2[3];


    return([UD1D2_0, UD1D2_1, UD1D2_2, vD1D2_0, vD1D2_1, vD1D2_2]);
}


double (D1, C)=
{
    u10 = D1[1];
    u11 = D1[2];
    u12 = D1[3];
    v10 = D1[4];
    v11 = D1[5];
    v12 = D1[6];
    
    f0 = C[1];
    f1 = C[2];
    f2 = C[3];
    h0 = C[4];
    h1 = C[5];
    h2 = C[6];
    h3 = C[7];

    \\ step 1
    \\  compute v1^3 + v1 *h - f =  v1 * ( v1^2 + h ) - f

    \\ First compute v1^2

    v1val1 = v12 + v11 + v10;
    v1valm1 = v12 - v11 + v10;
    v1val2 = mulconst (v12, 4) + mulconst (v11, 2) + v10;

    \\ 5SQ
    v2_0v2_1v2_2v2_3v2_4 = toom_cook_square_2 (v10,v11,v12);
    v2_0 = v2_0v2_1v2_2v2_3v2_4[1];
    v2_1 = v2_0v2_1v2_2v2_3v2_4[2];
    v2_2 = v2_0v2_1v2_2v2_3v2_4[3];
    v2_3 = v2_0v2_1v2_2v2_3v2_4[4];
    v2_4 = v2_0v2_1v2_2v2_3v2_4[5];

    \\ Compute the 4 highest coefficients of
    \\ v1 * (v1^2 + h) - f into a variable b.
    \\ For this S (4, 3  4), perform first Toom-Cook on the 3 highest
    \\ coefficients.

    \\ 6M
    b1b2b3b4b5b6 = toom_cook_2_3 (v10, v11, v12, v2_1 + h1, v2_2 + h2, v2_3 + h3, v2_4 );
    b1 = b1b2b3b4b5b6[1];
    b2 = b1b2b3b4b5b6[2];
    b3 = b1b2b3b4b5b6[3];
    b4 = b1b2b3b4b5b6[4];
    b5 = b1b2b3b4b5b6[5];
    b6 = b1b2b3b4b5b6[6];


    \\ correct:

    b4=  b4 - 1;

    \\ Perform an exact division by u1 to obtain w1.
    \\ 6M

    w13 = b6;
    w12 = b5 - w13 * u12;
    w11 = b4 - w13 * u11 - w12 * u12;
    w10 = b3 - w13 * u10 - w12 * u11 - w11 * u12;

    \\  compute the inverse t1 of w1 modulo u1
    \\  with the Euclidian algorithm between c = w1 && u1
    \\  16M + 2SQ

    t10t11t12res1 = modular_inverse_3n_3 (u10, u11, u12, w10, w11, w12, w13);
    t10 = t10t11t12res1[1];
    t11 = t10t11t12res1[2];
    t12 = t10t11t12res1[3];
    res1 = t10t11t12res1[4];


    \\  compute Remainder ( 3 v1^2 + h, u1)
    \\  4M

    Rem1_0Rem1_1Rem1_2 = remainder_4_3n ( mulconst( v2_0,3) + h0, mulconst( v2_1, 3) + h1 , mulconst( v2_2, 3) + h2, mulconst( v2_3, 3) + h3, mulconst( v2_4,3), u10, u11, u12);
    Rem1_0 = Rem1_0Rem1_1Rem1_2[1];
    Rem1_1 = Rem1_0Rem1_1Rem1_2[2];
    Rem1_2 = Rem1_0Rem1_1Rem1_2[3];

    \\ 5M
    c0c1c2c3c4 = toom_cook_2_2 (Rem1_0, Rem1_1 , Rem1_2 , t10 , t11 , t12);
    c0 = c0c1c2c3c4[1];
    c1 = c0c1c2c3c4[2];
    c2 = c0c1c2c3c4[3];
    c3 = c0c1c2c3c4[4];
    c4 = c0c1c2c3c4[5];


    \\ compute the remainder r of c by u1
    \\  4M
    \\ Vorsicht: mit den Zeichen r0,r1,r2: die habe ich aus den &&eren Implem. so uebernommen ...

    r2r1r0 = remainder_4_3n (c0, c1, c2, c3, c4, u10, u11, u12);
    r2 = r2r1r0[1];
    r1 = r2r1r0[2];
    r0 = r2r1r0[3];


    \\ Solve the linear equation

    \\  5M

    tmp1 = mulconst( res1 * v12, 2);
    tmp2 = v2_4;
    tmp3 = r0 * tmp2;
    tmp4 = r1 * tmp2;
    tmp5 = tmp3 * u12;
    gamma1_p =tmp1-(tmp4 - tmp5);
    Gamma1 = gamma1_p * res1;

    \\  compute the part s of w= res1 r0 y^2 +s y + t with interpolation by 0 && 1:
    \\  8M

    tmp6  = res1 * r0;
    tmp7  = r0 * tmp3;
    tmp8  = r2 * gamma1_p;
    tmp9  = tmp7 * u10;
    tmp10 = mulconst( tmp6 * v10, 2);
    sval0 = tmp8 - tmp9 - tmp10;

    rval1  = r0 + r1 + r2;
    u1val1 = 1 + u12 + u11 + u10;
    v1val1 = v12 + v11 + v10;
    tmp11 = rval1 * (tmp3 + gamma1_p);
    tmp12 = tmp7 * u1val1;
    tmp13 = mulconst( tmp6 * v1val1, 2);
    sval1  = tmp11 - tmp12 - tmp13;

    \\ s = s0 + s1 x
    s0 = sval0;
    s1 = sval1 - sval0;

    \\ Compute Coeff(w,x,3) 
    \\ 4M

    tmp14 = mulconst(tmp6 * v11, 2);
    tmp15 = v12 * (s1 + tmp14);
    tmp16 = tmp6 * tmp2;
    tmp17 = u12 * tmp16;
    coeff_w_x_3 = Gamma1 + tmp17 - tmp15;


    \\ compute simultaneous inverse, such that the computing of the reultant should be easy and monic ...
    \\ Cost :  if h3 == 0 ===> 2M + 1SQ + 1I
    \\         if h3 ne 0 ===> 8M + 2SQ + 1I

    tmp18    = tmp6 * coeff_w_x_3;
    tmp19    = tmp6^2;

    if(h3 == 0,
        inv      = invert(tmp18);
        coeffwx3_res1r0_m1 = inv;
        tmp20    = inv * coeff_w_x_3;
        res1r0m1 = tmp20,
    \\ else:
            tmp20    = h3 * tmp19;
            tmp21    = h3 * ( -mulconst(tmp18,2) + tmp20);
            tmp22    = coeff_w_x_3^2;
            B        = tmp22 + tmp21;
            tmp23    = tmp18 * B;
            inv      = invert(tmp23);
            tmp24    = B * inv;
            coeffwx3_res1r0_m1 = tmp24;
            tmp25    = coeff_w_x_3 * coeffwx3_res1r0_m1;
            res1r0m1 = tmp25 ;
            tmp26    = inv * tmp18;
            Bm1      = tmp26;
            tmp27    = tmp22 * Bm1;
            Am1      = tmp27
    );
       \\ we use Am1 later, as we need to make the resultant monic.


    \\  Compute sy = sy1 x + sy 0 ( = s /( res1 r0))
    \\  2M

    tmp28 = s0 * res1r0m1;
    tmp29 = s1 * res1r0m1;
    sy0 = tmp28;
    sy1 = tmp29;

    \\ Compute ty = u1 ( v12^2 x + gamma1) -v1 (v1 +sy) and tx = ty (res1 r0) ^2 ( res1 r0 coeff(w,x,3) )^(-1)
    \\ 12M

    tmp30 = Gamma1 * res1r0m1;
    gamma1 =  tmp30;
    tmp31 = u10 * gamma1;
    tmp32 = v10 * (v10 + sy0);
    tmp33 = (1+ u12 + u11 + u10) * (tmp2 + gamma1);
    tmp34 = (v12 + v11 + v10) * (v12 + v11 + v10 + sy1 + sy0);
    tmp35 = (-1 + u12 - u11 + u10) * (- tmp2 + gamma1);
    tmp36 = (v12 - v11 + v10) * (v12 - v11 + v10 - sy1 + sy0);
    tmp37 = coeff_w_x_3 * res1r0m1;
    tyval0   = tmp31 - tmp32;
    tyval1   = tmp33 - tmp34;
    tyvalm1  = tmp35 - tmp36;
    tyvalinf = tmp37;

        \\ for the following operation, no multiplication are necessary...

    ty0ty1ty2ty3 = interpolate(tyval0, tyval1, tyvalm1, tyvalinf);
    ty0 = ty0ty1ty2ty3[1];
    ty1 = ty0ty1ty2ty3[2];
    ty2 = ty0ty1ty2ty3[3];
    ty3 = ty0ty1ty2ty3[4];



    tmp38 = tmp19 * coeffwx3_res1r0_m1;
    tmp39 = ty0 * tmp38;
    tmp40 = ty1 * tmp38;
    tmp41 = ty2 * tmp38;
    tx0 = tmp39;
    tx1 = tmp40;
    tx2 = tmp41;

        \\ Step 2:

        \\ first coeff von tx^3 = x^9 + tx3_8 * x^8 + tx3_7 * x^7 + tx3_6 * x^6
    \\ 1SQ + 1M

    tmp42 = tx2^2;
    tmp43 = tx2 * ( mulconst ( mulconst ( tx1,2 ) , 3 ) + tmp42 );
    tx3_8 = mulconst (tx2,3);
    tx3_7 = mulconst ( tx1 + tmp42, 3 );
    tx3_6 = tmp43 + mulconst (tx0,3);


    \\ Berechnung von sy^2 = sy2_0 + sy2_1 * x + sy2_2 *x^2
    \\  3SQ

    tmp44 = sy1^2;
    tmp45 = (sy1 + sy0)^2;
    tmp46 = sy0^2;


    sy2_0 = tmp46;
    sy2_1 = (tmp45 - tmp44 - tmp46);
    sy2_2 = tmp44;

    \\ Computing of sy * ty = syty4 * x^4 + syty3 * x^3 + syty2 * x^2 + syty1 * x + syty0
    \\ 5M
    syty0syty1syty2syty3syty4 = toom_cook_1_3 (sy0,sy1,ty0,ty1,ty2,ty3);
    syty0 = syty0syty1syty2syty3syty4[1];
    syty1 = syty0syty1syty2syty3syty4[2];
    syty2 = syty0syty1syty2syty3syty4[3];
    syty3 = syty0syty1syty2syty3syty4[4];
    syty4 = syty0syty1syty2syty3syty4[5];



    \\ first coeff von sy^3 = sy3_3*x^3 + sy3_2 * x^2
    \\ 2M

    tmp47 = sy1 * sy2_2;
    tmp48 = mulconst( sy2_2 * sy0, 3);
    sy3_3 = tmp47;
    sy3_2 = tmp48;


    \\  Computing of H1= H18 * x^8  + H17 * x^7 + H16 * x^6 
    \\  1M

    tmp49 = (1 - mulconst( syty4 ,3 )) * f2;

    H18 = 1 - mulconst( syty4 ,3 );
    H17 = sy3_3 - mulconst( syty3 ,3 );
    H16 = f2 + sy3_2 - mulconst( syty2 ,3 ) + tmp49;

    \\ Computing H2:


    if ( h0 == 0 && h1 == 0 && h2 == 0 && h3 == 0,
        \\ Picard curves case:
        H23=0;
        H24=0;
        H25=0;
        H26=0;
        H27=0;
        H28=0;
        H29=0,
        
        \\ else
        if ( h1 == 0 && h2 == 0 && h3 == 0,
            \\ h is a constant (h3=h2=h1=0 && h0 ne 9)
            \\ 1M + 1SQ
            tmp50 = - mulconst( ty3^2 , 2 );
            tmp51 = h0 * tmp50;
            H23=0;
            H24=0;
            H25=0;
            H26=tmp51;
            H27=0;
            H28=0;
            H29=0,
            
            \\ else
            if ( h2 == 0 && h3 == 0,
                    \\ h is linear
                    \\ 6M
                P0P1P2 = toom_cook_1_1 (ty2, ty3, sy2_2 - mulconst( ty2, 2 ), - mulconst( ty3, 2 ));
                P0 = P0P1P2[1];
                P1 = P0P1P2[2];
                P2 = P0P1P2[3];

                H25H26H27 = toom_cook_1_1 (h0, h1, P1 + sy1, P2);
                H25 = H25H26H27[1];
                H26 = H25H26H27[2];
                H27 = H25H26H27[3];

                H23=0;
                H24=0;
                H28=0;
                H29=0,
                
                if ( h3 == 0,
                        \\ h is quadratic
                        \\ 10M
                    P0P1P2P3P4 = toom_cook_2_2 (ty1, ty2, ty3, sy2_1 - mulconst( ty1, 2 ) + h1, sy2_2 - mulconst( ty2, 2 ) + h2, - mulconst( ty3, 2 ));
                    P0 = P0P1P2P3P4[1];
                    P1 = P0P1P2P3P4[2];
                    P2 = P0P1P2P3P4[3];
                    P3 = P0P1P2P3P4[4];
                    P4 = P0P1P2P3P4[5];

                    H24H25H26H27H28 = toom_cook_2_2 (h0, h1, h2,  P2 +sy0,  P3 +sy1, P4);
                    H24 = H24H25H26H27H28[1];
                    H25 = H24H25H26H27H28[2];
                    H26 = H24H25H26H27H28[3];
                    H27 = H24H25H26H27H28[4];
                    H28 = H24H25H26H27H28[5];

                    H29=0,
                    
                    \\ else
                        \\ h is cubic
                        \\ 15M
                        P0P1P2P3P4P5P6 = toom_cook_3_3 (ty0, ty1, ty2, ty3, sy2_0 - mulconst(ty0,2) + h0, sy2_1 - mulconst(ty1,2) + h1, sy2_2 - mulconst(ty2,2) + h2, - mulconst(ty3, 2) + h3);
                        
                        P0 = P0P1P2P3P4P5P6[1];
                        P1 = P0P1P2P3P4P5P6[2];
                        P2 = P0P1P2P3P4P5P6[3];
                        P3 = P0P1P2P3P4P5P6[4];
                        P4 = P0P1P2P3P4P5P6[5];
                        P5 = P0P1P2P3P4P5P6[6];
                        P6 = P0P1P2P3P4P5P6[7];

                        tmp52 = f2 * sy1;
                        H23H24H25H26H27H28H29 = toom_cook_3_3 (h0, h1, h2, h3, P3 + tmp52 , P4 + sy0,  P5 + sy1, P6);
                        H23 = H23H24H25H26H27H28H29[1];
                        H24 = H23H24H25H26H27H28H29[2];
                        H25 = H23H24H25H26H27H28H29[3];
                        H26 = H23H24H25H26H27H28H29[4];
                        H27 = H23H24H25H26H27H28H29[5];
                        H28 = H23H24H25H26H27H28H29[6];
                        H29 = H23H24H25H26H27H28H29[7]

                )
            )
        )
    );


    \\ Berechnung der Resultante:
    \\ resWC= (k1^3 *( H1 + H2 ) + tx^3) * Am1 , notice that Am1 is ==ual to 1 if deg(h) < 3!!!

    \\   4M + 1SQ

    k1    = tmp38;
    tmp53 = k1^2;
    tmp54 = k1 * tmp53;
    tmp55 = tmp54 *( H16 + H26 );
    tmp56 = tmp54 *( H17 + H27 );
    tmp57 = tmp54 *( H18 + H28 );

    if (h3 == 0,
        ReswC6 = tmp55 + tx3_6;
        ReswC7 = tmp56 + tx3_7;
        ReswC8 = tmp57 + tx3_8;
        Am1=1,

        \\ else
            \\ 3M
            tmp58 = (tmp55 + tx3_6) * Am1;
            tmp59 = (tmp56 + tx3_7) * Am1;
            tmp60 = (tmp57 + tx3_8) * Am1;
            ReswC6 = tmp58;
            ReswC7 = tmp59;
            ReswC8 = tmp60
    );


    \\ Computing the first term of U=u1*u2 = x^6 + U5*x^5 + U4*x^4 + U3 * x^3 + ....

    \\ 1M + 1SQ

    tmp61 = u12^2;
    tmp62 = u12 * u11;
    tmp63 = tmp62;
    U5 = mulconst(u12, 2);
    U4 = mulconst(u11,2) + tmp61;
    U3 = mulconst(u10,2) + tmp62 + tmp63;

    \\ 3M
    Um0Um1Um2 = quotient_6n_3n (ReswC6, ReswC7, ReswC8, U3, U4, U5);
    Um0 = Um0Um1Um2[1];
    Um1 = Um0Um1Um2[2];
    Um2 = Um0Um1Um2[3];


    \\ 16M + 2SQ
    alpha10alpha11alpha12res2 = modular_inverse_3n_3 (Um0, Um1, Um2, ty0 - sy2_0 - h0, ty1 - sy2_1 - h1, ty2 - sy2_2 - h2, ty3 - h3);
    alpha10 = alpha10alpha11alpha12res2[1];
    alpha11 = alpha10alpha11alpha12res2[2];
    alpha12 = alpha10alpha11alpha12res2[3];
    res2 = alpha10alpha11alpha12res2[4];


    \\ 4M
    Rem10Rem11Rem12 = remainder_4_3n (syty0 - f0 , syty1 -f1, syty2 - f2, syty3, syty4 - 1 , Um0, Um1, Um2);
    Rem10 = Rem10Rem11Rem12[1];
    Rem11 = Rem10Rem11Rem12[2];
    Rem12 = Rem10Rem11Rem12[3];


    \\ 5M
    Rem20Rem21Rem22Rem23Rem24= toom_cook_2_2 (Rem10,Rem11,Rem12,alpha10,alpha11,alpha12);
    Rem20 = Rem20Rem21Rem22Rem23Rem24[1];
    Rem21 = Rem20Rem21Rem22Rem23Rem24[2];
    Rem22 = Rem20Rem21Rem22Rem23Rem24[3];
    Rem23 = Rem20Rem21Rem22Rem23Rem24[4];
    Rem24 = Rem20Rem21Rem22Rem23Rem24[5];


    \\ 4M
    fauxV12_0fauxV12_1fauxV12_2 = remainder_4_3n (Rem20 , Rem21, Rem22, Rem23, Rem24 , Um0, Um1, Um2);
    fauxV12_0 = fauxV12_0fauxV12_1fauxV12_2[1];
    fauxV12_1 = fauxV12_0fauxV12_1fauxV12_2[2];
    fauxV12_2 = fauxV12_0fauxV12_1fauxV12_2[3];


    \\ 5M + 1I

    tmp64   = (res2 *fauxV12_2);
    inv2    = invert(tmp64);
    tmp65   = inv2*fauxV12_2;
    invres2 = tmp65;


    vD1D2_0 = invres2 * fauxV12_0;
    vD1D2_1 = invres2 * fauxV12_1;
    vD1D2_2 = invres2 * fauxV12_2;

    \\ uD1D2=((inv2*res2^2*vD1D2)^3+(inv2*res2^2)^3*UnivariatePolynomial(f4)-(inv2*res2^2)^3*h1*vD1D2^2+(inv2*res2^2)^3*h2*vD1D2) div u_D1D2 

    \\  5M + 3SQ

    tmp66 = res2^2;
    tmp67 = inv2*tmp66;
    tmp68 = tmp67^2;
    tmp69 = tmp68*tmp67;
    tmp70 = tmp67 * vD1D2_1;
    tmp71 = tmp70^2;
    tmp72 = tmp67 * vD1D2_0;
    tmp73 = tmp70 * (tmp71 + mulconst ( mulconst ( tmp72 , 2), 3));

    if ( h0 == 0 && h1 == 0 && h2 == 0 && h3 == 0,
            \\ Picard curves case:
        tot5 = mulconst( tmp70, 3);
        tot4 = mulconst(tmp72 + tmp71, 3) - tmp69;
        tot3 = tmp73,

        if( h1 == 0 && h2 == 0 && h3 == 0,
                \\ h is a constant (h3=h2=h1=0 && h0 ne 9)
            tot5 = mulconst(tmp70, 3);
            tot4 = mulconst(tmp72+ tmp71, 3) - tmp69;
            tot3 = tmp73,
            
            if ( h2 == 0 && h3 == 0,
                    \\ h is linear
                    \\ 2M
                tmp74 = h1 * vD1D2_2;
                tmp75 = tmp69 * tmp74;      \\ coeff bei x^3

                tot5 = mulconst(tmp70,3);
                tot4 = mulconst(tmp72+ tmp71, 3) - tmp69;
                tot3 = tmp73 + tmp75,

                if ( h3 == 0,
                        \\ h is quadratic
                        \\ 5M
                    tmp74tmp75tmp76 = toom_cook_1_1 (vD1D2_1,vD1D2_2,h1,h2);
                    tmp74 = tmp74tmp75tmp76[1];
                    tmp75 = tmp74tmp75tmp76[2];
                    tmp76 = tmp74tmp75tmp76[3];

                    tmp77 = tmp69 * tmp75;      \\ coeff bei x^3
                    tmp78 = tmp69 * tmp76;     \\ coeff bei x^4

                    tot5 = mulconst( tmp70, 3);
                    tot4 = mulconst( tmp72+ tmp71, 3) - tmp69 + tmp78;
                    tot3 = tmp73 + tmp77,

                            \\ h is cubic
                            \\ 8M
                        tmp74tmp75tmp76tmp77tmp78 = toom_cook_2_2 (vD1D2_0, vD1D2_1, vD1D2_2 , h1,h2, h3 );
                        tmp74 = tmp74tmp75tmp76tmp77tmp78[1];
                        tmp75 = tmp74tmp75tmp76tmp77tmp78[2];
                        tmp76 = tmp74tmp75tmp76tmp77tmp78[3];
                        tmp77 = tmp74tmp75tmp76tmp77tmp78[4];
                        tmp78 = tmp74tmp75tmp76tmp77tmp78[5];

                        tmp79 = tmp69 * tmp76;      \\ coeff bei x^3
                        tmp80 = tmp69 * tmp77;      \\ coeff bei x^4
                        tmp81 = tmp69 * tmp78;      \\ coeff bei x^5

                        tot5 = mulconst( tmp70, 3) + tmp81;
                        tot4 = mulconst( tmp72+ tmp71, 3) - tmp69 +tmp80;
                        tot3 = tmp73 + tmp79
                )
            )
        )
    );


    \\ 3M
    UD1D2_0UD1D2_1UD1D2_2 = quotient_6n_3n (tot3, tot4, tot5, Um0, Um1, Um2);
    UD1D2_0 = UD1D2_0UD1D2_1UD1D2_2[1];
    UD1D2_1 = UD1D2_0UD1D2_1UD1D2_2[2];
    UD1D2_2 = UD1D2_0UD1D2_1UD1D2_2[3];


    return([UD1D2_0, UD1D2_1, UD1D2_2, vD1D2_0, vD1D2_1, vD1D2_2]);
}


makeRandomDiv(C,p)=
{
    uv = [0,0];
    uprod = 1;
    finished=False;
    while(finished==False, a = Mod(random(p),p); while(a==0, a = Mod(random(p),p);); b = Mod(random(p),p); c = Mod(random(p),p); bigu = x^6*a^3 + 3*x^5*a^2*b + 3*x^4*a*b^2 + 3*x^4*a^2*c + x^3*b^3 + 6*x^3*a*b*c + 3*x^2*b^2*c + 3*x^2*a*c^2 + 3*x*b*c^2 + c^3 - x^4 - C[3]*x^2 - C[2]*x - C[1]; factors = polrootsmod(bigu,p); if(length(factors)>=3, for(i=1,3, uprod*=(x-factors[i]); uv = [polcoeff(uprod,0),polcoeff(uprod,1),polcoeff(uprod,2),c,b,a]; finished=True;); uprod=1;););
    return(uv);
}

eqZero(N_1D,D,C,p)=
{
    Dinv = [0,0,0,0,0,0];
    Dinv[4] = D[4];
    Dinv[5] = D[5];
    Dinv[6] = D[6];
    a = D[6];
    b = D[5];
    c = D[4];
    B = isTypical(N_1D,C,p);
    
    subarr1 = [Dinv[4],Dinv[5],Dinv[6]];
    subarr2 = [N_1D[4],N_1D[5],N_1D[6]];
    
    if(subarr1 != subarr2 && B, return False;);
    
    bigu = x^6*a^3 + 3*x^5*a^2*b + 3*x^4*a*b^2 + 3*x^4*a^2*c + x^3*b^3 + 6*x^3*a*b*c + 3*x^2*b^2*c + 3*x^2*a*c^2 + 3*x*b*c^2 + c^3 - x^4 - C[3]*x^2 - C[2]*x - C[1];

    quopoly = 1/a^3 * bigu/(x^3 + D[3]*x^2 + D[2]*x + D[1]);
    
    Dinv[1] = polcoeff(quopoly,0);
    Dinv[2] = polcoeff(quopoly,1);
    Dinv[3] = polcoeff(quopoly,2);
    
    if(Dinv == N_1D && B, return True;);
    if(Dinv == N_1D && B==False, return True;);
    if(Dinv != N_1D && B, return False;);
    if(Dinv != N_1D && B==False, return True;);
}
isTypical(D,C,p)=
{
    polyu = x^3 + Mod(D[3],p)*x^2 + Mod(D[2],p)*x + Mod(D[1],p);
    polyv = Mod(D[6],p)*x^2 + Mod(D[5],p)*x + Mod(D[4],p);
    curve = polyv^3 - (C[1]+C[2]*x+C[3]*x^2+x^4);
    return (curve%polyu == 0);
}

Dpowers(l,D,C)=
{
    Dpows = List();
    listput(Dpows,D);
    for(i=1,l, Dtemp = Dpows[length(Dpows)]; D = double(Dtemp, C); listput(Dpows,D););
    return(Dpows);
}

multbyN(N,D,Dpows,C)=
{
    binaryN = binary(N);
    l = length(binaryN);
    Dtemp = Dpows[l];
    for(i=2,l, 
        if(binaryN[i]==1, 
           Dtemp = add(Dpows[l+1-i],Dtemp,C);
        )
    );
    return(Dtemp);
}

FermatPrimes = [5,17,257,65537];
nf = bnfinit(z^2+3);
zzeta = nfeltdiv(nf, nfeltadd(nf, -1, z), 2);

lift2mod3(p,cartierManin,C)=
{
    cartierManin = matrix(3,3,a,b,cartierManin[a][b]);
    A2modp = lift(Mod(-cartierManin[3,1]*cartierManin[1,3]-cartierManin[3,2]*cartierManin[2,3],p));

    jacvals = [0,0];

    f = x^4 + Mod(C[3],p)*x^2 + Mod(C[2],p)*x + Mod(C[1],p);

    mod3res = polisirreducible(f);

    candidates = 0;

    for(i=-1,2,
        A2 = A2modp + i*p;
        if(Mod(A2,3)==mod3res, 
            Jac = 1+A2+A2*p+p^3;
            candidates += 1;
            jacvals[candidates]=[Jac, [1,0,A2,0,p*A2,0,p^3]]
        )
    );

    numruns = 0;

    l = length(binary(p+1))*3+1;

    while(candidates>1,

        D = makeRandomDiv(C,p); 

        Dpows = Dpowers(l,D,C);

        numruns+=1; 
        if(numruns > 10 && (p in FermatPrimes), print("Use O(p^1/4) algorithm instead."); candidates=0;); 

        \\ Use a Monte-Carlo algorithm; this is more efficient in practice. This is correct with probability >=1-2^-10.

        if(numruns > 10 && jacvals[2][1]==(p+1)^3, jacvals[1]=0; candidates-=1;);

        Jac = jacvals[1]; 
        N_1D = multbyN(Jac[1]-1,D,Dpows,C);

        if(!eqZero(N_1D,D,C,p), 
           jacvals[1]=0;
           candidates-=1, 
           D2 = multbyN(3*(p+1)*p, D,Dpows,C);
           if(!eqZero(add(N_1D,D2,C),D,C,p), 
              jacvals[2]=0;
              candidates-=1
            )
        )
    );
    return(jacvals);
}

lift1mod3(p,cartierManin,t)=
{
    cartierManin = matrix(3,3,a,b,cartierManin[a][b]);
    
    v = sqrtint((4*p-t^2)/12);
    sqrt_3 = Mod(t/(2*v),p);
    zeta3 = (Mod(-1,p)+ sqrt_3)/(Mod(2,p));

    zeta3bar = zeta3^2;

    M = [Mod(1,p), zeta3; Mod(1,p), zeta3bar];
    r_p = Mod(cartierManin[1,1]+cartierManin[2,2], p);
    r_pbar = Mod(cartierManin[3,3], p);

    a_p = matsolve(M, Col([r_p, r_pbar]));
    a_p = [centerlift(a_p[1]), centerlift(a_p[2])];

    zeta6 = nfeltadd(nf, zzeta, 1);

    pi_ideal_element = nfeltadd(nf, -1*lift(zeta3), zzeta);
    pi_ideal_1 = idealhnf(nf, pi_ideal_element);
    pi_ideal_2 = idealhnf(nf, p);
    pi_ideal = idealadd(nf, pi_ideal_1, pi_ideal_2);

    pi_00 = bnfisprincipal(nf, pi_ideal)[2];
    pibar_00 = nfeltadd(nf, pi_00[1]-pi_00[2], -1*pi_00[2]*zzeta);

    Kp = [p, pi_00,1,1,[pibar_00[1],-pibar_00[2]; pi_00[2],-pi_00[1]]];

    sigma_zeta6 = nfmodpr(nf, zeta6, Kp);
    sigma_pibar = nfmodpr(nf, pibar_00, Kp);

    s = cartierManin[1,1] *  cartierManin[2,2] - cartierManin[1,2]*cartierManin[2,1];
    s_p = nfmodpr(nf, s, Kp);

    r_pbar_ffelt = nfmodpr(nf, cartierManin[3,3], Kp);
    sextic_val = s_p / (r_pbar_ffelt * sigma_pibar );

    my(sextic_power = 0);
    sextic_power = 0;

    for(i = 0, 5, if(sigma_zeta6^i == sextic_val, sextic_power = i));
    
    if(a_p[2]!=0,
       a_pbar = nfeltadd(nf, lift(a_p[1]-a_p[2]), -zzeta*lift(a_p[2]));
       a_p = nfeltadd(nf, lift(a_p[1]), lift(a_p[2])*zzeta);
    );
    if(a_p[2]==0,
       a_pbar = Col([a_p[1];0]);
       a_p = a_pbar;
    );
    c_p = nfeltmul(nf, nfeltmul(nf, nfeltpow(nf, zeta6, sextic_power), pibar_00), p);
    c_pbar = nfeltadd(nf, c_p[1]-c_p[2], -1*c_p[2]*zzeta);
    b_p = nfeltmul(nf, c_p, nfeltdiv(nf, a_pbar, p));
    b_pbar = nfeltadd(nf, b_p[1]-b_p[2], -1*b_p[2]*zzeta);
    
    A1 = -2*(a_p[1])+a_p[2];
    A2 = nfeltadd(nf, nfeltmul(nf, a_p, a_pbar), nfeltadd(nf, b_p, b_pbar));
    A3 = -nfeltadd(nf, nfeltadd(nf, c_p, c_pbar), nfeltadd(nf, nfeltmul(nf, a_p, b_pbar), nfeltmul(nf, a_pbar, b_p)));
    return([A1,A2,A3,p]);
}

liftLpolys(CMList,TrList,C)=
{
    l = length(CMList);
    for(i=30,l,
        p = CMList[i][1];
        CM = CMList[i][2];
        if(Mod(p,3)==1, t = TrList[i]; print(lift1mod3(p,CM,t)););
        if(Mod(p,3)==2, print(lift2mod3(p,CM,C)););
    )
}