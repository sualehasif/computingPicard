Ha := Sp(6,2)![1,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0];
Hb := Sp(6,2)![0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,0];
Gz := Sp(6,2)![1,1,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1];
Z3 := Sp(6,2)![1,1,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0];

H := sub<Sp(6,2)|Ha,Hb>;
assert IdentifyGroup(H) eq <648,533>;
G:=sub<Sp(6,2)|H,Gz>;
assert IdentifyGroup(G) eq <1296,2891>;
assert #[H: H in Subgroups(Sp(6,2):OrderEqual:=648)|IdentifyGroup(H`subgroup) eq <648,533>] eq 1;
assert #[H: H in Subgroups(Sp(6,2):OrderEqual:=1296)|IdentifyGroup(H`subgroup) eq <1296,2891>] eq 1;
{<CharacteristicPolynomial(g),Order(g)>:g in G|not g in H};
PG := quo<G|Z3>;
assert IdentifyGroup(PG) eq <432,734>;

T9 := TransitiveGroup(9,26);
assert IsIsomorphic(PG,T9);
assert #Subgroups(T9:IndexEqual:=2) eq 1;
H9 := Subgroups(T9:IndexEqual:=2)[1]`subgroup;

print {* CycleStructure(g):g in T9|not g in H9 *};

F<f0,f1,f2>:=FunctionField(Rationals(),3);
Rus<u0,u1,s0,s1>:=PolynomialRing(F,4);
I:=[-u1^3 + 2*s1, -3*u0*u1^2 + 2*s0 + s1^2 - f2, -3*u0^2*u1 + 2*s0*s1 - f1, -u0^3 + s0^2 - f0];
B:=GroebnerBasis(I);
Rx<x> := PolynomialRing(F);
psi := Evaluate(B[#B],[0,0,0,x]);
RF<f0,f1,f2>:=PolynomialRing(Integers(),3);
assert IsSquare(-3*Numerator(Discriminant(psi))); 

function a2parity(f0,f1,f2)
    p := Characteristic(Parent(f0));
    assert p gt 3 and p mod 3 eq 2;
    R := PolynomialRing(Parent(f0));
    h := R![Evaluate(c,[f0,f1,f2]):c in Coefficients(psi)];
    q := #Parent(f0);
    S<x> := quo<R|h>;
    return Degree(GCD(h,R!(x^q-x))) gt 0 select 1 else 0;
end function;

function a2value(f0,f1,f2)
    p := Characteristic(Parent(f0));
    assert p gt 3;
    P<x,y,z> := ProjectiveSpace(GF(p),2);
    C := Curve(P,y^3*z-x^4-f2*x^2*z^2-f1*x*z^3-f0*z^4);
    if IsSingular(C) then return 1000*p,false; end if;
    f := LPolynomial(C);
    if Degree(f) ne 6 then return 1000*p,false; end if;
    return Coefficient(f,2),true;
end function;

for f0,f1,f2 in GF(5) do a2,sts := a2value(f0,f1,f2); if sts then assert a2 mod 2 eq a2parity(f0,f1,f2); end if; end for;
for p in PrimesInInterval(5,30) do
    if p mod 3 ne 2 then continue; end if;
    F:=GF(p);
    for i:=1 to 100 do
        f0:= Random(F); f1:=Random(F); f2:=Random(F);
        a2,sts := a2value(f0,f1,f2);
        if sts then assert a2 mod 2 eq a2parity(f0,f1,f2); end if;
    end for;
    printf "Correctly computed a2 parity for 100 random curves over F_%o\n",p;
end for;

time printf "Distribution of a2-parities for y^3=x^4+x+1 over p < 2^16: %o ", {* a2parity(GF(p)!1,GF(p)!1,GF(p)!0):p in PrimesInInterval(5,2^16) | p mod 3 eq 2 *};
time printf "Distribution of a2-parities for y^3=x^4+x+1 over p < 2^20: %o ", {* a2parity(GF(p)!1,GF(p)!1,GF(p)!0):p in PrimesInInterval(5,2^20) | p mod 3 eq 2 *};
time printf "Distribution of a2-parities for y^3=x^4+x+1 over p < 2^24: %o ", {* a2parity(GF(p)!1,GF(p)!1,GF(p)!0):p in PrimesInInterval(5,2^24) | p mod 3 eq 2 *};
