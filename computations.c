// List of useful procedures
< "procedures.c";
< "search_conics.c"


// Section 2

ring RF = (0,e),(x,y,z),dp;
minpoly = e4 + 1;

// LF - maximal tangency lines to the Fermat quartic
list LF;
LF[1] = x + (e3)*y;
LF[2] = x + (-e)*y;
LF[3] = x + (e)*y;
LF[4] = x + (-e3)*y;
LF[5] = x + (e3)*z;
LF[6] = x + (-e)*z;
LF[7] = x + (e)*z;
LF[8] = x + (-e3)*z;
LF[9] = y + (e3)*z;
LF[10] = y + (-e)*z;
LF[11] = y + (e)*z;
LF[12] = y + (-e3)*z;

// Computing the t-vector
t_vector(LF);

// Ideal I of 48 double points and its free resolution
ideal I = intersect_all(P_operator(2, LF));
resolution Ri = mres(I, 0);
print(Ri);

// PF - maximal tangency points to the Fermat quartic
list PF = dualize_lines(LF);

// Ideal of points PF as a complete intersection
mstd(intersect_all(PF))[2];


// Section 3

ring RK = (0,i),(x,y,z),dp;
minpoly = i2 + 1;

// LK - maximal tangency lines to the Komiya-Kuribayashi quartic
list LK;
LK[1] = -2i*x - y + z;
LK[2] = -2i*x + y + z;
LK[3] = 2i*x - y + z;
LK[4] = 2i*x + y + z;
LK[5] = -2i*y - x + z;
LK[6] = -2i*y + x + z;
LK[7] = 2i*y - x + z;
LK[8] = 2i*y + x + z;
LK[9] = -2i*z - y + x;
LK[10] = -2i*z + y + x;
LK[11] = 2i*z - y + x;
LK[12] = 2i*z + y + x;

// Intersections of LK
t_vector(LK);

// Ideal J of 66 double points and its free resolution
ideal J = intersect_all(P_operator(2, LK));
resolution Rj = mres(J, 0);
print(Rj);

// PK - maximal tangency points to the Komiya-Kuribayashi quartic
list PK;
PK[1] = make_point(i, 1, -1);
PK[2] = make_point(-i, 1, 1);
PK[3] = make_point(-i, 1, -1);
PK[4] = make_point(i, 1, 1);
PK[5] = make_point(-1, -i, 1);
PK[6] = make_point(1, -i, 1);
PK[7] = make_point(-1, i, 1);
PK[8] = make_point(1, i, 1);
PK[9] = make_point(1, -1, -i);
PK[10] = make_point(1, 1, -i);
PK[11] = make_point(1, -1, i);
PK[12] = make_point(1, 1, i);

// Lines dual to PK
list PK' = dualize_points(PK);
t_vector(PK');

// Ideal of points of multiplicity 4 as a complete intersection
mstd(intersect_all(P_operator(4, PK')))[2];

// Ideal K of 30 double points and its free resolution
ideal K = intersect_all(P_operator(2, PK'));
resolution Rk = mres(K, 0);
print(Rk);


// Section 4

// Polynomial Q
poly Q = 1056x12 - 19278x10y2 - 75207x8y4 - 111042x6y6 - 75207x4y8 - 19278x2y10 + 1056y12 - 19278x10z2 - 137198x8y2z2 - 287194x6y4z2 - 287194x4y6z2 - 137198x2y8z2 - 19278y10z2 - 75207x8z4 - 287194x6y2z4 - 413110x4y4z4 - 287194x2y6z4 - 75207y8z4 - 111042x6z6 - 287194x4y2z6 - 287194x2y4z6 - 111042y6z6 - 75207x4z8 - 137198x2y2z8 - 75207y4z8 - 19278x2z10 - 19278y2z10 + 1056z12;

// Saturate MTPs from Q
ideal I = sat(ideal(Q, x4 + y4 + z4 + 3x2y2 + 3y2z2 + 3z2x2), intersect_all(PK))[1];
mstd(I)[2];


// Section 5

// Find the conics for sextactic points on the Fermat quartic
ring RFc = (0,u),(x,y,z,a,b,c,d,e,f),dp;
minpoly = u8 + 4u6 + 2u4 + 28u2 + 1;

number s = 5/24u7+19/24u5+5/24u3+151/24u; // 2^(1/4)
number i = u - s;
number q = (1/2)*(s^2)*(i+1); // 8th root of unity

// List of sextactic points on the Fermat quartic
list S;
for(int j = 1; j <= 7; j = j + 2)
{
	for(int k = 1; k <= 7; k = k + 2)
	{
		S = S + list(1, q^(j-1), s*q^k);
		S = S + list(1, s*q^j, q^(k-1));
		S = S + list(1, q^j/s, q^k/s);
	}
}

// List of conics with pairs of points
list CON = find_conics(x4 + y4 + z4, S);

// First group of 24 intersection points
ring R1 = (0,u),(x,y,z),dp;
minpoly = u4 - 4u2 + 16;

number i = (1/8)*u3;
number p = u - i; // sqrt(3)

list P;
for(int j = 1; j <= 4; j++)
{
	P = P + list(make_point(2*i^j, i*p+1, 0));
	P = P + list(make_point(i*p+1, 2*i^j, 0));
	P = P + list(make_point(2*i^j, 0, i*p+1));
	P = P + list(make_point(i*p+1, 0, 2*i^j));
	P = P + list(make_point(0, 2*i^j, i*p+1));
	P = P + list(make_point(0, i*p+1, 2*i^j));
}
ideal I = intersect_all(P);
mstd(I)[2];

// Second group of 24 intersection points
ring R2 = (0,u),(x,y,z),dp;
minpoly = u8 + 4u6 + 2u4 + 28u2 + 1;

number s = 5/24u7+19/24u5+5/24u3+151/24u; // 2^(1/4)
number i = u - s;
number q = (1/2)*(s^2)*(i+1); // 8th root of unity

list P;
for(int j = 1; j <= 4; j++)
{
	P = P + list(make_point(i^j, q*s, 0));
	P = P + list(make_point(q*s, i^j, 0));
	P = P + list(make_point(i^j, 0, q*s));
	P = P + list(make_point(q*s, 0, i^j));
	P = P + list(make_point(0, i^j, q*s));
	P = P + list(make_point(0, q*s, i^j));
}
ideal I = intersect_all(P);
mstd(I)[2];


// Find the conics for sextactic points on the Komiya-Kuribayashi quartic
ring RKc = (0,u),(x,y,z,a,b,c,d,e,f),dp;
minpoly = u4 - 8u2 + 36;

number i = 1/12u3 - 1/6u;
number p = u - i; // sqrt(5)

// List of sextactic points on the Komiya-Kuribayashi quartic
list S;
S = S + list(0, 2, i*(p+1));
S = S + list(0, 2, -i*(p+1));
S = S + list(i*(p+1), 0, 2);
S = S + list(-i*(p+1), 0, 2);
S = S + list(2, i*(p+1), 0);
S = S + list(2, -i*(p+1), 0);
S = S + list(0, 2, i*(p-1));
S = S + list(0, 2, -i*(p-1));
S = S + list(i*(p-1), 0, 2);
S = S + list(-i*(p-1), 0, 2);
S = S + list(2, i*(p-1), 0);
S = S + list(2, -i*(p-1), 0);
S = S + list(1, 1, i*p);
S = S + list(1, 1, -i*p);
S = S + list(1, i*p, 1);
S = S + list(1, -i*p, 1);
S = S + list(i*p, 1, 1);
S = S + list(-i*p, 1, 1);
S = S + list(1, -1, i*p);
S = S + list(1, -1, -i*p);
S = S + list(-1, i*p, 1);
S = S + list(-1, -i*p, 1);
S = S + list(i*p, 1, -1);
S = S + list(-i*p, 1, -1);

// List of conics with pairs of points
list CON = find_conics(x4 + y4 + z4 + 3x2y2 + 3y2z2 + 3z2x2, S);
