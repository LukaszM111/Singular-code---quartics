// Coefficient at a given monomial
proc coef_at(poly K, poly vars, poly m)
{
	matrix M = coef(K, vars);
	poly c = 0;
	for(int i = 1; i <= ncols(M); i++)
	{
		if(M[1,i] == m)
		{
			c = M[2,i];
			break;
		}
	}
	return(c);
}

// All monomials of variables m1 and m2, of degree k
proc all_monomials(poly m1, poly m2, int k)
{
	list L;
	for(int i = 0; i <= k; i++)
	{
		L = L + list(m1^(k-i) * m2^i);
	}
	return(L);
}

// Find conditions for multiplicity index >= 4 for a given point
proc conditions(poly K, number cx, number cy, number cz)
{
	poly T = ax2 + by2 + cz2 + dxy + exz + fyz;
	poly m(1..2);
	if(cx != 0)
	{
		cy = cy/cx;
		cz = cz/cx;
		cx = 1;
		K = subst(K, x, 1, y, y + cy, z, z + cz);
		T = subst(T, x, 1, y, y + cy, z, z + cz);
		m(1) = y;
		m(2) = z;	
	}
	else
	{
		if(cy != 0)
		{
			cx = cx/cy;
			cz = cz/cy;
			cy = 1;
			K = subst(K, x, x + cx, y, 1, z, z + cz);
			T = subst(T, x, x + cx, y, 1, z, z + cz);
			m(1) = x;
			m(2) = z;
		}
		else
		{
			cx = cx/cz;
			cy = cy/cz;
			cz = 1;
			K = subst(K, x, x + cx, y, y + cy, z, 1);
			T = subst(T, x, x + cx, y, y + cy, z, 1);
			m(1) = x;
			m(2) = y;
		}
	}
	poly C(1) = subst(T, m(1), 0, m(2), 0);
	T = T - C(1);
	int k = 0;
	poly cf;
	cf = subst(subst(K, m(2), 0)/m(1), m(1), 0);
	if(cf != 0)
	{
		K = K / subst(subst(K, m(2), 0)/m(1), m(1), 0);
		k = 1;
	}
	cf = subst(subst(K, m(1), 0)/m(2), m(2), 0);
	if(cf != 0)
	{
		K = K / subst(subst(K, m(1), 0)/m(2), m(2), 0);
		k = 2;
	}
	if(k == 0)
	{
		return(list());
	}
	int l = 3 - k;
	int i, j;
	list L;
	for(i = 1; i < 4; i++)
	{
		L = all_monomials(m(1), m(2), i - 1);

		for(j = size(L); j >= 1; j--)
		{
			T = T - L[j]*K*subst(T/(L[j]*m(k)), m(l), 0, m(k), 0);
		}
		poly C(i+1) = subst(subst(T, m(k), 0)/(m(l)^i), m(l), 0);
		T = T - C(i+1) * m(l)^i;
	}
	ideal I;
	for(i = 1; i <= 4; i++){ I[i] = C(i); }
	return(I);
}

// Find all conics and pairs of points
proc find_conics(poly K, list P)
{
	int p = size(P) div 3;
	poly C = ax2 + by2 + cz2 + dxy + exz + fyz;
	list CON;
	list C1;
	list C2;
	ideal I;

	for(int i = 1; i <= p-1; i++)
	{
		for(int j = i+1; j <= p; j++)
		{
			I = conditions(K, P[3*i-2], P[3*i-1], P[3*i]), conditions(K, P[3*j-2], P[3*j-1], P[3*j]), C;
			I = std(I);
			if(size(I) == 6)
			{
				if(deg(I[6]) > 1)
				{
					CON = CON + list(factorize(I[6])[1][2]);
					C1 = C1 + list(i);
					C2 = C2 + list(j);
				}
			}
		}
	}
	return(list(CON, C1, C2))
}
