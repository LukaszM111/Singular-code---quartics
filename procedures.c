option(redSB);
option(noredefine);
option(noloadLib);
LIB "elim.lib";

// Make a point ideal from its coordinates
proc make_point(poly xc, poly yc, poly zc)
{
	ideal I = xc*y - yc*x, xc*z - zc*x, yc*z - zc*y;
	return(std(I));	
}

// Show point coordinates from a point ideal
proc point_coordinates(ideal P)
{
	ideal I;
	I = P, x - 1;
	I = std(I);
	if(size(I) != 3)
	{
		I = P, y - 1;
		I = std(I);
		if(size(I) != 3)
		{
			I = P, z - 1;
			I = std(I);
		}
	}
	poly xc = -subst(I[3]/(subst(I[3],x,1) - subst(I[3],x,0)), x, 0);
	poly yc = -subst(I[2]/(subst(I[2],y,1) - subst(I[2],y,0)), y, 0);
	poly zc = -subst(I[1]/(subst(I[1],z,1) - subst(I[1],z,0)), z, 0);
	return(list(xc, yc, zc));
}

// Dualize a point to a line
proc point_to_line(ideal P)
{
	list L = point_coordinates(P);
	return(x*L[1] + y*L[2] + z*L[3]);
}

// Dualize a list of points
proc dualize_points(list L)
{
	int l = size(L);
	list M;
		
	for(int i = 1; i <= l; i++)
	{
		M = M + list(point_to_line(L[i]));
	}
	
	return(M);
}

// Dualize a line to a point
proc line_to_point(poly L)
{
	poly xc = subst(L, x, 1, y, 0, z, 0);
	poly yc = subst(L, x, 0, y, 1, z, 0);
	poly zc = subst(L, x, 0, y, 0, z, 1);	
	ideal I = x*yc - y*xc, x*zc - z*xc, y*zc - z*yc;
	return(std(I));
}

// Dualize a list of lines
proc dualize_lines(list L)
{
	int l = size(L);
	list M;
	for(int i = 1; i <= l; i++)
	{
		M = M + list(line_to_point(L[i]));
	}	
	return(M);
}

// Does a given point lie on a given line?
proc is_on_line(ideal P, poly L)
{
	return(size(reduce(L, std(P))) == 0);
}

// Are these points the same?
proc same_points(ideal P1, ideal P2)
{
	return(size(reduce(P1, std(P2))) == 0);
}

// Intersection points and their multiplicities
proc line_arrangement(list L)
{
	list P;
	int new, m;
	int p = 0;
	int l = size(L);
	ideal I;
	for(int i = 1; i <= l-1; i++)
	{
		for(int j = i+1; j <= l; j++)
		{
			I = L[i], L[j];
			I = std(I);
			if(p > 0)
			{
				new = 1;
				for(int k = 1; k <= p; k++)
				{
					if(same_points(P[k], I))
					{
						new = 0;
						break;
					}
				}
				if(new)
				{
					P = P + list(I);
					p++;
				}
			}
			else
			{
				P = P + list(I);
				p++;
			}
		}
	}
	list M;
	for(int i = 1; i <= p; i++)
	{
		m = 0;
		for(int j = 1; j <= l; j++)
		{
			if(is_on_line(P[i], L[j]))
			{
				m++;
			}
		}
		M = M + list(m);
	}
	return(list(P, M));
}

// T-vector
proc t_vector(list L)
{
	list A = line_arrangement(L)[2];
	int m = 2;
	int k;
	int p = 0;
	int l = size(A);
	string s;
	while(p < l)
	{
		k = 0;
		for(int i = 1; i <= l; i++)
		{
			if(A[i] == m){ k++; }
		}
		p = p + k;
		if(k > 0)
		{
			s = "t_" + string(m) + " = " + string(k);
			s;
		}
		m++;
	}
}

// Roulleau P-operator
proc P_operator(list n, list L)
{
	list A = line_arrangement(L);
	int l = size(A[1]);
	int ln = size(n);
	list M;
	int i, j, ok;
	for(i = 1; i <= l; i++)
	{
		ok = 0;
		for(j = 1; j <= ln; j++)
		{
			if(A[2][i] == n[j])
			{
				ok = 1;
				break;
			}
		}
		if(ok)
		{
			M = M + list(A[1][i]);
		}
	}
	return(M);
}

// Intersection of all point ideals from the list
proc intersect_all(list L)
{
	int l = size(L);
	ideal I = L[1];
	for(int i = 2; i <= l; i++)
	{
		I = intersect(I, L[i]);
	}
	return(std(I));
}
