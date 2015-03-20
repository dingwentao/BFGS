#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define n 1000

int ludcmp(double a[n][n], int order[n]);
void solvlu(double a[n][n], double b[n], double x[n], int order[n]);
void swaprows(double a[n][n], int i, int j);
static int pivot(double a[n][n], int order[n], int jcol);
double gradient(double x[n], int i);
void output(double x[n]);
double diff(double x0[n], double x1[n]);

int ludcmp(double a[n][n], int order[n])
{
	int i, j, k, nm1;
	int flag = 1;    /* changes sign with each row interchange */
	double sum, diag;
	
	for (i=0; i<n; i++) order[i] = i;

	if (pivot(a,order,0)) flag = -flag;
	diag = 1.0/a[0][0];
	for (i=1; i<n; i++) a[0][i] *= diag;

	nm1 = n-1;
	for (j = 1; j<nm1; j++)
	{
		for (i=j; i<n; i++)
		{
			sum = 0.0;
			for (k=0; k<j; k++)
			  sum += a[i][k]*a[k][j];
			a[i][j] -= sum;
		}
		if (pivot(a,order,j)) flag = -flag;
		diag = 1.0/a[j][j];
		for (k=j+1; k<n; k++)
		{
			sum = 0.0;
			for (i=0; i<j; i++)
			  sum += a[j][i]*a[i][k];
			a[j][k] = (a[j][k]-sum)*diag;
		}
	}

	sum = 0.0;
	for (k=0; k<nm1; k++)
	  sum += a[nm1][k]*a[k][nm1];
	a[nm1][nm1] -= sum;
	return flag;
}

int pivot(double a[n][n], int order[n], int jcol)
{
	int i, ipvt;
	double big, anext;
	ipvt = jcol;
	big = fabs(a[ipvt][ipvt]);
	for (i = ipvt+1; i<n; i++)
	{
		anext = fabs(a[i][jcol]);
		if (anext>big)
		{
			big = anext;
			ipvt = i;
		}
	}

	if (ipvt == jcol) return 0;
	swaprows(a, jcol, ipvt);
	i = order[jcol];
	order[jcol] = order[ipvt];
	order[ipvt] = i;
	return 1;
}

void swaprows(double a[n][n], int i, int j)
{
	for (int k = 0; k < n; k++)
	{
		double temp;
		temp = a[i][k];
		a[i][k] = a[j][k];
		a[j][k] = temp;
	}
	return;
}

void solvlu(double a[n][n], double b[n], double x[n], int order[n])
{
	int i,j;
	double sum;
	for (i=0; i<n; i++)
	{
		j = order[i];
		x[i] = b[j];
	}

	x[0] /= a[0][0];
	for (i=1; i<n; i++)
	{
		sum = 0.0;
		for (j=0; j<i; j++)
		  sum += a[i][j]*x[j];
		x[i] = (x[i]-sum)/a[i][i];
	}

	for (i=n-2; i>=0; i--)
	{
		sum = 0.0;
		for (j=i+1; j<n; j++)
		  sum += a[i][j]*x[j];
		x[i] -= sum;
	}
	return;
}

double diff(double x0[n], double x1[n])
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	  sum += fabs(x0[i] - x1[i]);
	return sum;
}

void output(double x[n])
{
	int i;
	for (i=0; i<10; i++)
	  printf ("%f ", x[i]);
	for (i=n-10; i<n; i++)
	  printf ("%f ", x[i]);
	printf ("\n");
	return;
}

double gradient(double x[n], int i)
{
	double m = 0.0;
	for (int j=0; j<n; j++)
	  m += x[j]*x[j];
	return (x[i]*m-1);
}

int main()
{
	int i, j, k, l;
	int itr = 1;
	double B[n][n], T[n][n];
	double x0[n], x1[n], s[n], neggrad[n], y[n], t[n];
	int order[n];
	double difference;

	for (i=0; i<n; i++)
	  for (j=0; j<n; j++)
		B[i][j] = 0.0;
	for (i=0; i<n; i++)
	  B[i][i] = 1.0;
	for (i=0; i<n; i++)
	  x1[i] = 0.0;

	clock_t start, end;
	start = clock();
	do{
		for (i=0; i<n; i++)
		  x0[i] = x1[i];
		
		for (i=0; i<n; i++)
		  neggrad[i] = -gradient(x0, i);
		
		for (i=0; i<n; i++)
		  for (j=0; j<n; j++)
			T[i][j] = B[i][j];
		
		ludcmp(B, order);
		solvlu(B, neggrad, s, order);
		
		for (i=0; i<n; i++)
		  for (j=0; j<n; j++)
			B[i][j] = T[i][j];

		for (i=0; i<n; i++)
		  x1[i] = x0[i] + s[i];
		
		for (i=0; i<n; i++)
		  y[i] = gradient(x1, i) - gradient(x0,i);

		double c1, c2;
		c1 = 0.0;
		c2 = 0.0;
		for (i=0; i<n; i++)
		{
			double sum = 0.0;
			for (j=0; j<n; j++)
			  sum += B[i][j]*s[j];
			t[i] = sum;
			c1 += y[i]*s[i];
			c2 += s[i]*t[i];
		}

		for (i=0; i<n; i++)
		  for (j=0; j<n; j++)
		  {
			  double sum = 0.0;
			  for (k=0; k<n; k++)
				sum += t[i]*s[k]*B[k][j];
			  T[i][j] = sum;
		  }
		
		for (i=0; i<n; i++)
		  for (j=0; j<n; j++)
			B[i][j] += y[i]*y[j]/c1 - T[i][j]/c2;
		
		difference = diff(x0, x1);
		printf ("iteration = %d, diff = %f\n", itr++, difference);
	} while (difference >= 1e-8);
	end = clock();
	printf ("Time = %f\n", (double)(end-start)/CLOCKS_PER_SEC);
	output(x1);
	return 0;
}
