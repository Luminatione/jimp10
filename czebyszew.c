#include "makespl.h"
#include "matrix.h"
#include "piv_ge_solver.h"

#include <math.h>
#include <stdlib.h>

struct chebychew
{
	double a, b;
	int coefficientsAmount;
	double* coefficients;
};

//we r assuming x_j -> PI*(k-0.5)/N (?)
double getJthC(int j, int N, double y_j, double x_j)
{
	double c = 0;
	for (int k = 0; k < N; k++)
	{
		c += y_j * cos(j * x_j);
	}
	return c * 2 / N;
}
double getA(points_t* points)
{
	double a = points->x[0];
	for (int i = 1; i < points->n; i++)
	{
		a = min(a, points->x[i]);
	}
	return a;
}
double getB(points_t* points)
{
	double b = points->x[0];
	for (int i = 1; i < points->n; i++)
	{
		b = max(b, points->x[i]);
	}
	return b;
}
int getN(int pointsAmount)
{
	int	N = pointsAmount - 3 > 10 ? 10 : pointsAmount - 3;
	char* approximationBaseSizeStr = getenv("APPROX_BASE_SIZE");
	int approximationBaseSize = -1;
	if (approximationBaseSizeStr != NULL)
	{
		approximationBaseSize = atoi(approximationBaseSizeStr);
	}
	return approximationBaseSize > 0 ? approximationBaseSize : N;
}
double computeTx_n(double x, int n)
{
	double previous = 1;
	return cos(n * acos(x));
}
double castXtoU(double a, double b, double x)
{
	return (2 * x - a - b) / (b - a);
}

double chebyshewN(double u, int n, double kind)
{
	if (n < 0)
	{
		return 0;
	}
	if (n == 0)
	{
		return 1;
	}
	double previous = 1;
	double current = kind * u;
	for (int i = 1; i < n; ++i)
	{
		double temp = 2 * u * current - previous;
		previous = current;
		current = temp;
	}
	return current;
}
double T_n(double u, int n)
{
	return chebyshewN(u, n, 1);
}
double U_n(double u, int n)
{
	return chebyshewN(u, n, 2);
}
double f(double x, const struct chebychew polynomial)
{
	double u = castXtoU(polynomial.a, polynomial.b, x);
	double value = 0;
	for (int i = 0; i < polynomial.coefficientsAmount; i++)
	{
		value += polynomial.coefficients[i] * T_n(u, i);
	}
	return value;
}
double d1f(double x, const struct chebychew polynomial)
{
	double u = castXtoU(polynomial.a, polynomial.b, x);
	double value = 0;
	for (int i = 0; i < polynomial.coefficientsAmount; i++)
	{
		value += i * polynomial.coefficients[i] * U_n(u, i - 1);
	}
	return value;
}
double d2f(double x, const struct chebychew polynomial)
{
	double u = castXtoU(polynomial.a, polynomial.b, x);
	double value = polynomial.coefficients[0];
	for (int i = 1; i < polynomial.coefficientsAmount; i++)
	{
		value += i * i * polynomial.coefficients[i] * chebyshewN(u, i - 2, 3);//(i * T_n(u, i) - u * U_n(u, i - 1)) / (u * u - 1);
	}
	return value;
}
double d3f(double x, const struct chebychew polynomial)
{
	double u = castXtoU(polynomial.a, polynomial.b, x);
	double value = polynomial.coefficients[0];
	for (int i = 1; i < polynomial.coefficientsAmount; i++)
	{
		value += i * i * i * polynomial.coefficients[i] * chebyshewN(u, i - 3, 4);
		//((u * u - 1) * i * i * U_n(u, i - 1) - u * (i * T_n(u, i) - u *
		//U_n(u, i - 1)) - 2 * u * i * T_n(u, i) + (u * u + 1) * U_n(u, i - 1)) / (u * u - 1) / (u * u - 1);
	//(u * u + 1) comes from https://www.wolframalpha.com/input/?i2d=true&i=D%5Ba*n*Divide%5B%5C%2840%29n*T0%5C%2840%29x%5C%2841%29-x*U1%5C%2840%29x%5C%2841%29%5C%2841%29%2CPower%5Bx%2C2%5D-1%5D%2Cx%5D but why is it here?
	}
	return value;
}

struct chebychew createChebyshevPolynomial(double a, double b, int N, matrix_t* T)
{
	struct chebychew polynomial;
	polynomial.coefficientsAmount = T->rn;
	polynomial.coefficients = malloc(sizeof(*polynomial.coefficients) * T->rn);
	polynomial.a = a;
	polynomial.b = b;
	for (int i = 0; i < polynomial.coefficientsAmount; i++)
	{
		polynomial.coefficients[i] = get_entry_matrix(T, i, N);
	}
	return polynomial;
}

void fillSpline(spline_t* spl, struct chebychew polynomial)
{
	for (int i = 0; i < polynomial.coefficientsAmount; i++)
	{
		spl->x[i] = polynomial.a + i * (polynomial.b - polynomial.a) / (spl->n - 1);
		spl->f[i] = f(spl->x[i], polynomial);
		spl->f1[i] = d1f(spl->x[i], polynomial);
		spl->f2[i] = d2f(spl->x[i], polynomial);
		spl->f3[i] = d3f(spl->x[i], polynomial);
	}
}

matrix_t* pointsToSymetricMatrix(points_t* pts, double* u, int N)
{
	matrix_t* T = make_matrix(N, N + 1);
	//matrix is build like T_0(u_0) T_1(u_0)...
	//					   T_0(u_1) T_1(u_1)...
	// and so on
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < pts->n; k++)
			{
				add_to_entry_matrix(T, i, j, computeTx_n(u[k], j) * computeTx_n(u[k], i));
			}
		}
		for (int k = 0; k < pts->n; k++)
		{
			add_to_entry_matrix(T, i, N, pts->y[k] * computeTx_n(u[k], i));
		}
	}
	return T;
}
//we can only compute chebyshev polynomials for x in <-1, 1> so we have to keep our points in this interval
double* transformXToChebyshevInterval(points_t* points, double a, double b)
{
	double* u = malloc(points->n * sizeof(*u));
	for (int i = 0; i < points->n; i++)
	{
		u[i] = castXtoU(a, b, points->x[i]);
	}
}
void make_spl(points_t* pts, spline_t* spl)
{
	double a = getA(pts), b = getB(pts);
	int N = getN(pts->n) + 1;
	matrix_t* T = pointsToSymetricMatrix(pts, transformXToChebyshevInterval(pts, a, b), N);
#ifdef DEBUG
	write_matrix(T, stdout);
#endif
	if (piv_ge_solver(T) || alloc_spl(spl, N))
	{
		spl->n = 0;
		return;
	}
#ifdef DEBUG
	write_matrix(T, stdout);
#endif
	fillSpline(spl, createChebyshevPolynomial(a, b, N, T));
}
