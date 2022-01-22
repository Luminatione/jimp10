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

double getA(points_t* points)
{
	double a = points->x[0];
	for (int i = 1; i < points->n; i++)
	{
		a = fmin(a, points->x[i]);
	}
	return a;
}
double getB(points_t* points)
{
	double b = points->x[0];
	for (int i = 1; i < points->n; i++)
	{
		b = fmax(b, points->x[i]);
	}
	return b;
}
int getN(int pointsAmount)
{
	int	N = pointsAmount > 5 ? 5 : pointsAmount;
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
double chebyshewN1d(double u, int n)
{
	if (n < 1)
	{
		return 0;
	}
	if (n == 1)
	{
		return 1;
	}
	double previous = 1;
	double current = 4 * u;
	for (int i = 2; i < n; ++i)
	{
		double temp = 2 * u * current - previous + T_n(u, i);
		previous = current;
		current = temp;
	}
	return current;
}
double chebyshewN2d(double u, int n)
{
	if (n < 2)
	{
		return 0;
	}
	if (n == 2)
	{
		return 4;
	}
	double previous = 4;
	double current = 24 * u;
	for (int i = 3; i < n; ++i)
	{
		double temp = 2 * u * current - previous + 4 * chebyshewN1d(u, i);
		previous = current;
		current = temp;
	}
	return current;
}
double chebyshewN3d(double u, int n)
{
	if (n < 3)
	{
		return 0;
	}
	if (n == 3)
	{
		return 24;
	}
	double previous = 24;
	double current = 192 * u;
	for (int i = 4; i < n; ++i)
	{
		double temp = 2 * u * current - previous + 6 * chebyshewN2d(u, i);
		previous = current;
		current = temp;
	}
	return current;
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
//https://en.wikipedia.org/wiki/Chebyshev_polynomials#Differentiation_and_integration
double d1f(double x, const struct chebychew polynomial)
{
	double u = castXtoU(polynomial.a, polynomial.b, x);
	double value = 0;
	for (int i = 0; i < polynomial.coefficientsAmount; i++)
	{
		value += polynomial.coefficients[i] * chebyshewN1d(u, i);
	}
	return value;
}
double d2f(double x, const struct chebychew polynomial)
{
	double u = castXtoU(polynomial.a, polynomial.b, x);
	double value = 0;
	if (u == 1.0 || u == -1.0)
	{
		for (int i = 0; i < polynomial.coefficientsAmount; i++)
		{
			value += polynomial.coefficients[i] * pow(u, i) * (pow(i, 4) - pow(i, 2)) / 3;
		}
		return value;
	}

	for (int i = 0; i < polynomial.coefficientsAmount; i++)
	{
		value += polynomial.coefficients[i] * chebyshewN2d(u, i);
	}
	return value;
}

double d3f(double x, const struct chebychew polynomial)
{
	double u = castXtoU(polynomial.a, polynomial.b, x);
	double value = 0;
	if (u == 1.0 || u == -1.0)
	{
		for (int i = 0; i < polynomial.coefficientsAmount; i++)
		{
			value += polynomial.coefficients[i] * pow(u, i + 1) * (pow(i, 6) - 5 * pow(i, 4) + 4 * pow(i, 2)) / 15;
		}
		return value;
	}
	for (int i = 0; i < polynomial.coefficientsAmount; i++)
	{
		value += polynomial.coefficients[i] * chebyshewN3d(u, i);
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
	//matrix is built like T_0(u_0) T_1(u_0)...
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
	return u;
}
void make_spl(points_t* pts, spline_t* spl)
{
	double a = getA(pts), b = getB(pts);
	int N = getN(pts->n);
	double* u = transformXToChebyshevInterval(pts, a, b);
	matrix_t* T = pointsToSymetricMatrix(pts, u, N);
	free(u);

	if (piv_ge_solver(T) || alloc_spl(spl, N))
	{
		spl->n = 0;
		return;
	}

	struct chebychew polynomial = createChebyshevPolynomial(a, b, N, T);
	fillSpline(spl, polynomial);

#ifdef _DEBUG
	FILE* debugOutput = fopen("dOutput", "w");
	for (int i = 0; i < spl->n; i++)
	{
		fprintf(debugOutput, "%g %g\n", spl->x[i], spl->f[i]);
	}
	fclose(debugOutput);
#endif

	free(polynomial.coefficients);
	freeMatrix(T);
	free(T);
}
