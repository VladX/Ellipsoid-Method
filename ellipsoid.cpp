#include <iostream>
#include <math.h>
#include <string.h>
using namespace std;

class Vector {
private:
	double * data;
	const size_t n;
public:
	inline Vector (size_t n) : n(n) {
		data = new double[n]();
	}
	
	inline Vector (size_t n, const double * d) : n(n) {
		data = new double[n];
		memcpy(data, d, sizeof(double) * n);
	}
	
	inline Vector (const Vector & v) : n(v.n) {
		data = new double[v.n];
		memcpy(data, v.data, sizeof(double) * v.n);
	}
	
	inline ~Vector () {
		delete[] data;
	}
	
	inline double operator[] (size_t i) const {
		return data[i];
	}
	
	inline double & operator[] (size_t i) {
		return data[i];
	}
	
	inline double operator* (const Vector & v) const {
		double res = 0;
		for (size_t i = 0; i < n; ++i)
			res += data[i] * v[i];
		return res;
	}
	
	inline void operator-= (const Vector & v) {
		for (size_t i = 0; i < n; ++i)
			data[i] -= v[i];
	}
	
	inline void operator/= (double s) {
		for (size_t i = 0; i < n; ++i)
			data[i] /= s;
	}
	
	inline void dump () const {
		for (size_t i = 0; i < n; ++i)
			cout << data[i] << ' ';
		cout << endl;
	}
};

class Matrix {
private:
	double * data;
	const size_t n;
public:
	inline Matrix (size_t n) : n(n) {
		data = new double[n*n]();
	}
	
	inline Matrix (const Matrix * m) : n(m->n) {
		data = new double[m->n * m->n];
		memcpy(data, m->data, sizeof(double) * m->n * m->n);
	}
	
	inline Matrix (const Matrix & m) : n(m.n) {
		data = new double[m.n * m.n];
		memcpy(data, m.data, sizeof(double) * m.n * m.n);
	}
	
	inline ~Matrix () {
		delete[] data;
	}
	
	inline size_t size () const {
		return n;
	}
	
	inline const double * operator[] (size_t i) const {
		return data + i * n;
	}
	
	inline double * operator[] (size_t i) {
		return data + i * n;
	}
	
	//%* Определитель *)
	inline double det () const {
		double det = 1;
		double ** a = new double *[n];
		for (size_t i = 0; i < n; ++i)
			a[i] = new double[n];
		const double EPS = 1E-9;
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < n; ++j)
				a[i][j] = data[i * n + j];
		for (size_t i = 0; i < n; ++i) {
			size_t k = i;
			for (size_t j=i+1; j<n; ++j)
				if (fabs(a[j][i]) > fabs(a[k][i]))
					k = j;
			if (fabs(a[k][i]) < EPS) {
				det = 0;
				break;
			}
			swap(a[i], a[k]);
			if (i != k)
				det = -det;
			det *= a[i][i];
			for (size_t j = i + 1; j < n; ++j)
				a[i][j] /= a[i][i];
			for (size_t j = 0; j < n; ++j)
				if (j != i && fabs(a[j][i]) > EPS)
					for (size_t k=i+1; k<n; ++k)
						a[j][k] -= a[i][k] * a[j][i];
		}
		for (size_t i = 0; i < n; ++i)
			delete[] a[i];
		delete[] a;
		return det;
	}
	
	inline void operator*= (double s) {
		for (size_t i = 0; i < n*n; ++i)
			data[i] *= s;
	}
	
	inline void operator-= (const Matrix & m) {
		for (size_t i = 0; i < n*n; ++i)
			data[i] -= m.data[i];
	}
};

//%* Единичная матрица *)
class IdentityMatrix : public Matrix {
public:
	IdentityMatrix (size_t n) : Matrix(n) {
		for (size_t i = 0; i < n; ++i)
			this->operator[](i)[i] = 1;
	}
};

//%* Объем эллипсоида с точностью до константы *)
double volume (const Matrix & B) {
	return sqrt(B.det());
}

//%* Радиус сферы, сод. многогранник *)
double ComputeInitialRadius (const Matrix & A, const Vector & b) {
	const size_t n = A.size();
	double R = ceil(sqrt(n)) * ceil(1.0/pow(2.0, n*n));
	for (size_t i = 0; i < n; ++i)
		R *= fabs(b[i]) + 1.0;
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
			R *= fabs(A[i][j]) + 1.0;
	return R;
}

int main () {
	size_t n; //%* Размер матрицы A (n*n) *)
	cin >> n;
	IdentityMatrix I(n);
	Matrix A(n);
	Vector b(n);
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
			cin >> A[i][j]; //%* Считываем матрицу A *)
	for (size_t i = 0; i < n; ++i)
		cin >> b[i]; //%* Считываем вектор b *)
	double R = ComputeInitialRadius(A, b); //%* Считаем начальный радиус сферы, целиком сод. многогранник *)
	Vector x(n);
	Matrix B(I); //%* Начальный эллипсоид *)
	B *= R*R;
	const double L = 1e-2;
	while (volume(B) > L) {
		size_t i = 0;
		for (i = 0; i < n; ++i) { //%* Находим первое неравенство, которое не выполняется *)
			double s = 0;
			for (size_t j = 0; j < n; ++j)
				s += A[i][j] * x[j];
			if (s > b[i])
				break;
		}
		if (i == n) { //%* Все неравенства выполнены, решение найдено *)
			cout << "%*Решение:*)\n";
			for (i = 0; i < n; ++i)
				cout << x[i] << ' ';
			cout << endl;
			return 0;
		}
		Vector a(n, A[i]);
		Vector Bka(n);
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < n; ++j)
				Bka[i] += B[i][j] * a[j];
		double aTBka = a * Bka;
		{ //%* Пересчитываем x (центр эллипсоида) *)
			Vector tmp(Bka);
			tmp /= (n + 1) * sqrt(aTBka);
			x -= tmp;
		}
		{ //%* Пересчитываем B (матрица, задающая эллипсоид) *)
			Matrix tmp(n);
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < n; ++j)
					tmp[i][j] = Bka[i] * Bka[j];
			tmp *= 2.0 / ((n + 1.0) * aTBka);
			B -= tmp;
			B *= ((double) n*n) / (n*n - 1.0);
		}
	}
	cout << "%*Решений не найдено*)" << endl;
	return 0;
}