#include <iostream>
#include <sstream>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <cstdlib>
#include <assert.h>
#include <complex>
#include <vector>
#include <ctime>
#include <cmath>
#define M_PI 3.14159265358979323846

using namespace std;

bool operator==(const complex<double>& lhs, const complex<double>& rhs){
	const double epsilon = 1e-7;
	return (abs(lhs - rhs) < epsilon);
}

template<class T>
class Matrix{
protected:
	size_t r, c;
	T **mtr;
public:

	Matrix(){
		r = c = 0;
		mtr = nullptr;
	}

	Matrix(size_t n, size_t m, T val){
		r = n, c = m;
		mtr = new T*[n];
		for (size_t i = 0; i < n; ++i)
			mtr[i] = new T[m];
		for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			mtr[i][j] = val;
	}

	Matrix(size_t n, size_t m, T** a){
		r = n, c = m;
		mtr = new T*[n];
		for (size_t i = 0; i < n; ++i)
			mtr[i] = new T[m];
		for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			mtr[i][j] = a[i][j];
	}

	Matrix(Matrix& a){
		r = a.r, c = a.c;
		mtr = new T*[r];
		for (size_t i = 0; i < r; ++i)
			mtr[i] = new T[c];
		for (size_t i = 0; i < r; ++i)
		for (size_t j = 0; j < c; ++j)
			mtr[i][j] = a[i][j];
	}

	Matrix& operator=(const Matrix& rhs){
		if (this == &rhs) return *this;

		for (size_t i = 0; i < r; ++i) delete[] mtr[i];
		delete[] mtr;
		r = 0, c = 0;

		if (rhs.mtr != nullptr){
			r = rhs.r;
			c = rhs.c;
			mtr = new T*[r];
			for (size_t i = 0; i < r; ++i)
				mtr[i] = new T[c];
			for (size_t i = 0; i < rhs.r; ++i)
			for (size_t j = 0; j < rhs.c; ++j)
				mtr[i][j] = rhs[i][j];
		}
		return *this;
	}

	T* operator [](size_t i){ return mtr[i]; }
	const T* const operator [](size_t i) const { return mtr[i]; }

	Matrix operator +(const Matrix& rhs){
		Matrix res = *this;
		assert(r == rhs.r && c == rhs.c);
		for (size_t i = 0; i < r; ++i)
		for (size_t j = 0; j < c; ++j)
			res[i][j] += rhs[i][j];
		return res;
	}

	Matrix& operator +=(const Matrix& rhs){
		assert(r == rhs.r && c == rhs.c);
		for (size_t i = 0; i < r; ++i)
		for (size_t j = 0; j < c; ++j)
			mtr[i][j] += rhs[i][j];
		return *this;
	}

	Matrix operator -(const Matrix& rhs){
		Matrix res = *this;
		assert(r == rhs.r && c == rhs.c);
		for (size_t i = 0; i < r; ++i)
		for (size_t j = 0; j < c; ++j)
			res[i][j] -= rhs[i][j];
		return res;
	}

	Matrix& operator -=(const Matrix& rhs){
		assert(r == rhs.r && c == rhs.c);
		for (size_t i = 0; i < r; ++i)
		for (size_t j = 0; j < c; ++j)
			mtr[i][j] -= rhs[i][j];
		return *this;
	}

	Matrix operator *(const Matrix& rhs){
		Matrix res(r, rhs.c, (T)0);
		assert(c == rhs.r);
		for (size_t i = 0; i < r; ++i)
		for (size_t k = 0; k < rhs.c; ++k)
		for (size_t j = 0; j < c; ++j)
			res[i][k] += mtr[i][j] * rhs[j][k];
		return res;
	}

	Matrix& operator *=(const Matrix& rhs){
		assert(c == rhs.r);
		*this = (*this) * rhs;
		return *this;
	}

	bool operator ==(const Matrix<T> &rhs){
		if (c != rhs.c || r != rhs.r) return false;
		for (size_t i = 0; i < r; ++i)
		for (size_t j = 0; j < c; ++j)
		if (mtr[i][j] != rhs[i][j]) return false;
		return true;
	}

	bool operator !=(const Matrix<T> &rhs){
		return !(*this == rhs);
	}

	Matrix& transpose(){
		if (r == c){
			for (size_t i = 0; i < r; ++i)
			for (size_t j = i + 1; j < r; ++j)
				swap((*this)[i][j], (*this)[j][i]);
		}
		else{
			Matrix res(c, r, (T)0);
			for (size_t i = 0; i < r; ++i)
			for (size_t j = 0; j < c; ++j)
				res[j][i] = (*this)[i][j];
			*this = res;
		}
		return *this;
	}

	size_t get_rows(){ return r; }
	size_t get_cols(){ return c; }

	~Matrix(){
		for (size_t i = 0; i < r; ++i) delete[] mtr[i];
		delete[] mtr;
	}

};


template<class T>
void reorderCells(Matrix<T>& a){
	vector<size_t> cells;
	size_t sz = 1;
	size_t n = a.get_cols();
	for (size_t i = 0; i < n; ++i)
	if (i + 1 < n && a[i][i + 1] == (T)1)
		++sz;
	else{
		cells.push_back(sz);
		sz = 1;
	}
	sort(cells.begin(), cells.end(), greater<size_t>());
	for (size_t i = 0; i < n - 1; ++i)
		a[i][i + 1] = 0;
	size_t k = 0;
	for (size_t i = 0; i < cells.size();++i){
		size_t csz = cells[i];
		while (--csz){
			a[k][k + 1] = 1;
			++k;
		}
		++k;
	}
}

template<class T>
void restoreUnits(Matrix<T>& A, Matrix<T>& J){
	assert(A.get_rows() == J.get_rows());
	assert(A.get_cols() == J.get_cols());
	assert(A.get_cols() == A.get_rows());
	size_t n = A.get_cols();
	Matrix<T> I(n, n, (T)0);                                //identity matrix
	for (size_t i = 0; i < n; ++i)
		I[i][i] = (T)1;

	for (size_t k = 0; k < n - 1; ++k){                     //iterate by units in J ( Jordan canonical form)
		if (J[k][k + 1] != (T)1) continue;

		Matrix<T> B = I;                                    //matrix B
		B[k + 1][k + 1] = A[k][k + 1];
		Matrix<T> RB = I;                                   //inverse matrix to B
		RB[k + 1][k + 1] = (T)1 / A[k][k + 1];

		A = (B * A) * RB;									//similiarity matrix transform - get one in A[k][k + 1]

		for (size_t i = 0; i < n; ++i){
			if (i == k + 1) continue;
			Matrix<T> C = I;
			C[k + 1][i] = A[k][i];
			Matrix<T> RC = I;
			RC[k + 1][i] = -A[k][i];

			A = (C * A) * RC;								//similiarity matrix transform - get zero in A[k][i]
		}
	}

	vector<pair<size_t, size_t> > cells;					//squares bounds [li,ri]
	size_t sz = 1;
	for (size_t k = 0; k < n; ++k){
		if (k != n - 1 && J[k][k + 1] == (T)1)
			++sz;
		else{
			cells.push_back(make_pair(k + 1 - sz, k));
			sz = 1;
		}
	}
	//reverse(cells.begin(), cells.end());

	//код ниже попытка перенести остатки из широкой части в узкую, но если на одной вертикали нужно переносить в более чем одном
	//квадрате , то тогда при умножении мы опять получаем на некоторых позициях остатки
	for (size_t ver = 0; ver < cells.size(); ++ver){
		for (size_t hor = 0; hor < cells.size(); ++hor){
			if (cells[hor].second - cells[hor].first <= cells[ver].second - cells[ver].first) continue;

			int narrow_side = cells[ver].second - cells[ver].first + 1;

			for (int j = cells[hor].second; j >(int)cells[hor].first; --j){
				int diag = min(narrow_side, j - (int)cells[hor].first);
				for (int i = cells[ver].second; i >(int)cells[ver].second - diag; --i){
					int nj = j - ((int)cells[ver].second - i);
					Matrix<T> D = I;
					D[i][nj - 1] = -A[i][nj];
					Matrix<T> RD = I;
					RD[i][nj - 1] = A[i][nj];

					A = (D * A) * RD;
				}
			}

			cout << "\nAfter:\n\n";
			output(&A);
		}
	}

}

complex<double> randomComplex(double eps){
	double r = eps * (double)rand() / (double)RAND_MAX;
	double x = (double)rand() * 2 * M_PI / (double)RAND_MAX;
	complex<double> c(r*cos(x), r*sin(x));
	return c;
}

template<class T>
void output(Matrix<complex<T> >* M){
	for (size_t i = 0; i < M->get_rows(); ++i)
	for (size_t j = 0; j < M->get_cols(); ++j){
		complex<T> to = (*M)[i][j];
		if (to == (complex<T>)0) to = (complex<T>)0;
		if (to == (complex<T>)1) to = (complex<T>)1;
		int ac = 10000;
		T fc = round(to.real() * ac) / ac, sc = round(to.imag() * ac) / ac;
		stringstream ss;
		ss << "(" << fc << "," << sc << ")";
		string s;
		ss >> s;
		while (s.size() < 17) s = ' ' + s;
		cout << s << (j == M->get_cols() - 1 ? '\n' : ' ');
	}
}

int main(){
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	srand(time(0));
	int n;
	double eps;
	cin >> n >> eps;
	complex<double> **a;
	a = new complex <double>*[n];
	for (int i = 0; i < n; ++i)
		a[i] = new complex<double>[n];
	for (int i = 0; i < n; ++i)
	for (int j = 0; j < n; ++j){
		int num;
		cin >> num;
		a[i][j] = num;
	}
	Matrix<complex<double> > J((size_t)n, (size_t)n, a);
	reorderCells(J);
	Matrix <complex<double> > E(n, n, (complex<double>)0);
	for (int i = 0; i < n; ++i)
	for (int j = 0; j < n; ++j)
		E[i][j] = randomComplex(eps);
	Matrix <complex<double> > A = E + J;
	cout << "before\n\n";
	output(&A);
	restoreUnits(A, J);
	cout << "\nafter\n\n";
	output(&A);
	return 0;
}