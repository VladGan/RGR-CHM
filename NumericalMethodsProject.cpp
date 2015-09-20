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

bool operator==(const complex<double>& lhs, const complex<double>& rhs){  // precision of calculation - (1e-13 ~ 0)
	const double epsilon = 1e-12;
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

	T* operator [](size_t i){return mtr[i];}
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

	size_t get_rows(){return r;}
	size_t get_cols(){return c;}

	~Matrix(){	
		for (size_t i = 0; i < r; ++i) delete[] mtr[i];
		delete[] mtr;
	}

};


template<class T>
bool isJord(Matrix<T>& a){
	size_t n = a.get_cols();
	for (size_t i = 0; i < n; ++i){
		for (size_t j = 0; j < n; ++j){
			if (j != i + 1 && a[i][j] != (T)0) return false;
			else if (j == i + 1 && a[i][j] != (T)0 && a[i][j] != (T)1) return false;
		}
	}
	return true;
}

template<class T>
void reorderCells(Matrix<T>& A, Matrix<T>& J){
	size_t sz = 1;
	size_t n = J.get_cols();
	vector<pair<size_t, pair<vector<size_t>, T> > >  cells;					//size of cell, vector indices of rows in this cell, eigen value of cell
	Matrix<T> B(n, n, (T)0);                                                // sorted by jordan blocks Matrix A
	vector<size_t> rows;
	for (size_t i = 0; i < n; ++i){
		if (i + 1 < n && J[i][i + 1] == (T)1){
			rows.push_back(i);
			++sz;

		}
		else{
			rows.push_back(i);
			cells.push_back(make_pair(sz, make_pair(rows, J[i][i])));
			rows.clear();
			sz = 1;
		}
	}
	sort(cells.begin(), cells.end(), [](const pair<size_t, pair<vector<size_t>, T> > &lhs, const pair<size_t, pair<vector<size_t>, T> > &rhs){
		return lhs.first > rhs.first;
	});
	for (size_t i = 0; i < n - 1; ++i)
		J[i][i + 1] = 0;
	for (size_t i = 0; i < n; ++i)
		J[i][i] = 0;
	size_t k = 0;
	for (auto csz : cells){
		for (int i = k; i < k + csz.first; ++i)
			for (int j = 0; j < n; ++j)
				B[i][j] = A[csz.second.first[i - k]][j];
		//
		while (--csz.first){
			J[k][k] = csz.second.second;
			J[k][k + 1] = (T)1;
			++k;
		}
		J[k][k] = csz.second.second;
		++k;
	}

	A = B;
}

template<class T> 
void restoreUnits(Matrix<T>& A, Matrix<T>& J){
	assert(A.get_rows() == J.get_rows());
	assert(A.get_cols() == J.get_cols());
	assert(A.get_cols() == A.get_rows());
	size_t n = A.get_cols();								//matrix size
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
		
		cout << "\nПолучить единицу в позиции A["<<k<<"]["<<k + 1<<"]\n\n";
		output(&A);
		
		for (size_t i = 0; i < n; ++i){						//iterate by columns: A[k][i], for all (i != k + 1)
			if (i == k + 1) continue;
			
			if (i == k){
				T residual = A[k][k] - J[k][k];

				Matrix<T> C = I;
				C[k + 1][i] = residual;
				Matrix<T> RC = I;
				RC[k + 1][i] = -residual;

				A = (C * A) * RC;							//similiarity matrix transform - get jordan eigen value in A[k][k]
			}
			else{
				Matrix<T> C = I;
				C[k + 1][i] = A[k][i];
				Matrix<T> RC = I;
				RC[k + 1][i] = -A[k][i];

				A = (C * A) * RC;							//similiarity matrix transform - get zero in A[k][i]
			}
		}

		cout << "\nПеренести остатки с текущей строки\n\n";
		output(&A);
	}

	vector<pair<size_t, size_t> > cells;					//cells bounds [li,ri]
	size_t sz = 1;
	for (size_t k = 0; k < n; ++k){
		if (k != n - 1 && J[k][k + 1] == (T)1)
			++sz;
		else{
			cells.push_back(make_pair(k + 1 - sz, k));
			sz = 1;
		}
	}
	reverse(cells.begin(), cells.end());
	
	for (size_t hor = 0; hor < cells.size(); ++hor){                                
		for (size_t ver = 0; ver < cells.size(); ++ver){
			if (cells[ver].second - cells[ver].first == cells[hor].second - cells[hor].first){
				int left = cells[hor].first, right = cells[hor].second;
				
				for (int j = right; j > left; --j){
					for (int i = n - 1; i > cells[ver].second; --i){
						Matrix<T> D = I;
						D[i][j - 1] = -A[i][j];
						Matrix<T> RD = I;
						RD[i][j - 1] = A[i][j];

						A = (D * A) * RD;					  //similiarity matrix transform - transfer remains to the narrow side]
					}
				}
				
				if (cells[ver].second != n - 1){
					cout << "\nОбнулить текущие столбцы снизу ввверх(до нужного уровня) в пределах горионтального деления:\n\n";
					output(&A);
				}

				break;
			}
		}
	}

}

template<class T>
void output(Matrix<complex<T> >* M){
	for (size_t i = 0; i < M->get_rows(); ++i)
		for (size_t j = 0; j < M->get_cols(); ++j){
			complex<T> to = (*M)[i][j];
			if (to == (complex<T>)0) to = (complex<T>)0;
			if (to == (complex<T>)1) to = (complex<T>)1;
			int ac = 1000000;
			T fc = round(to.real() * ac) / ac, sc = round(to.imag() * ac) / ac;
			stringstream ss;
			ss << "(" << fc << "," << sc << ")";
			string s;
			ss >> s;
			while (s.size() < 23) s = ' ' + s;
			cout << s << (j == M->get_cols() - 1 ? '\n' : ' ');
		}
}

int main(){
	freopen("inputJ.txt", "r", stdin);
	freopen("outputJ.txt", "w", stdout);
	int n;
	cin >> n;
	complex<double> **a = new complex<double>*[n];
	complex<double> *diag = new complex<double>[n];
	for (int i = 0; i < n; ++i) 
		a[i] = new complex<double>[n];
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			cin >> a[i][j];
	for (int i = 0; i < n; ++i)
		cin >> diag[i];
	//end of reading
	Matrix<complex<double> > J((size_t)n, (size_t)n, complex<double>(0, 0));
	Matrix<complex<double> > A((size_t)n,(size_t)n, a);
	for (int i = 0; i < n; ++i) J[i][i] = diag[i]; 
	for (int i = 0; i < n - 1; ++i) 
		if (round(A[i][i + 1].real()) == 1.0)
			J[i][i + 1] = complex<double>(1.0, 0);
	reorderCells(A, J);
	std::cout << "Начальная матрица A = E + J\n\n";
	output(&A);
	restoreUnits(A, J);
	std::cout << "\nКонечное преобразование\n\n";
	output(&A);
	return 0;
}
