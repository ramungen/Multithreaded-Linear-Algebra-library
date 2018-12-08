#include "matrix.h"


namespace la {

	matrix id(size_t dim) {
		matrix identity(dim, dim);

		for (size_t row = 0; row < dim; ++row) {
			for (size_t col = 0; col < dim; ++col) {
				identity.at(row, col) = (size_t)(row == col);
			}
		}
		return identity;
	}

	void init_mat(double**& mat, size_t rows, size_t cols) {
		mat = new double*[rows];
		for (size_t row = 0; row < rows; ++row) {
			mat[row] = new double[cols];
		}
	}

	// matrix multiplication helper
	double product_ij(const matrix& lhs, const matrix& rhs, size_t i, size_t j) {
		double sum = 0;
		for (size_t k = 0; k < rhs.dims().first; ++k) {
			sum += (lhs.at(i, k) * rhs.at(k, j));
		}
		return sum;
	}

	matrix elementwise(const matrix& lhs, const matrix& rhs) {
		_ASSERT(lhs.dims() == rhs.dims());

		matrix newmat(lhs.dims().first, lhs.dims().second);

		for (size_t row = 0; row < lhs.dims().first; ++row) {
			for (size_t col = 0; col < lhs.dims().second; ++col) {
				newmat.at(row, col) = lhs.at(row, col) * rhs.at(row, col);
			}
		}
		return newmat;
	}

	inline std::pair<size_t, size_t> matrix::dims() const {
		return { rows_, columns_ };
	}

	matrix matrix::T() {
		matrix transpose(columns_, rows_);
		for (size_t row = 0; row < rows_; ++row) {
			for (size_t col = 0; col < columns_; col++) {
				transpose.at(col, row) = mat_[row][col];
			}
		}
		return transpose;
	}


	matrix::matrix(std::vector<std::vector<double>>& data, size_t rows, size_t cols)
		: mat_(nullptr), rows_(rows), columns_(cols) {

		init_mat(mat_, rows_, columns_);
		// copy the elements of the vector
		for (size_t row = 0; row < rows_; ++row) {
			for (size_t column = 0; column < columns_; ++column) {
				mat_[row][column] = data[row][column];
			}
		}
	}

	matrix::matrix(size_t rows, size_t cols) : columns_(cols), rows_(rows), mat_(nullptr) {
		init_mat(mat_, rows, cols);
	}

	matrix::matrix(const matrix& oth) noexcept :
	mat_(nullptr), 
	rows_(oth.rows_),
	columns_(oth.columns_)
	{
		init_mat(mat_, rows_, columns_);

		for (size_t row = 0; row < rows_; ++row) {
			for (size_t column = 0; column < columns_; ++column) {
				mat_[row][column] = oth.at(row, column);
			}
		}
	}

	matrix::matrix(matrix && oth) noexcept : 
		mat_(nullptr), 
		rows_(0), 
		columns_(0) 
	{
		std::swap(mat_, oth.mat_);
		std::swap(rows_, oth.rows_);
		std::swap(columns_, oth.columns_);
		
	}

	la::matrix::matrix(double** mat, size_t rows, size_t cols) :
		mat_(mat), rows_(rows), columns_(cols) {

	}

	const double& matrix::at(size_t row, size_t col) const {
		return mat_[row][col];
	}

	double& matrix::at(size_t row, size_t col) {
		return mat_[row][col];
	}

	matrix::~matrix() {
		clear();
	}

	std::vector<std::vector<double>> matrix::to_stl() const {
		std::vector<std::vector<double>> stl_matrix(rows_, std::vector<double>(columns_));
		for (size_t row = 0; row < rows_; ++row) {
			for (size_t col = 0; col < columns_; ++col) {
				stl_matrix[row][col] = mat_[row][col];
			}
		}
		return stl_matrix;
	}

	void matrix::print() const {
		for (size_t row = 0; row < rows_; ++row) {
			for (size_t col = 0; col < columns_; ++col) {
				std::cout << mat_[row][col] << ' ';
			}
			std::cout << '\n';
		}

	}


	void la::matrix::perform_elementwise(double** mat,
		std::function<double(size_t, size_t)> action) {
		for (size_t row = 0; row < rows_; ++row) {
			for (size_t col = 0; col < columns_; ++col) {
				mat[row][col] = action(row, col);
			}
		}
	}

	void matrix::clear() {
		for (size_t row = 0; row < rows_; ++row) {
			delete[] mat_[row];
		}
		delete[] mat_;
		mat_ = nullptr;
		rows_ = 0;
		columns_ = 0;
	}

	void matrix::swap(matrix& oth) {
		std::swap(mat_, oth.mat_);
		std::swap(rows_, oth.rows_);
		std::swap(columns_, oth.columns_);
	}

	matrix& matrix::operator=(matrix rhs) noexcept {
		// handling self assignment
		if (rhs.mat_ != this->mat_) {
			clear();
			swap(rhs);
			return *this;
		}
	}

	matrix operator*(matrix lhs, const matrix& rhs) {
		lhs *= rhs;
		return lhs;
	}

	vector matrix::operator*(const vector& rhs) {
		_ASSERT(columns_ == rhs.size());

		vector result(rows_);

		for (size_t row = 0; row < rows_; ++row) {
			double elem = 0;
			for (size_t col = 0; col < columns_; ++col) {
				elem += mat_[row][col] * rhs.at(col);
			}
			result[row] = elem;
		}
		return result;
	}

	matrix operator*(matrix lhs, double number) {
		lhs *= number;
		return lhs;
	}

	matrix operator*(double number, matrix rhs) {
		rhs *= number;
		return rhs;
	}

	matrix operator+(matrix lhs, const matrix& rhs) {
		lhs += rhs;
		return lhs;
	}

	matrix operator+(matrix lhs, const vector& rhs) {
		lhs += rhs;
		return lhs;
	}

	matrix operator+(matrix lhs, double number) {
		lhs += number;
		return lhs;
	}

	matrix operator+(double number, const matrix& rhs) {
		matrix newmat(rhs.rows_, rhs.columns_);
		for (size_t row = 0; row < rhs.dims().first; ++row) {
			for (size_t col = 0; col < rhs.dims().second; ++col) {
				newmat.at(row, col) = rhs.at(row, col) + number;
			}
		}
		return newmat;
	}

	matrix operator-(matrix lhs, const matrix& rhs) {
		lhs -= rhs;
		return lhs;
	}

	matrix operator-(matrix lhs, const vector& rhs) {
		lhs -= rhs;
		return lhs;
	}

	matrix operator-(matrix lhs, double number) {
		lhs -= number;
		return lhs;
	}

	matrix operator-(double number, const matrix& rhs) {
		matrix newmat(rhs.rows_, rhs.columns_);
		for (size_t row = 0; row < rhs.rows_; ++row) {
			for (size_t col = 0; col < rhs.columns_; ++col) {
				newmat.at(row, col) = number - rhs.at(row, col);
			}
		}
		return newmat;
	}

	matrix& matrix::operator*=(const matrix& rhs) {
		_ASSERT(this->dims().second == rhs.dims().first);

		double** newmat_ = nullptr;

		// for self-multiplication 
		size_t newrows = rows_;
		size_t newcols = rhs.columns_;

		init_mat(newmat_, rows_, rhs.columns_);

		for (size_t row = 0; row < rows_; ++row) {
			for (size_t col = 0; col < rhs.columns_; ++col) {
				newmat_[row][col] = product_ij(*this, rhs, row, col);
			}
		}

		clear();
		mat_ = newmat_;
		rows_ = newrows;
		columns_ = newcols;
		return *this;
	}

	matrix& matrix::operator*=(double number) {
		for (size_t row = 0; row < rows_; ++row) {
			for (size_t col = 0; col < columns_; ++col) {
				mat_[row][col] *= number;
			}
		}
		return *this;
	}

	matrix& matrix::operator+=(const matrix& rhs) {
		_ASSERT(this->dims() == rhs.dims());

		for (size_t row = 0; row < rhs.dims().first; ++row) {
			for (size_t col = 0; col < rhs.dims().second; ++col) {
				mat_[row][col] += rhs.at(row, col);
			}
		}
		return *this;
	}

	matrix& matrix::operator+=(const vector& rhs) {
		_ASSERT(rhs.size() == this->rows_);

		for (size_t row = 0; row < this->rows_; ++row) {
			for (size_t col = 0; col < this->columns_; ++col) {
				mat_[row][col] += rhs.at(row);
			}
		}
		return *this;
	}

	matrix& matrix::operator+=(double number) {
		
		for (size_t row = 0; row < rows_; ++row) {
			for (size_t col = 0; col < columns_; ++col) {
				 mat_[row][col] += number;
			}
		}
		return *this;
	}

	matrix& matrix::operator-=(const matrix& rhs) {
		_ASSERT(this->dims() == rhs.dims());

		for (size_t row = 0; row < rhs.dims().first; ++row) {
			for (size_t col = 0; col < rhs.dims().second; ++col) {
				mat_[row][col] -= rhs.at(row, col);
			}
		}
		return *this;
	}

	matrix& matrix::operator-=(const vector& rhs) {
		_ASSERT(rhs.size() == this->rows_);

		for (size_t row = 0; row < this->rows_; ++row) {
			for (size_t col = 0; col < this->columns_; ++col) {
				mat_[row][col] -= rhs.at(row);
			}
		}
		return *this;
	}

	matrix & matrix::operator-=(double number) {
		for (size_t row = 0; row < rows_; ++row) {
			for (size_t col = 0; col < columns_; ++col) {
				mat_[row][col] -= number;
			}
		}
		return *this;
	}

}