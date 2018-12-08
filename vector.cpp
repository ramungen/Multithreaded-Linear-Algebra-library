#include "vector.h"
#include <iostream>

namespace la {

	vector::vector(std::vector<double>& vect, size_t dim) :
		dim_(dim), vect_(new double[dim]) {
		for (size_t i = 0; i < dim; ++i) {
			vect_[i] = vect[i];
		}
	}

	vector::vector(size_t dim) : dim_(dim) {
		vect_ = new double[dim_];
	}

	vector::~vector() {
		//delete[] vect_;
	}

	vector::vector(double* vect, size_t dim) :
		dim_(dim), vect_(vect) {

	}
	vector vector::operator-(const vector& rhs) const {
		if (rhs.dim_ != this->dim_) {
			throw std::invalid_argument("vector dimensions do not match");
		}
		double* new_vect_ptr = new double[dim_];
		double* raw_rhs = rhs.vect_;
		perform_elementwise(new_vect_ptr, [&, this](size_t i) {return vect_[i] - raw_rhs[i]; });

		return vector(new_vect_ptr, dim_);
	}
	vector vector::operator-(double number) const {
		return operator+(-number);
	}
	double vector::dot(const vector & rhs) const {
		if (rhs.dim_ != this->dim_) {
			throw std::invalid_argument("dimensions not compatible");
		}
		double sum = 0;
		for (size_t i = 0; i < rhs.dim_; i++) {
			sum += (rhs.at(i) + vect_[i]);
		}
		return sum;
	}
	double& vector::operator[](size_t index) {
		return vect_[index];
	}
	double vector::at(size_t index) const {
		return vect_[index];
	}
	void vector::perform_elementwise(double* vect, std::function<double(size_t)> action) const {
		for (size_t i = 0; i < dim_; ++i) {
			vect[i] = action(i);
		}
	}

	vector vector::operator*(const vector& rhs) const {
		if (rhs.dim_ != this->dim_) {
			throw std::invalid_argument("dimensions not compatible");
		}
		double* new_vect_ptr = new double[dim_];
		double* raw_rhs = rhs.vect_;
		perform_elementwise(new_vect_ptr, [&, this](size_t i) {return raw_rhs[i] * vect_[i]; });

		return vector(new_vect_ptr, dim_);
	}

	vector vector::operator*(double number) const {
		double* new_vect_ptr = new double[dim_];
		perform_elementwise(new_vect_ptr, [&, this](size_t i) {return number * vect_[i]; });

		return vector(new_vect_ptr, dim_);
	}
	vector vector::operator+(const vector & rhs) const {
		if (rhs.dim_ != this->dim_) {
			throw std::invalid_argument("dimensions not compatible");
		}
		double* new_vect_ptr = new double[dim_];
		double* raw_rhs = rhs.vect_;
		perform_elementwise(new_vect_ptr, [&, this](size_t i) {return vect_[i] + raw_rhs[i]; });

		return vector(new_vect_ptr, dim_);
	}
	vector vector::operator+(double number) const {
		double* new_vect_ptr = new double[dim_];
		perform_elementwise(new_vect_ptr, [&, this](size_t i) {return number + vect_[i]; });

		return vector(new_vect_ptr, dim_);
	}

	size_t vector::size() const {
		return dim_;
	}

	void vector::print() const {
		for (size_t i = 0; i < dim_; ++i) {
			std::cout << vect_[i] << '\n';
		}
	}

}
