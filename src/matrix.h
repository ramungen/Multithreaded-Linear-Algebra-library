#pragma once

#include <vector>
#include <functional>
#include "vector.h"
#include <iostream>

/*
	TODO:
	1. implement row vector addition
	3. remove print()
*/
namespace la {

	class matrix {
	public:
		// transpose
		matrix T();
		std::pair<size_t, size_t> dims() const;

		matrix(std::vector<std::vector<double>>& data, size_t rows, size_t cols);
		matrix(size_t rows, size_t cols);
		// copy constructor
		matrix(const matrix& oth) noexcept;
		matrix(matrix&& oth) noexcept;
		// copy assignment
		matrix& operator=(matrix rhs) noexcept;

		// matrix multiplicaiton 
		friend matrix operator*(matrix lhs, const matrix& rhs);
		// vector multiplication
		vector operator*(const vector& rhs);

		// elementwise operations 
		matrix& operator*=(const matrix& rhs);
		matrix& operator*=(double number);
		
		matrix& operator+=(const matrix& rhs);
		matrix& operator+=(const vector& rhs);
		matrix& operator+=(double number);

		matrix& operator-=(const matrix& rhs);
		matrix& operator-=(const vector& rhs);
		matrix& operator-=(double number);

		// elementwise friend operations
		friend matrix operator*(matrix lhs, double number);
		friend matrix operator*(double number, matrix rhs);

		friend matrix operator+(matrix lhs, const matrix& rhs);
		friend matrix operator+(matrix lhs, const vector& rhs);
		friend matrix operator+(matrix lhs, double number);
		friend matrix operator+(double number, const matrix& rhs);

		friend matrix operator-(matrix lhs, const matrix& rhs);
		friend matrix operator-(matrix lhs, const vector& rhs);
		friend matrix operator-(matrix lhs, double number);
		friend matrix operator-(double number, const matrix& rhs);

		const double& at(size_t row, size_t col) const;
		double& at(size_t row, size_t col);
		friend matrix elementwise(const matrix& lhs, const matrix& rhs);
		~matrix();

		// convert to stl matrix
		std::vector<std::vector<double> > to_stl() const;

		// for debugging
		void print() const;
	private:
		matrix(double** mat, size_t rows, size_t cols);
		void perform_elementwise(double** mat, std::function<double(size_t, size_t)> action);

		// utility functions
		void clear();
		void swap(matrix& oth);

	// private members
	private:
		double** mat_;
		size_t rows_;
		size_t columns_;
	};


	// matrix functions
	matrix id(size_t dim);
}
