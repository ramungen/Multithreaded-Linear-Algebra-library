#pragma once

#include <functional>
#include <vector>
// a column vector

/*
TODO:

1. make the dot() function a friend instead
2. Implement transposition of vectors
3. remove print()
4. make destructor work
5. define copy constructor myself
*/
namespace la {

	class vector {
	public:
		vector(std::vector<double>& vect, size_t dim);
		vector(size_t dim);
		~vector();

		// elementwise operations
		vector operator*(const vector& rhs) const;
		vector operator*(double number) const;

		vector operator+(const vector& rhs) const;
		vector operator+(double number) const;

		vector operator-(const vector& rhs) const;
		vector operator-(double number) const;

		// dot product
		double dot(const vector& rhs) const;
		
		double& operator[](size_t index);
		double at(size_t index) const;

		size_t size() const;

		// for debugging
		void print() const;

	private:
		void perform_elementwise(double* vect, std::function<double(size_t)> action) const;
		vector(double* vect, size_t dim);

		double* vect_;
		size_t dim_;
	};
}