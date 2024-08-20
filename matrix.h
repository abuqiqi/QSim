#pragma once

#include <bits/stdc++.h>
using namespace std;

#define ll long long int
#define DTYPE complex<double>

template <typename T>
class Matrix {
private:
    void clear(); // Clear the matrix
public:
    ll row, col;
    T** data;
    static map<string, shared_ptr<Matrix>> MatrixDict; // A global matrix dictionary

    //
    // Constructors
    //
    Matrix(); // Default constructor
    Matrix(ll r, ll c); // Initialize a all-zero matrix
    Matrix(ll r, ll c, T** temp);
    Matrix(const Matrix& matrx); // Copy constructor
    Matrix(Matrix&& matrx); // Move constructor

    // 
    // Operations
    // 
    Matrix& operator=(const Matrix& matrx); // Copy assignment
    Matrix& operator=(Matrix&& matrx); // Move assignment

    Matrix operator+(const Matrix& matrx) const; // Matrix addition
    Matrix& operator+=(const Matrix& matrx); // Matrix addition
    Matrix operator*(const Matrix& matrx) const; // Matrix multiplication
    Matrix tensorProduct(const Matrix& matrx) const; // Tensor product

    void rotationX(double theta); // Rotation X gate matrix
    void rotationY(double theta); // Rotation Y gate matrix
    void rotationZ(double theta); // Rotation Z gate matrix
    void u1(double lambda); // U1 gate matrix

    void identity(ll r); // Set the matrix to be an identity matrix
    void zero(ll r, ll c); // Set the matrix to be a zero matrix

    bool isZero() const; // Check if the matrix is a zero matrix

    //
    // Utility functions
    //
    ll getSize() const; // Get the size of the matrix
    void print(string fname = "") const; // Print a Matrix
    static void initMatrixDict();
    static ll getMatrixDictSize(); // Get the number of elements in the matrix dictionary
    static void printMatrixDict(); // Print the matrix dictionary

    //
    // Destructor
    //
    ~Matrix();
};


// Static initialization helper class
class StaticInitializer {
public:
    StaticInitializer() {
        Matrix<DTYPE>::initMatrixDict();
    }
};
