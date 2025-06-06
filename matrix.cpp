#include "matrix.h"

//
// Constructors of Matrix
//

// Default constructor
template<typename T>
Matrix<T>::Matrix() {
    data = nullptr;
    row = 0;
    col = 0;
}

// Initialize a all-zero matrix
template<typename T>
Matrix<T>::Matrix(ll r, ll c) { 
    row = r;
    col = c;
    data = new T *[row];
    for (ll i = 0; i < row; i++) {
        data[i] = new T[col];
        for (ll j = 0; j < col; j++) {
            data[i][j] = 0;
        }
    }
}

// Initialize a matrix with a given 2D array
template<typename T>
Matrix<T>::Matrix(ll r, ll c, T **temp) {
    row = r;
    col = c;
    data = new T *[row];
    for (ll i = 0; i < row; i++) {
        data[i] = new T[col];
        for (ll j = 0; j < col; j++) {
            data[i][j] = *((T *) temp + i * col + j);
        }
    }
}

// Copy constructor
template<typename T>
Matrix<T>::Matrix(const Matrix<T> &matrx) {
    row = matrx.row;
    col = matrx.col;
    data = new T *[row];
    for (ll i = 0; i < row; i++) {
        data[i] = new T[col];
        memcpy(data[i], matrx.data[i], col * sizeof(T));
    }
}

// Move constructor
template<typename T>
Matrix<T>::Matrix(Matrix<T> &&matrx) {
    row = matrx.row;
    col = matrx.col;
    data = matrx.data;
    matrx.row = 0;
    matrx.col = 0;
    matrx.data = nullptr;
}

//
// Operations
//

template<typename T>
void Matrix<T>::clear() {
    if (data != nullptr) {
        for (ll i = 0; i < row; i++) {
            delete[] data[i];
        }
        delete[] data;
        data = nullptr;
    }
    row = 0;
    col = 0;
}

// Copy assignment
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& matrx) {
    if (this != &matrx) {
        clear();
        row = matrx.row;
        col = matrx.col;
        data = new T*[row];
        for (ll i = 0; i < row; i++) {
            data[i] = new T[col];
            memcpy(data[i], matrx.data[i], col * sizeof(T));
        }
    }
    return *this;
}

// Move assignment
template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>&& matrx) {
    if (this != &matrx) {
        clear();
        row = matrx.row;
        col = matrx.col;
        data = matrx.data;
        matrx.row = 0;
        matrx.col = 0;
        matrx.data = nullptr;
    }
    return *this;
}

// Matrix addition C = A + B
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& matrx) const {
    if (row != matrx.row || col != matrx.col) {
        cout << "[ERROR] Matrix +: row != matrx.row || col != matrx.col. " << endl;
        exit(1);
    }
    Matrix<T> temp(row, col);
    for (ll i = 0; i < row; i++) {
        for (ll j = 0; j < col; j++) {
            temp.data[i][j] = data[i][j] + matrx.data[i][j];
        }
    }
    return temp;
}

// Matrix addition A += B
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& matrx) {
    if (row != matrx.row || col != matrx.col) {
        cout << "[ERROR] Matrix +=: row != matrx.row || col != matrx.col. " << endl;
        exit(1);
    }
    for (ll i = 0; i < row; i++) {
        for (ll j = 0; j < col; j++) {
            data[i][j] += matrx.data[i][j];
        }
    }
    return *this;
}

// Matrix multiplication C = A * B
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& matrx) const {
    if (col != matrx.row) {
        cout << "[ERROR] Matrix *: col != matrx.row. " << endl;
        this->print();
        // matrx.print();
        exit(1);
    }
    Matrix<T> temp(row, matrx.col);
    for (ll i = 0; i < row; i++) {
        for (ll j = 0; j < matrx.col; j++) {
            for (ll k = col - 1; k >= 0; k--) {
                temp.data[i][j] += data[i][k] * matrx.data[k][j];
            }
        }
    }
    return temp;
}

template<typename T>
vector<T> Matrix<T>::operator*(const vector<T>& vec) const {
    if ((size_t)col != vec.size()) {
        cout << "[ERROR] Matrix * vector: col != vec.size(). " << endl;
        exit(1);
    }
    vector<T> temp(row, 0);
    for (ll i = 0; i < row; i++) {
        for (ll j = 0; j < col; j++) {
            temp[i] += data[i][j] * vec[j];
        }
    }
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& scalar) const {
    Matrix<T> temp(row, col);
    for (ll i = 0; i < row; i++) {
        for (ll j = 0; j < col; j++) {
            temp.data[i][j] = data[i][j] * scalar;
        }
    }
    return temp;
}

// Tensor product C = A tensorProduct B
template<typename T>
Matrix<T> Matrix<T>::tensorProduct(const Matrix<T>& matrx) const {
    Matrix<T> temp(row * matrx.row, col * matrx.col);
    for (ll ar = 0; ar < row; ar++){
        for (ll ac = 0; ac < col; ac++){
            if (data[ar][ac] == 0.0) continue;
            for (ll br = 0; br < matrx.row; br++){
                for (ll bc = 0; bc < matrx.col; bc++){
                    temp.data[ar * matrx.row + br][ac * matrx.col + bc] = data[ar][ac] * matrx.data[br][bc];
                }
            }
        }
    }
    return temp;
}

// Rotation X
template<typename T>
void Matrix<T>::rotationX(double theta) {
    T rx[2][2] = {{{cos(theta/2), 0}, {0, -sin(theta/2)}},
                  {{0, -sin(theta/2)}, {cos(theta/2), 0}}};
    clear();
    row = 2;
    col = 2;
    data = new T*[2];
    for (ll i = 0; i < 2; i++) {
        data[i] = new T[2];
        for (ll j = 0; j < 2; j++) {
            data[i][j] = rx[i][j];
        }
    }
}

// Rotation Y
template<typename T>
void Matrix<T>::rotationY(double theta) {
    T ry[2][2] = {{cos(theta/2), -sin(theta/2)},
                  {sin(theta/2), cos(theta/2)}};
    clear();
    row = 2;
    col = 2;
    data = new T*[2];
    for (ll i = 0; i < 2; i++) {
        data[i] = new T[2];
        for (ll j = 0; j < 2; j++) {
            data[i][j] = ry[i][j];
        }
    }
}

// Rotation Z
template<typename T>
void Matrix<T>::rotationZ(double theta) {
    T rz[2][2] = {{exp(-T(0, 1) * theta / 2.0), 0},
                  {0, exp(T(0, 1) * theta / 2.0)}};
    clear();
    row = 2;
    col = 2;
    data = new T*[2];
    for (ll i = 0; i < 2; i++) {
        data[i] = new T[2];
        for (ll j = 0; j < 2; j++) {
            data[i][j] = rz[i][j];
        }
    }
}

// U1 gate
template<typename T>
void Matrix<T>::u1(double lambda) {
    T u1[2][2] = {{1, 0},
                  {0, exp(T(0, 1) * lambda)}};
    clear();
    row = 2;
    col = 2;
    data = new T*[2];
    for (ll i = 0; i < 2; i++) {
        data[i] = new T[2];
        for (ll j = 0; j < 2; j++) {
            data[i][j] = u1[i][j];
        }
    }
}

// Set the matrix to be an identity matrix
template<typename T>
void Matrix<T>::identity(ll r) {
    clear();
    row = r;
    col = r;
    data = new T*[row];
    for (ll i = 0; i < row; i++) {
        data[i] = new T[col];
        for (ll j = 0; j < col; j++) {
            data[i][j] = (i == j) ? 1 : 0;
        }
    }
}

// Set the matrix to be a zero matrix
template<typename T>
void Matrix<T>::zero(ll r, ll c) {
    clear();
    row = r;
    col = c;
    data = new T*[row];
    for (ll i = 0; i < row; i++) {
        data[i] = new T[col];
        for (ll j = 0; j < col; j++) {
            data[i][j] = 0;
        }
    }
}

// Check if the matrix is a zero matrix
template<typename T>
bool Matrix<T>::isZero() const {
    for (ll i = 0; i < row; ++ i) {
        for (ll j = 0; j < col; ++ j) {
            if (data[i][j] != T{0}) {
                return false;
            }
        }
    }
    return true;
}

//
// Utility functions
//

// Get the size of the matrix
template<typename T>
ll Matrix<T>::getSize() const {
    // check if is a diagonal matrix
    if (row == col) {
        bool isDiagonal = true;
        for (ll i = 0; i < row; i++) {
            for (ll j = 0; j < col; j++) {
                if (i != j && data[i][j] != (DTYPE)0) {
                    isDiagonal = false;
                    break;
                }
            }
        }
        if (isDiagonal) {
            return row * sizeof(T);
        }
    }
    return row * col * sizeof(T);
}

// Print the matrix
template<typename T>
void Matrix<T>::print(string fname) const {
    ostream* out;
    ofstream file;

    if (fname.empty()) {
        out = &cout;
    } else {
        file.open(fname);
        if (!file.is_open()) {
            cerr << "[ERROR] Cannot open " << fname << endl;
            exit(1);
        }
        out = &file;
    }

    *out << "----- Matrix: [" << row << "] * [" << col << "] -----" << endl;
    out -> setf(std::ios::left);
    if (data == nullptr) {
        *out << "Matrix is Empty!" << endl;
    } else {
        for (ll i = 0; i < row; i++) {
            for (ll j = 0; j < col; j++) {
                out -> width(13);
                *out << fixed << setprecision(4) << data[i][j];
            }
            *out << endl;
        }
    }
    if (file.is_open()) {
        file.close();
    }
}

//
// Destructor
//
template<typename T>
Matrix<T>::~Matrix() {
    clear();
    // cout << "~Matrix Class Destruct!" << endl;
}

//
// Global notations
//

template<typename T>
map<string, shared_ptr<Matrix<T>>> Matrix<T>::MatrixDict; // A global matrix dictionary

template<typename T>
void Matrix<T>::initMatrixDict() {
    // T mark[1][1] = {{1}}; // placeholder
    // MatrixDict["-"] = make_shared<Matrix<T>>(1, 1, (T**)mark);

    // T zeros[2][1] = {{1}, {0}};
    // MatrixDict["ZEROS"] = make_shared<Matrix<T>>(2, 1, (T**)zeros);

    // T ones[2][1] = {{0}, {1}};
    // MatrixDict["ONES"] = make_shared<Matrix<T>>(2, 1, (T**)ones);

    // T plus[2][1] = {{1.0 / sqrt(2)}, {1.0 / sqrt(2)}};
    // MatrixDict["PLUS"] = make_shared<Matrix<T>>(2, 1, (T**)plus);

    // T minus[2][1] = {{1.0 / sqrt(2)}, {-1.0 / sqrt(2)}};
    // MatrixDict["MINUS"] = make_shared<Matrix<T>>(2, 1, (T**)minus);

    // -------------- Gates -----------------
    MatrixDict.clear();
    T ide[2][2] = {{1, 0},
                   {0, 1}};
    MatrixDict["I"] = make_shared<Matrix<T>>(2, 2, (T**)ide);
    MatrixDict["-"] = MatrixDict["I"];

    T h[2][2] = {{1.0 / sqrt(2), 1.0 / sqrt(2)},
                 {1.0 / sqrt(2), -1.0 / sqrt(2)}};
    MatrixDict["H"] = make_shared<Matrix<T>>(2, 2, (T**)h);

    T x[2][2] = {{0, 1},
                 {1, 0}};
    MatrixDict["X"] = make_shared<Matrix<T>>(2, 2, (T**)x);

    T y[2][2] = {{0, -T(0, 1)},
                 {T(0, 1), 0}};
    MatrixDict["Y"] = make_shared<Matrix<T>>(2, 2, (T**)y);

    T z[2][2] = {{1, 0},
                 {0, -1}};
    MatrixDict["Z"] = make_shared<Matrix<T>>(2, 2, (T**)z);

    T s[2][2] = {{1, 0},
                 {0, T(0, 1)}};
    MatrixDict["S"] = make_shared<Matrix<T>>(2, 2, (T**)s);

    T t[2][2] = {{1, 0},
                 {0, exp(T(0, 1) * T(acos(-1), 0) / T(4, 0))}};
    MatrixDict["T"] = make_shared<Matrix<T>>(2, 2, (T**)t);

    T swap[4][4] = {{1, 0, 0, 0},
                    {0, 0, 1, 0},
                    {0, 1, 0, 0},
                    {0, 0, 0, 1}};
    MatrixDict["SWAP"] = make_shared<Matrix<T>>(4, 4, (T**)swap);
}

template<typename T>
ll Matrix<T>::getMatrixDictSize() {
    ll dictsize = 0;
    for (auto it = MatrixDict.begin(); it != MatrixDict.end(); ++ it) {
        if (it->second != nullptr) {
            dictsize += it->second->getSize();
        }
    }
    return dictsize;
}

template<typename T>
void Matrix<T>::printMatrixDict() {
    // cout << "[INFO] The number of elements in MatrixDict: [" << getMatrixDictSize() << "]" << endl;
    for (auto it = MatrixDict.begin(); it != MatrixDict.end(); ++ it) {
        cout.width(13);
        cout << fixed << it->first << ": " << it->second << endl;
        if (it->second != nullptr) {
            it->second->print();
        }
    }
}

template class Matrix<DTYPE>;

// Create a static instance of StaticInitializer
static StaticInitializer staticInitializer;