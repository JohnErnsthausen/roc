#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include "exceptions.hpp"

// ============
// Declarations
// ============

template <typename T>
class matrix
{
public:
    // === Constructors ===

    // Constructor to create an UNINITIALIZED matrix
    // WARNING: Good for performance, but make sure to never use any uninitialized elements!
    // First argument: number of rows
    // Second argument: number of columns
    matrix(const size_t &, const size_t &);

    // Copy constructor to create a new matrix with the same elements as an existing matrix
    matrix(const matrix &);

    // Move constructor to move the elements of an existing matrix to a new matrix
    matrix(matrix &&);

    // === Member functions ===

    // Overloaded operator = to assign the elements of one matrix to another matrix
    matrix &operator=(const matrix &);

    // Overloaded operator = to move the elements of one matrix to another matrix
    matrix &operator=(matrix &&);

    // Member function used to obtain (but not modify) the number of rows in the matrix
    size_t get_rows() const;

    // Member function used to obtain (but not modify) the number of columns in the matrix
    size_t get_cols() const;

    // Overloaded operator () used to access matrix elements WITHOUT range checking
    // The indices start from 0: m(0, 1) would be the element at row 1, column 2
    // First version: returns a reference, thus allows modification of the element
    T &operator()(const size_t &, const size_t &);

    // Overloaded operator () used to access matrix elements WITHOUT range checking
    // The indices start from 0: m(0, 1) would be the element at row 1, column 2
    // Second version: does not return a reference and declared as const,
    // does not allow modification of the element
    T operator()(const size_t &, const size_t &) const;

    // Static member function used to set the character width of the elements
    // when printing a matrix (will be used with std::setw)
    static void set_output_width(const int &);

    // === Friend functions ===

    // Overloaded binary operator << used to easily print out a matrix to a stream
    template <typename U>
    friend std::ostream &operator<<(std::ostream &, const matrix<U> &);

private:
    // The number of rows
    size_t rows{0};

    // The number of columns
    size_t cols{0};

    // A pointer to an array storing the elements of the matrix in flattened (1-dimensional) form
    T *elements{nullptr};

    // A smart pointer to manage the memory allocated for the elements
    std::unique_ptr<T[]> smart{nullptr};

    // The character width of the elements when printing a matrix (will be used with std::setw)
    static int output_width;
};

// Initialize output_width to have a default value of 5
template <typename T>
int matrix<T>::output_width{5};

// ==============
// Implementation
// ==============

// === Constructors ===

// Uninitialized constructor
template <typename T>
matrix<T>::matrix(const size_t &input_rows, const size_t &input_cols)
    : rows(input_rows), cols(input_cols)
{
    if (rows == 0 or cols == 0)
      throw sayMessage{"Dimentions of matrix must be positive.\n"};
    smart.reset(new T[rows * cols]);
    elements = smart.get();
}

// Copy constructor
template <typename T>
matrix<T>::matrix(const matrix<T> &m)
    : rows(m.rows), cols(m.cols)
{
    smart.reset(new T[rows * cols]);
    elements = smart.get();
    std::memcpy(elements, m.elements, rows * cols * sizeof(double));
}

// Move constructor
template <typename T>
matrix<T>::matrix(matrix<T> &&m)
    : rows(m.rows), cols(m.cols)
{
    smart = move(m.smart);
    elements = smart.get();
    m.rows = 0;
    m.cols = 0;
    m.elements = nullptr;
}

// === Member functions ===

// Copy assignment
template <typename T>
matrix<T> &matrix<T>::operator=(const matrix<T> &m)
{
    if (m.rows != rows or m.cols != cols)
      throw sayMessage{"Copy assignment requires equal matrix dimensions.\n"};
    if(elements == m.elements) return *this;
    rows = m.rows;
    cols = m.cols;
    smart.reset(new T[rows * cols]);
    elements = smart.get();
    std::memcpy(elements, m.elements, rows * cols * sizeof(double));
    return *this;
}

// Move assignment
template <typename T>
matrix<T> &matrix<T>::operator=(matrix<T> &&m)
{
    if (m.rows != rows or m.cols != cols)
      throw sayMessage{"Move assignment requires equal matrix dimensions.\n"};
    if(elements == m.elements) return *this;
    rows = m.rows;
    cols = m.cols;
    smart = move(m.smart);
    elements = smart.get();
    m.rows = 0;
    m.cols = 0;
    m.elements = nullptr;
    return *this;
}

template <typename T>
inline size_t matrix<T>::get_rows() const
{
    return rows;
}

template <typename T>
inline size_t matrix<T>::get_cols() const
{
    return cols;
}

// With reference
template <typename T>
inline T &matrix<T>::operator()(const size_t &row, const size_t &col)
{
    return elements[(cols * row) + col];
}

// No reference
template <typename T>
inline T matrix<T>::operator()(const size_t &row, const size_t &col) const
{
    return elements[(cols * row) + col];
}

// For pretty printing
template <typename T>
inline void matrix<T>::set_output_width(const int &w)
{
    output_width = w;
}

// === Friend functions ===

template <typename T>
std::ostream &operator<<(std::ostream &out, const matrix<T> &m)
{
    if (m.rows == 0 and m.cols == 0)
        out << "()\n";
    else
    {
        for (size_t i{0}; i < m.rows; i++)
        {
            out << "( ";
            for (size_t j{0}; j < m.cols; j++)
                out << std::setw(m.output_width) << m(i, j) << ' ';
            out << ")\n";
        }
        out << '\n';
    }
    return out;
}
