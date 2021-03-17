#pragma once

#include <cstring>
#include <iostream>
#include <vector>
#include "exceptions.hpp"

// ============
// Declarations
// ============

template <typename T>
class vectorf
{
 public:
  // === Constructors ===

  // Constructor to create an UNINITIALIZED vector
  // WARNING: Good for performance, but make sure to never use any uninitialized
  // elements! First argument: number of elements
  vectorf(const size_t &);

  // Constructor to create an initialized vectorf
  // Argument: a vector containing the elements of vectorf
  // Number of elements is inferred automatically
  vectorf(const std::vector<T> &);

  // Copy constructor to create a new vector with the same elements as an
  // existing vector
  vectorf(const vectorf &);

  // Move constructor to move the elements of an existing vector to a new vector
  vectorf(vectorf &&);

  // === Member functions ===

  // Overloaded operator = to assign the elements of one vector to another
  // vector
  vectorf &operator=(const vectorf &);

  // Overloaded operator = to move the elements of one vector to another vector
  vectorf &operator=(vectorf &&);

  // Member function used to obtain (but not modify) the size or length of the
  // vector
  size_t get_size() const;

  // Overloaded operator () used to access vector elements WITHOUT range
  // checking The indices start from 1: v(2) would be the element at row 2 First
  // version: returns a reference, thus allows modification of the element
  T &operator()(const size_t &);

  // Overloaded operator () used to access vector elements WITHOUT range
  // checking The indices start from 1: v(2) would be the element at row 2
  // Second version: does not return a reference and declared as const,
  // does not allow modification of the element
  T operator()(const size_t &) const;

  // Static member function used to set the character width of the elements
  // when printing a vector (will be used with std::setw)
  static void set_output_width(const int &);

  // === Friend functions ===

  // Overloaded binary operator << used to easily print out a matrix to a stream
  template <typename U>
  friend std::ostream &operator<<(std::ostream &, const vectorf<U> &);

 private:
  // The number of rows
  size_t size{0};

  // A pointer to an array storing the elements of the vector
  T *elements{nullptr};

  // A smart pointer to manage the memory allocated for the elements
  std::unique_ptr<T[]> smart{nullptr};

  // The character width of the elements when printing a matrix (will be used
  // with std::setw)
  static int output_width;
};

// Initialize output_width to have a default value of 5
template <typename T>
int vectorf<T>::output_width{5};

// ==============
// Implementation
// ==============

// === Constructors ===

// Uninitialized constructor
template <typename T>
vectorf<T>::vectorf(const size_t &input_size) : size(input_size)
{
  if (size == 0) throw sayMessage{"Size of vector must be positive.\n"};
  smart.reset(new T[ size ]);
  elements = smart.get();
}

template <typename T>
vectorf<T>::vectorf(const std::vector<T> &v)
{
  size = v.size();
  smart.reset(new T[ size ]);
  elements = smart.get();
  std::memcpy(elements, &v[ 0 ], size * sizeof(T));
}

// Copy constructor
template <typename T>
vectorf<T>::vectorf(const vectorf<T> &v) : size(v.size)
{
  smart.reset(new T[ size ]);
  elements = smart.get();
  std::memcpy(elements, v.elements, size * sizeof(T));
}

// Move constructor
template <typename T>
vectorf<T>::vectorf(vectorf<T> &&v) : size(v.size)
{
  smart = move(v.smart);
  elements = smart.get();
  v.size = 0;
  v.elements = nullptr;
}

// === Member functions ===

// Copy assignment
template <typename T>
vectorf<T> &vectorf<T>::operator=(const vectorf<T> &v)
{
  if (v.size != size)
    throw sayMessage{"Copy assignment requires equal vector size.\n"};
  if (elements == v.elements) return *this;
  size = v.size;
  smart.reset(new T[ size ]);
  elements = smart.get();
  std::memcpy(elements, v.elements, size * sizeof(T));
  return *this;
}

// Move assignment
template <typename T>
vectorf<T> &vectorf<T>::operator=(vectorf<T> &&v)
{
  if (v.size != size)
    throw sayMessage{"Move assignment requires equal vector size.\n"};
  if (elements == v.elements) return *this;
  size = v.size;
  smart = move(v.smart);
  elements = smart.get();
  v.size = 0;
  v.elements = nullptr;
  return *this;
}

template <typename T>
inline size_t vectorf<T>::get_size() const
{
  return size;
}

// With reference
template <typename T>
inline T &vectorf<T>::operator()(const size_t &row)
{
  return elements[ row - 1 ];
}

// No reference
template <typename T>
inline T vectorf<T>::operator()(const size_t &row) const
{
  return elements[ row - 1 ];
}

// For pretty printing
template <typename T>
inline void vectorf<T>::set_output_width(const int &w)
{
  output_width = w;
}

// === Friend functions ===

template <typename T>
std::ostream &operator<<(std::ostream &out, const vectorf<T> &v)
{
  if (v.size == 0)
    out << "()\n";
  else
  {
    out << "( ";
    for (size_t i{1}; i <= v.size; i++)
    {
      out << std::setw(v.output_width) << v(i) << ' ';
    }
    out << ")\n";
  }
  return out;
}
