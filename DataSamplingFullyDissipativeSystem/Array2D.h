//Copied and modified from HiFi
#pragma once

#include <iostream>

// #define DEBUG

namespace TEST {

template <typename T> class Array2D {

public:
  // #### constructors ####

  // default constructor

  Array2D();

  // constructor 1

  Array2D(int in_dim_0, int in_dim_1, int in_dim_2 = 1);

  // copy constructor

  Array2D(const Array2D<T> &in_Array2D);

  // assignment

  Array2D<T> &operator=(const Array2D<T> &in_Array2D);

  Array2D<T> &operator=(const T val);

  // destructor

  ~Array2D();

  // setup & free

  void setup(int in_dim_0, int in_dim_1, int in_dim_2 = 1);

  void free();

  Array2D<T> &cite(const Array2D<T> &in_Array2D);

  // access/set 1d

  T &operator()(int in_pos_0);

  // access/set 2d

  T &operator()(int in_pos_0, int in_pos_1);

  // access/set 3d

  T &operator()(int in_pos_0, int in_pos_1, int in_pos_2);

  // return dimension

  // T sum();

  int get_dim(int in_dim);

  // method to get maximum value of Array2D

  T get_max(void);

  // method to get minimum value of Array2D

  T get_min(void);

  // print

  void print(void);

  /*! Initialize Array2D to zero - Valid for numeric data types (int, float,
   * double) */
  void initialize_to_zero();

  /*! Initialize Array2D to given value */
  void initialize_to_value(const T val);

  bool owned;

protected:
  int dim_0;
  int dim_1;
  int dim_2;

public:
  T *_data;
};

// definitions

// #### constructors ####

// default constructor

template <typename T> Array2D<T>::Array2D() {
  dim_0 = 1;
  dim_1 = 1;
  dim_2 = 1;
  _data = new T[dim_0 * dim_1 * dim_2];
  owned = true;
}

// constructor 1

template <typename T>
Array2D<T>::Array2D(int in_dim_0, int in_dim_1, int in_dim_2) {
  dim_0 = in_dim_0;
  dim_1 = in_dim_1;
  dim_2 = in_dim_2;
  _data = new T[dim_0 * dim_1 * dim_2];
  owned = true;
}

// copy constructor

template <typename T> Array2D<T>::Array2D(const Array2D<T> &in_Array2D) {

  dim_0 = in_Array2D.dim_0;
  dim_1 = in_Array2D.dim_1;
  dim_2 = in_Array2D.dim_2;

  _data = new T[dim_0 * dim_1 * dim_2];

  for (int i = 0; i < dim_0 * dim_1 * dim_2; i++) {
    _data[i] = in_Array2D._data[i];
  }
  owned = true;
}

// assignment

template <typename T>
Array2D<T> &Array2D<T>::operator=(const Array2D<T> &in_Array2D) {
  int i;

  if (this == &in_Array2D) {
    return (*this);
  } else {
    delete[] _data;

    dim_0 = in_Array2D.dim_0;
    dim_1 = in_Array2D.dim_1;
    dim_2 = in_Array2D.dim_2;

    _data = new T[dim_0 * dim_1 * dim_2];
    // NOTE: THIS COPIES POINTERS; NOT VALUES
    for (i = 0; i < dim_0 * dim_1 * dim_2; i++) {
      _data[i] = in_Array2D._data[i];
    }

    owned = true;
    return (*this);
  }
}

template <typename T> Array2D<T> &Array2D<T>::operator=(const T val) {
  for (int i = 0; i < dim_0 * dim_1 * dim_2; i++)
    _data[i] = val;
  return (*this);
}

// destructor

template <typename T> Array2D<T>::~Array2D() {
  if (owned)
    delete[] _data;
  // do we need to deallocate gpu memory here as well?
}

// #### methods ####

// setup

template <typename T>
void Array2D<T>::setup(int in_dim_0, int in_dim_1, int in_dim_2) {
  delete[] _data;

  dim_0 = in_dim_0;
  dim_1 = in_dim_1;
  dim_2 = in_dim_2;

  _data = new T[dim_0 * dim_1 * dim_2];
  owned = true;
}

template <typename T> void Array2D<T>::free() {
  dim_0 = 0;
  dim_1 = 0;
  dim_2 = 0;
  if (owned) {
    delete[] _data;
    _data = nullptr;
  }
}

template <typename T>
Array2D<T> &Array2D<T>::cite(const Array2D<T> &in_Array2D) {
  dim_0 = in_Array2D.dim_0;
  dim_1 = in_Array2D.dim_1;
  dim_2 = in_Array2D.dim_2;
  _data = in_Array2D._data;
  owned = false;
  return (*this);
}

template <typename T> T &Array2D<T>::operator()(int in_pos_0) {
// if debug mode, check the bounds
#ifdef DEBUG
  if (in_pos_0 >= dim_0 || in_pos_0 < 0)
    throw("Array2D excell bound!");
#endif
  return _data[in_pos_0 * dim_1 * dim_2]; // column major with matrix indexing
}

template <typename T> T &Array2D<T>::operator()(int in_pos_0, int in_pos_1) {
// if debug mode, check the bounds
#ifdef DEBUG
  if (in_pos_0 >= dim_0 || in_pos_0 < 0 || in_pos_1 >= dim_1 || in_pos_1 < 0)
    throw("Array2D excell bound!");
#endif
  return _data[in_pos_0 * dim_1 * dim_2 +
               in_pos_1 * dim_2]; // row major with 0 indexing
}

template <typename T>
T &Array2D<T>::operator()(int in_pos_0, int in_pos_1, int in_pos_2) {
#ifdef DEBUG
  if (in_pos_0 >= dim_0 || in_pos_1 >= dim_1 || in_pos_2 >= dim_2 ||
      in_pos_0 < 0 || in_pos_1 < 0 || in_pos_2 < 0)
    throw("Array2D excell bound!");
#endif
  return _data[in_pos_0 * dim_1 * dim_2 + in_pos_1 * dim_2 +
               in_pos_2]; // column major with matrix indexing
}

// obtain dimension

template <typename T> int Array2D<T>::get_dim(int in_dim) {
  if (in_dim == 0) {
    return dim_0;
  } else if (in_dim == 1) {
    return dim_1;
  } else if (in_dim == 2) {
    return dim_2;
  } else {
    std::cout << "ERROR: Invalid dimension ... " << std::endl;
    return 0;
  }
}

/* template <typename T>
T Array2D<T>::sum(int s0_0, int s0_1, int s1_0=1, int s1_1=1)
{
    T   ss = 0;

    for(i=0; i<dim_0*dim_1*dim_2*dim_3; i++)
        ss += _data[i];

    return ss;
}  */

// method to calculate maximum value of Array2D
// Template specialization
template <typename T> T Array2D<T>::get_max(void) {
  int i;
  T max = -__DBL_MAX__;

  for (i = 0; i < dim_0 * dim_1 * dim_2; i++) {
    if (_data[i] > max)
      max = _data[i];
  }
  return max;
}

// method to calculate minimum value of Array2D
// Template specialization
template <typename T> T Array2D<T>::get_min(void) {
  int i;
  T min = __DBL_MAX__;

  for (i = 0; i < dim_0 * dim_1 * dim_2; i++) {
    if (_data[i] < min)
      min = _data[i];
  }
  return min;
}
// print

template <typename T> void Array2D<T>::print(void) {
  int i, j, k;
  bool threeD = (dim_2 == 1 ? false : true);
  for (i = 0; i < dim_0; ++i) {
    if (threeD)
      std::cout << std::endl << "ans(" << i << ",:,:) = " << std::endl;
    for (j = 0; j < dim_1; ++j) {
      for (k = 0; k < dim_2; ++k) {
        std::cout << " " << (*this)(i, j, k) << " ";
      }
      std::cout << std::endl;
    }
    if (threeD)
      std::cout << std::endl;
  }
}

// Initialize values to zero (for numeric data types)
template <typename T> void Array2D<T>::initialize_to_zero() {

  for (int i = 0; i < dim_0 * dim_1 * dim_2; i++) {
    _data[i] = 0;
  }
}

// Initialize Array2D to given value
template <typename T> void Array2D<T>::initialize_to_value(const T val) {
  for (int i = 0; i < dim_0 * dim_1 * dim_2; i++) {
    _data[i] = val;
  }
}
} // namespace MFast