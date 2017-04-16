#ifndef _ANPI_MATRIX_HPP
#define _ANPI_MATRIX_HPP

#include <cstddef>
#include <cstring>

namespace anpi
{
  /**
   * Row-major simple matrix class
   */
  template<typename T>
  class Matrix {
  private:
    /// All matrix data
    T* _data;
    /// Number of rows
    size_t _rows;
    /// Number of columns
    size_t _cols;

  public:
    /// Construct an empty matrix
    Matrix();
    
    /**
     * Construct a matrix rows x cols and initialize all
     * elements with the given value     
     */
    Matrix(const size_t rows,
	   const size_t cols,
	   const T initVal=T());

    /**
     * Construct a matrix rows x cols and initialize all
     * elements with the memory content at the given pointer
     */
    Matrix(const size_t rows,
	   const size_t cols,
	   const T *const initMem);

    /**
     * Copy constructor will do a deep copy on the given matrix
     */
    Matrix(const Matrix<T>& other);
    
    /** 
     * Release all memory
     */
    ~Matrix();

    /**
     * Deep copy another matrix
     */
    Matrix<T>& operator=(const Matrix<T>& other);
    
    /// Return pointer to a given row
    inline T* operator[](const size_t row) {return _data + row*_cols;}

    /// Return read-only pointer to a given row
    const T* operator[](const size_t row) const {return _data + row*_cols;}

    /// Return reference to the element at the r row and c column
    T& operator()(const size_t row,const size_t col) {
      return *(_data + (row*_cols + col));
    }

    /// Return const reference to the element at the r row and c column
    const T& operator()(const size_t row,const size_t col) const {
      return *(_data + (row*_cols + col));
    }

    /**
     * Allocate memory for the given number of rows and cols
     */
    void allocate(const size_t row,
		  const size_t col);

    /**
     * Fill all elements of the matrix with the given value
     */
    void fill(const T val);

    /**
     * Fill all elements of the matrix with the given memory block
     *
     * The user must ensure that the given memory block has enough elements
     */
    void fill(const T* mem);

    /**
     * Number of rows
     */
    inline size_t rows() const { return _rows; }

    /**
     * Number of columns
     */
    inline size_t cols() const { return _cols; }

    /**
     * Pointer to data block
     */
    inline T* data() { return _data; }

    /**
     * Pointer to data block
     */
    inline const T* data() const { return _data; }
    
  }; // class Matrix

  // --------------
  // Implementation
  // --------------

  template<typename T>
  Matrix<T>::Matrix() : _data(0),_rows(0),_cols(0) {}
 
  template<typename T>
  Matrix<T>::Matrix(const size_t r,
		    const size_t c,
		    const T initVal)
    : _data(0),_rows(0),_cols(0)  {
    allocate(r,c);
    fill(initVal);
  }

  /**
   * Construct a matrix rows x cols and initialize all
   * elements with the memory content at the given pointer
   */
  template<typename T>
  Matrix<T>::Matrix(const size_t r,
		    const size_t c,
		    const T *const initMem)
    : _data(0),_rows(0),_cols(0)  {
    allocate(r,c);
    fill(initMem);
  }
  

  template<typename T>
  Matrix<T>::Matrix(const Matrix<T>& other)
    : _data(0),_rows(0),_cols(0) {
    allocate(other._rows,other._cols);
    fill(other.data());
  }
  
  template<typename T>
  Matrix<T>::~Matrix() {
    delete[] _data;
    _data=0;
    _rows=0;
    _data=0;
  }
  
  template<typename T>
  Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
    allocate(other._rows,other._cols);
    fill(other.data());
    return *this;
  }

  template<typename T>
  void Matrix<T>::allocate(const size_t r,
		      const size_t c) {
    delete[] _data;
    _data = new T[r*c];
    _rows=r;
    _cols=c;
  }

  template<typename T>
  void Matrix<T>::fill(const T val) {
    T* end = _data + (_rows*_cols);
    for (T* ptr = _data;ptr!=end;++ptr) {
      *ptr = val;
    }
  }

  template<typename T>
  void Matrix<T>::fill(const T* mem) {
    std::memcpy(_data,mem,sizeof(T)*_rows*_cols);
  }

} // namespace ANPI


#endif
