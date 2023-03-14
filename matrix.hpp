#include <algorithm>
#include <iostream>
#include <vector>

template <size_t N, size_t M, typename T = int64_t>
class Matrix {
 private:
  std::vector<std::vector<T>> matrix_;

 public:
  Matrix() {
    std::vector<T> fill(M, T());
    std::vector<std::vector<T>> matrix(N, fill);
    matrix_ = matrix;
  }
  
  Matrix(std::vector<std::vector<T>> matrix) { matrix_ = matrix; }
  
  Matrix(T elem) {
    std::vector<T> fill(M, elem);
    std::vector<std::vector<T>> matrix(N, fill);
    matrix_ = matrix;
  }
  
  Matrix& operator+=(const Matrix& other) {
    for (size_t index = 0; index < N; ++index) {
      for (size_t jin = 0; jin < M; ++jin) {
        matrix_[index][jin] += other.matrix_[index][jin];
      }
    }
    return *this;
  }
  
  Matrix& operator-=(const Matrix& other) {
    for (size_t index = 0; index < N; ++index) {
      for (size_t jin = 0; jin < M; ++jin) {
        matrix_[index][jin] -= other.matrix_[index][jin];
      }
    }
    return *this;
  }
  
  Matrix& operator*=(const T& multiplier) {
    for (size_t index = 0; index < N; ++index) {
      for (size_t jin = 0; jin < M; ++jin) {
        matrix_[index][jin] *= multiplier;
      }
    }
    return *this;
  }
  
  Matrix<M, N, T> Transposed() {
    Matrix<M, N, T> result;
    for (size_t index = 0; index < M; ++index) {
      for (size_t jin = 0; jin < N; ++jin) {
        result(index, jin) = matrix_[jin][index];
      }
    }
    return result;
  }
  
  T& operator()(const int64_t& ind, const int64_t& jin) { return matrix_[ind][jin]; }
  
  const T& operator()(const int64_t& ind, const int64_t& jin) const { return matrix_[ind][jin]; }
  
  bool operator==(const Matrix& other) const { return matrix_ == other.matrix_; }
};

template <size_t N, typename T>
class Matrix<N, N, T> {
 private:
  std::vector<std::vector<T>> matrix_;

 public:
  Matrix() {
    std::vector<T> fill(N, T());
    std::vector<std::vector<T>> matrix(N, fill);
    matrix_ = matrix;
  }
  
  Matrix(std::vector<std::vector<T>> matrix) { matrix_ = matrix; }
  Matrix(T elem) {
    std::vector<T> fill(N, elem);
    std::vector<std::vector<T>> matrix(N, fill);
    matrix_ = matrix;
  }
  
  Matrix& operator+=(const Matrix& other) {
    for (size_t index = 0; index < N; ++index) {
      for (size_t jin = 0; jin < N; ++jin) {
        matrix_[index][jin] += other.matrix_[index][jin];
      }
    }
    return *this;
  }
  
  Matrix& operator-=(const Matrix& other) {
    for (size_t index = 0; index < N; ++index) {
      for (size_t jin = 0; jin < N; ++jin) {
        matrix_[index][jin] -= other.matrix_[index][jin];
      }
    }
    return *this;
  }
  
  Matrix& operator*=(const T& multiplier) {
    for (size_t index = 0; index < N; ++index) {
      for (size_t jin = 0; jin < N; ++jin) {
        matrix_[index][jin] *= multiplier;
      }
    }
    return *this;
  }
  
  Matrix<N, N, T> Transposed() {
    Matrix<N, N, T> result;
    for (size_t index = 0; index < N; ++index) {
      for (size_t jin = 0; jin < N; ++jin) {
        result(index, jin) = matrix_[jin][index];
      }
    }
    return result;
  }
  
  T& operator()(const int64_t& ind, const int64_t& jin) { return matrix_[ind][jin]; }
  
  const T& operator()(const int64_t& ind, const int64_t& jin) const { return matrix_[ind][jin]; }
  
  bool operator==(const Matrix& other) const { return matrix_ == other.matrix_;}

  T Trace() {
    T result = T();
    for (size_t index = 0; index < N; ++index) {
      result += matrix_[index][index];
    }
    return result;
  }
};

template <size_t N, size_t M, typename T = int64_t>
Matrix<N, M, T> operator+(const Matrix<N, M, T>& left,
                          const Matrix<N, M, T>& right) {
  Matrix<N, M, T> result = left;
  result += right;
  return result;
}

template <size_t N, size_t M, typename T = int64_t>
Matrix<N, M, T> operator-(const Matrix<N, M, T>& left,
                          const Matrix<N, M, T>& right) {
  Matrix<N, M, T> result = left;
  result -= right;
  return result;
}

template <size_t N, size_t M, typename T = int64_t>
Matrix<N, M, T> operator*(const Matrix<N, M, T>& left, const T& right) {
  Matrix<N, M, T> result = left;
  result *= right;
  return result;
}

template <size_t N, size_t M, size_t K, typename T = int64_t>
Matrix<N, K, T> operator*(const Matrix<N, M, T>& left,
                          const Matrix<M, K, T>& right) {
  Matrix<N, K, T> result;
  for (size_t index = 0; index < N; ++index) {
    for (size_t jin = 0; jin < K; ++jin) {
      for (size_t kin = 0; kin < M; ++kin) {
        result(index, jin) += left(index, kin) * right(kin, jin);
      }
    }
  }
  return result;
}
