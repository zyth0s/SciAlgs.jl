
#include <valarray>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/const_array.hpp"

// take Julia array as argument and return type
void test_array_set(jlcxx::ArrayRef<double> a, const int64_t i, const double v)
{
  a[i] = v;
}

void test_matrix_set(jlcxx::ArrayRef<double,2> a, const int64_t i, const int64_t j, const double v)
{
  // Assumes ? x 3 matrix
  a[i*3+j] = v;
}

// Constant 1D array
const double* const_vector()
{
  // You need to declare static to prolong the life of d after exit,
  // otherwise you would get memory garbage.
  static double d[] = {1., 2., 3};
  return &d[0];
}

// Constant Multidimensional
const double* const_matrix()
{
  // You need to declare static to prolong the life of d after exit
  static double d[2][3] = {{1., 2., 3}, {4., 5., 6.}};
  return &d[0][0];
}

// 1D array
//double* mutable_vector()
//{
//  //double d[] = {1., 2., 3};
//  double d[3];
//  d[0] = 1.0;
//  d[1] = 2.0;
//  d[2] = 3.0;
//  double* out = d;
//  return out;
//}

// Multidimensional
//double* mutable_matrix()
//{
//  double d[2][3] = {{1., 2., 3}, {4., 5., 6.}};
//  return &d[0][0];
//}

void eval_rho(double& rho) {
    rho = 4.0;
}


JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("eval_rho", &eval_rho );
  mod.method("test_array_set", &test_array_set);
  mod.method("test_matrix_set", &test_matrix_set);
  mod.method("return_array", []() {
        // You need to declare static to prolong the life of a after exit
        static double a[2][3] = {{1., 2., 3}, {4., 5., 6.}};
        return jlcxx::make_julia_array(&a[0][0], 3, 2);
      });
  mod.method("const_vector",   []() { return jlcxx::make_const_array(const_vector(), 3); });
  mod.method("const_matrix",   []() { return jlcxx::make_const_array(const_matrix(), 3, 2); });
  //mod.method("mutable_vector", []() { return jlcxx::make_julia_array(const_vector(), 3); });
  //mod.method("mutable_matrix", []() { return jlcxx::make_julia_array(const_matrix(), 3, 2); });

  //mod.method("mutable_vector", []() { return jlcxx::make_julia_array(mutable_vector(), 3); });
  //mod.method("mutable_matrix", []() { return jlcxx::make_julia_array(mutable_matrix(), 3, 2); });
  mod.method("return_valarray", []() {
        double _a[] = {1., 2., 3, 4., 5., 6.};
        // You need to declare static to prolong the life of a after exit
        static std::valarray<double> a(_a, 6);
        return jlcxx::make_julia_array(&a[0], 3,2);
      });

}

