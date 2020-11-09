
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"

// take Julia array as argument and return type
void test_array_set(jlcxx::ArrayRef<double> a, const int64_t i, const double v)
{
  a[i] = v;
}

// Constant 1D array
const double* const_vector()
{
  static double d[] = {1., 2., 3};
  return d;
}

// Constant Multidimensional
const double* const_matrix()
{
  static double d[2][3] = {{1., 2., 3}, {4., 5., 6.}};
  return &d[0][0];
}

// 1D array
double* mutable_vector()
{
  //double d[] = {1., 2., 3};
  double d[3];
  d[0] = 1.0;
  d[1] = 2.0;
  d[2] = 3.0;
  double* out = d;
  return out;
}

// Multidimensional
//double* mutable_matrix()
//{
//  double d[2][3] = {{1., 2., 3}, {4., 5., 6.}};
//  return &d[0][0];
//}


JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("test_array_set", &test_array_set);
  //mod.method("const_vector", []() { return jlcxx::make_const_array(const_vector(), 3); });
  //mod.method("const_matrix", []() { return jlcxx::make_const_array(const_matrix(), 3, 2); });
  //mod.method("mutable_vector", []() { return jlcxx::make_julia_array(mutable_vector(), 3); });
  //mod.method("mutable_matrix", []() { return jlcxx::make_julia_array(mutable_matrix(), 3, 2); });
}

