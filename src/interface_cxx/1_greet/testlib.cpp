
#include "jlcxx/jlcxx.hpp"

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("greet", []() {return "Hello";});
}

