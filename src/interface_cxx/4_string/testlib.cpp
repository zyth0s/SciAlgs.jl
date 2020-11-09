#include <string>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"

int read_string(std::string s) {
  return s.length();
}
std::string echo_string(std::string s) {
  return s;
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("read_string", &read_string);
  //mod.method("echo_string", &echo_string);
}

