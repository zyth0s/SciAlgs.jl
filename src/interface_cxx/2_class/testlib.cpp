
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"

class Foo
{
  public:
    Foo(int i = 0) : m_value(i) {}
    int add(int i) const { return m_value + i; }
    Foo* thisptr() { return this; }
    Foo& thisref() { return *this; }
    Foo thiscopy() { return *this; }
    const int* valueptr() { return &m_value; }
  private:
    int m_value;
};

struct MyStruct { MyStruct() {} };

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.add_type<MyStruct>("MyStruct");
  mod.add_type<Foo>("Foo")
    .constructor<int>()
    .method("add", &Foo::add)
    .method("thisptr", &Foo::thisptr)
    .method("thisref", &Foo::thisref)
    .method("thiscopy", &Foo::thiscopy)
    .method("valueptr", &Foo::valueptr);
  mod.method("sumfoos", [] (const std::vector<Foo>& vec)
  {
    int result = 0;
    for (auto foo: vec)
    {
      result += *foo.valueptr();
    }
    return result;
  });
}

