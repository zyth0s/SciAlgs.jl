
You need to adapt build.sh to reflect the location of your
libcxxwrap, that can be inquired with

```julia
using CxxWrap
CxxWrap.prefix_path()
```
