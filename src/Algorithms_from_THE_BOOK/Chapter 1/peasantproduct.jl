""" Computes the integer product c = a * b"""
function peasantproduct(a::T, b::T) where T <: Integer
  c = zero(T)
  while a > one(T)
    if isodd(a)
      c = c + b
    end
    a = a >> 1 # divide a by 2
    b = b << 1 # multiply b by 2
  end
  return c + b
end    

c = peasantproduct(10, 33)
