
abstract type ISBN end

struct ISBN10 <: ISBN
   group
   publisher
   title
   checkdigit
end

struct ISBN13 <: ISBN
   prefix
   group
   publisher
   title
   checkdigit
end

"""
Construct 10-digit ISBN object from its parts
"""
function ISBN(group, publisher, title, checkdigit)

   ISBN10(group, publisher, title, checkdigit)
end

"""
Construct 13-digit ISBN object from its parts
"""
function ISBN(prefix, group, publisher, title, checkdigit)

   ISBN13(prefix, group, publisher, title, checkdigit)
end

"""
Construct any ISBN object from a string with its parts separated either by hyphens or whitespaces
"""
function ISBN(isbn::T) where T <: AbstractString
   if occursin('-', isbn)
      isbn = ISBN(parse.(Int,split(isbn,'-'))...)
   elseif occursin(' ', isbn)
      isbn = ISBN(parse.(Int,split(isbn,' '))...)
   end
   isbn
end

# Properties

"""
Full ISBN number from its parts
"""
function number(isbn::T) where T <: ISBN
   number = ""
   for field in fieldnames(typeof(isbn))
      number *= string(getfield(isbn, field))
   end
   parse(Int, number)
end

import Base.digits

digits(isbn::T) where T <: ISBN = digits(number(isbn))

"""
Calculate the check digit of a 10-digit ISBN.
"""
function checkdigit(isbn::ISBN10)
   # control digit = ∑ᵢ₌₁⁹ i aᵢ   (mod 11)
   suma = sum(collect(1:9) .* reverse(digits(isbn))[1:9])
   mod( suma, 11)
end

"""
Calculate the check digit of a 13-digit ISBN.
"""
function checkdigit(isbn::ISBN13)
   suma = sum(repeat([1,3],6) .* reverse(digits(isbn))[1:12])
   r = rem(suma, 10)
   r = r == 0 ? r : 10 - r
end

function isvalid(isbn::T) where T <: ISBN
   checkdigit(isbn) == isbn.checkdigit
end
