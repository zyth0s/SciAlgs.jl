
module m
  contains
  ! NOTE: integer size may not be C-compatible: use iso_c_binding
  ! NOTE: name will be mangled to __m_MOD_five: use bind(C, name="five")
  integer function five()
    five = 5
  end function five
end module m
