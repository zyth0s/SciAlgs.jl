
module basic_math
contains

  subroutine myfun(val, n) bind(C, name='myfun') ! to avoid name mangling

    use iso_c_binding ! C-compatible types
    integer(c_long_long),                 intent(in)    :: n
    integer(c_long_long), dimension(n,n), intent(inout) :: val

    val(1,2) = 4
    print *, sum(val)
    print *, n, " x ", n

  end subroutine myfun

end module basic_math
