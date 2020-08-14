
module basic_math
contains

  subroutine mypow(val) bind(C, name='mypow') ! to avoid name mangling

    use iso_c_binding ! C-compatible types
    integer(c_long_long), intent(inout) :: val

    val = val ** 2

  end subroutine mypow

end module basic_math
