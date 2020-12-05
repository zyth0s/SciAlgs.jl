
module basic_math
contains

  subroutine mypow(val) ! mangled

    integer, intent(inout) :: val

    val = val ** 2

  end subroutine mypow

end module basic_math
