function quartic ( x )

!*****************************************************************************80
!
!! quartic() evaluates a function defined by a sum of fourth powers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Input:
!
!    real ( kind = rk ) X(10), the argument.
!
!  Output:
!
!    real ( kind = rk ) QUARTIC, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) quartic
  real ( kind = rk ) x(10)

  quartic = sum ( x(1:10) ** 4 )

  return
end
