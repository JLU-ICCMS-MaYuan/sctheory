function powell ( x )

!*****************************************************************************80
!
!! powell() evaluates the Powell quartic function.
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
!    real ( kind = rk ) X(4), the argument.
!
!  Output:
!
!    real ( kind = rk ) POWELL, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ) fx3
  real ( kind = rk ) fx4
  real ( kind = rk ) powell
  real ( kind = rk ) x(4)

  fx1 = x(1) + 10.0D+00 * x(2)
  fx2 = x(3) - x(4)
  fx3 = x(2) - 2.0D+00 * x(3)
  fx4 = x(1) - x(4)

  fx =            fx1 * fx1 &
     +  5.0D+00 * fx2 * fx2 &
     +            fx3 * fx3 * fx3 * fx3 &
     + 10.0D+00 * fx4 * fx4 * fx4 * fx4

  powell = fx

  return
end
