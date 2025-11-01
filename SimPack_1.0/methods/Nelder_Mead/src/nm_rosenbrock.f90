function rosenbrock ( x )

!*****************************************************************************80
!
!! rosenbrock() evaluates the Rosenbrock parabolic value function.
!
!  Discussion:
!
!    Thanks to Vivek Rao for pointing out a discrepancy in the local
!    dimensioning of X.
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
!    real ( kind = rk ) X(2), the argument.
!
!  Output:
!
!    real ( kind = rk ) ROSENBROCK, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ) rosenbrock
  real ( kind = rk ) x(2)

  fx1 = x(2) - x(1) * x(1)
  fx2 = 1.0D+00 - x(1)

  fx = 100.0D+00 * fx1 * fx1 &
     +             fx2 * fx2

  rosenbrock = fx

  return
end
