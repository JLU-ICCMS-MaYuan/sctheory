function helical ( x )

!*****************************************************************************80
!
!! helical() evaluates the Fletcher-Powell helical valley function.
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
!    real ( kind = rk ) X(3), the argument.
!
!  Output:
!
!    real ( kind = rk ) HELICAL, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ) fx3
  real ( kind = rk ) helical
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) x(3)

  if ( 0.0D+00 < x(1) ) then
    theta = atan2 ( x(2), x(1) ) / 2.0D+00 / r8_pi
  else if ( x(1) < 0.0D+00 ) then
    theta = 0.5D+00 + atan2 ( x(2), x(1) ) / 2.0D+00 / r8_pi
  else if ( x(1) == 0.0D+00 ) then
    theta = 0.25D+00
  end if

  fx1 = x(3) - 10.0D+00 * theta
  fx2 = sqrt ( x(1) * x(1) + x(2) * x(2) )
  fx3 = x(3)

  fx = 100.0D+00 * fx1 * fx1 &
     +             fx2 * fx2 &
     +             fx3 * fx3

  helical = fx

  return
end
