subroutine gapEq(n, x, fvec)

  implicit none

  integer n
  real(16) x(n), fvec(n)
  !x(1)
  real(16) rho
  real(16) lam1,lam2,lam3,lam4,lam5,lam6,lam7,c,kappa

  common /gapPara/ lam1,lam2,lam3,lam4,lam5,lam6,lam7,c,kappa

  rho=x(1)

  fvec(1)=lam1 - c/(Sqrt(2.Q+0)*Sqrt(rho)) + lam2*(-kappa + rho) + &
          (lam3*(-kappa + rho)**2)/2.Q+0 + (lam4*(-kappa + rho)**3)/6.Q+0 + &
          (lam5*(-kappa + rho)**4)/24.Q+0 + (lam6*(-kappa + rho)**5)/120.Q+0 + &
          (lam7*(-kappa + rho)**6)/720.Q+0

end
