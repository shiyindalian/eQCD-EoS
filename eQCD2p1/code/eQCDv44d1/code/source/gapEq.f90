subroutine gapEq(n, x, fvec)

  implicit none

  integer n
  real(8) x(n), fvec(n)                     !x(1)
  real(8) rho
  real(8) lam1,lam2,lam3,lam4,lam5,lam6,lam7,c,kappa

  common /gapPara/ lam1,lam2,lam3,lam4,lam5,lam6,lam7,c,kappa


  rho=x(1)

  fvec(1)=lam1 - c/(Sqrt(2.)*Sqrt(rho)) + lam2*(-kappa + rho) + &
          (lam3*(-kappa + rho)**2)/2. + (lam4*(-kappa + rho)**3)/6. + &
          (lam5*(-kappa + rho)**4)/24. + (lam6*(-kappa + rho)**5)/120. + &
          (lam7*(-kappa + rho)**6)/720.

end
