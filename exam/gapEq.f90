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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gapEq_s(n, x, fvec)

  implicit none

  integer n
  real(8) x(n), fvec(n)                     !x(1)
  real(8) rho,rho_s
  real(8) drho1V_s
  real(8) lam1,lam2,lam3,lam4,lam5,lam6,lam7,c_sigma_s,kappa


  common /gapPara_s/ lam1,lam2,lam3,lam4,lam5,lam6,lam7,c_sigma_s,kappa



  rho_s=x(1)
  rho=2.*rho_s

  drho1V_s=lam1 + lam2*(-kappa + rho) + (lam3*(-kappa + rho)**2)/2. + &
           (lam4*(-kappa + rho)**3)/6. + (lam5*(-kappa + rho)**4)/24. + &
           (lam6*(-kappa + rho)**5)/120. + (lam7*(-kappa + rho)**6)/720.

  fvec(1)=sqrt(2.)*sqrt(2.*rho_s)*drho1V_s-c_sigma_s

end

