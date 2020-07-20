subroutine nfdx(l,lb,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
!Calculate the derivatives of fermi distribution functions
  implicit none

  real(16) l,lb,nf0,nf1,nf2
  real(16) nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x
  real(16) T,mu 
  !temperature and chemical potential

  common /Tmu/ T,mu

  nfd0x=nf0 + lb*nf1 + l*nf2
  nfd1x=(-2*lb*nf1 - l*nf2 + 3*(nf0**2 + (lb*nf1 + l*nf2)**2 +                  &
       nf0*(-1 + 2*lb*nf1 + 2*l*nf2)))/T
  nfd2x=(18*nf0**3 + 27*nf0**2*(-1 + 2*lb*nf1 + 2*l*nf2) +                      &
    (-1 + 3*lb*nf1 + 3*l*nf2)*(6*lb**2*nf1**2 + 4*lb*nf1*(-1 + 3*l*nf2) +       &
       l*nf2*(-1 + 6*l*nf2)) +                                                  &
    9*nf0*(1 + 6*lb**2*nf1**2 + 2*l*nf2*(-2 + 3*l*nf2) +                        &
       lb*nf1*(-5 + 12*l*nf2)))/T**2
  nfd3x=(162*nf0**4 + 162*lb**4*nf1**4 + 324*nf0**3*(-1 + 2*lb*nf1 + 2*l*nf2) +   &
    216*lb**3*nf1**3*(-1 + 3*l*nf2) +                                           &
    8*lb*nf1*(-1 + 3*l*nf2)*(1 + 9*l*nf2*(-1 + 3*l*nf2)) +                      &
    l*nf2*(-1 + 3*l*nf2)*(1 + 18*l*nf2*(-1 + 3*l*nf2)) +                        &
    12*lb**2*nf1**2*(7 + 9*l*nf2*(-5 + 9*l*nf2)) +                              &
    27*nf0**2*(7 + 36*lb**2*nf1**2 + 4*l*nf2*(-7 + 9*l*nf2) +                   &
       8*lb*nf1*(-4 + 9*l*nf2)) +                                               &
    3*nf0*(-9 + 216*lb**3*nf1**3 + 36*lb**2*nf1**2*(-7 + 18*l*nf2) +            &
       8*lb*nf1*(11 + 27*l*nf2*(-2 + 3*l*nf2)) +                                &
       2*l*nf2*(29 + 18*l*nf2*(-5 + 6*l*nf2))))/T**3
  nfd4x=(1944*nf0**5 + 4860*nf0**4*(-1 + 2*lb*nf1 + 2*l*nf2) +                  &
    810*nf0**3*(5 + 24*lb**2*nf1**2 + 4*l*nf2*(-5 + 6*l*nf2) +                  &
       2*lb*nf1*(-11 + 24*l*nf2)) +                                             &
    45*nf0**2*(-27 + 432*lb**3*nf1**3 + 108*lb**2*nf1**2*(-5 + 12*l*nf2) +      &
       2*l*nf2*(83 + 216*l*nf2*(-1 + l*nf2)) +                                  &
       2*lb*nf1*(107 + 162*l*nf2*(-3 + 4*l*nf2))) +                             &
    (-1 + 3*lb*nf1 + 3*l*nf2)*(648*lb**4*nf1**4 +                               &
       864*lb**3*nf1**3*(-1 + 3*l*nf2) +                                        &
       l*nf2*(-1 + 6*l*nf2)*(1 + 36*l*nf2*(-1 + 3*l*nf2)) +                     &
       8*lb*nf1*(-1 + 3*l*nf2)*(2 + 27*l*nf2*(-1 + 4*l*nf2)) +                  &
       12*lb**2*nf1**2*(26 + 9*l*nf2*(-19 + 36*l*nf2))) +                       &
    3*nf0*(27 + 3240*lb**4*nf1**4 + 1620*lb**3*nf1**3*(-3 + 8*l*nf2) +          &
       60*lb**2*nf1**2*(41 + 108*l*nf2*(-2 + 3*l*nf2)) +                        &
       10*l*nf2*(-26 + 3*l*nf2*(43 + 108*l*nf2*(-1 + l*nf2))) +                 &
       5*lb*nf1*(-95 + 12*l*nf2*(61 + 27*l*nf2*(-7 + 8*l*nf2)))))/T**4
  nfd5x=(29160*nf0**6 + 29160*lb**6*nf1**6 +                                    &
    87480*nf0**5*(-1 + 2*lb*nf1 + 2*l*nf2) +                                    &
    58320*lb**5*nf1**5*(-1 + 3*l*nf2) +                                         &
    9*nf0*(-27 + 2*lb*nf1*(412 +                                                &
          45*lb*nf1*(-77 + 12*lb*nf1*(22 + 3*lb*nf1*(-11 + 6*lb*nf1)))) +       &
       374*l*nf2 + 900*l*lb*nf1*(-2 + 3*lb*nf1)*                                &
        (5 + 12*lb*nf1*(-2 + 3*lb*nf1))*nf2 +                                   &
       540*l**2*(-5 + 9*lb*nf1*(9 + 4*lb*nf1*(-9 + 10*lb*nf1)))*nf2**2 +        &
       2160*l**3*(5 + 6*lb*nf1*(-8 + 15*lb*nf1))*nf2**3 +                       &
       3240*l**4*(-7 + 30*lb*nf1)*nf2**4 + 19440*l**5*nf2**5) +                 &
    l*nf2*(-1 + 3*l*nf2)*(1 + 90*l*nf2*(1 - 6*l*nf2)**2*(-1 + 3*l*nf2)) +       &
    3240*lb**4*nf1**4*(13 + 27*l*nf2*(-3 + 5*l*nf2)) +                          &
    7290*nf0**4*(13 + 60*lb**2*nf1**2 + 4*l*nf2*(-13 + 15*l*nf2) +              &
       8*lb*nf1*(-7 + 15*l*nf2)) +                                              &
    12960*lb**3*nf1**3*(-1 + 3*l*nf2)*(1 + l*nf2*(-7 + 15*l*nf2)) +             &
    4*lb*nf1*(-1 + 3*l*nf2)*(8 +                                                &
       45*l*nf2*(-1 + 3*l*nf2)*(5 + 36*l*nf2*(-1 + 3*l*nf2))) +                 &
    1620*nf0**3*(-27 + 360*lb**3*nf1**3 +                                       &
       36*lb**2*nf1**2*(-13 + 30*l*nf2) +                                       &
       lb*nf1*(197 + 216*l*nf2*(-4 + 5*l*nf2)) +                                &
       4*l*nf2*(41 + 9*l*nf2*(-11 + 10*l*nf2))) +                               &
    27*nf0**2*(279 + 16200*lb**4*nf1**4 +                                       &
       12960*lb**3*nf1**3*(-2 + 5*l*nf2) +                                      &
       90*lb**2*nf1**2*(163 + 72*l*nf2*(-11 + 15*l*nf2)) +                      &
       20*lb*nf1*(-172 + 9*l*nf2*(133 + 360*l*nf2*(-1 + l*nf2))) +              &
       20*l*nf2*(-119 + 9*l*nf2*(53 + 18*l*nf2*(-6 + 5*l*nf2)))) +              &
    6*lb**2*nf1**2*(248 + 45*l*nf2*                                             &
        (-85 + 9*l*nf2*(59 + 12*l*nf2*(-14 + 15*l*nf2)))))/T**5

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nbdx(T,nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x)
!Calculate the derivatives of fermi distribution functions
  implicit none

  real(16) nb
  real(16) nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x
  real(16) T

  nbd0x=nb
  nbd1x=-((nb*(1 + nb))/T)
  nbd2x=(nb*(1 + nb)*(1 + 2*nb))/T**2
  nbd3x=-((nb*(1 + nb)*(1 + 6*nb*(1 + nb)))/T**3)
  nbd4x=(nb*(1 + nb)*(1 + 2*nb)*(1 + 12*nb*(1 + nb)))/T**4
  nbd5x=-((nb*(1 + nb)*(1 + 30*nb*(1 + nb)*(1 + 2*nb)**2))/T**5)

end
