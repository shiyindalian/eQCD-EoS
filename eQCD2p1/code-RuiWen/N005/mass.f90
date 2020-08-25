subroutine mass2(y,dmbo2drho)

  implicit none

  integer NMAX 
  !maximal number of differential equations
  parameter(NMAX=50)
  real(16) y(NMAX)
  integer N_str(5) 
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng
  !number of lam h Z ck
  real(16) lam10,lam20,lam30,lam40,lam50,lam60,lam70
  real(16) lam01,lam11,lam21,lam31,lam41,lam51
  real(16) lam02,lam12,lam22,lam32
  real(16) lam03,lam13
  !constant of Taylor expansion of potential V
  real(16) ck
  !anomaly breaking of axial U(1)A symmetry constant (KMT term)
  real(16) Sl,Ss
  real(16) rho,rho2
  !expansion points
  real(16) dmbo2drho(2,12)
  real(16) hd00s19,hd10s19,hd20s19,hd30s19,hd40s19,hd50s19
  real(16) hd01s19,hd11s19,hd21s19,hd31s19,hd02s19,hd12s19
  real(16) hd00p22,hd10p22,hd20p22,hd30p22,hd40p22,hd50p22
  real(16) hd01p22,hd11p22,hd21p22,hd31p22,hd02p22,hd12p22

  common /strucFun/ N_str

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lam10=y(1)
  lam20=y(2)
  lam30=y(3)
  lam40=y(4)
  lam50=y(5)
  lam60=y(6)
  lam70=y(7)
  lam01=y(8)
  lam11=y(9)
  lam21=y(10)
  lam31=y(11)
  lam41=y(12)
  lam51=y(13)
  lam02=y(14)
  lam12=y(15)
  lam22=y(16)
  lam32=y(17)
  lam03=y(18)
  lam13=y(19)
  Sl=y((Nv+1)+(Nh+2)+Nz+3)
  Ss=y((Nv+1)+(Nh+2)+Nz+4)
  ck=y((Nv+1)+(Nh+2)+Nz+5)

  rho=(Sl**2 + Ss**2)/2.Q+00
  rho2=(Sl**2 - 2*Ss**2)**2/24.Q+00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  hd00p22=lam10 + (ck*Sqrt(rho + Sqrt(6.Q+0)*Sqrt(rho2)))/Sqrt(3.Q+0) -                 &
    Sqrt(2.Q+0/3.Q+0)*lam01*Sqrt(rho2)
  hd10p22=lam20 + (ck/Sqrt(rho + Sqrt(6.Q+0)*Sqrt(rho2)) -                          &
       2*Sqrt(2.Q+0)*lam11*Sqrt(rho2))/(2.Q+0*Sqrt(3.Q+0))
  hd20p22=lam30 - ck/(4.Q+0*Sqrt(3.Q+0)*(rho + Sqrt(6.Q+0)*Sqrt(rho2))**1.5) -             &
    Sqrt(2.Q+0/3.Q+0)*lam21*Sqrt(rho2)
  hd30p22=lam40 + (Sqrt(3.Q+0)*ck)/(8.Q+0*(rho + Sqrt(6.Q+0)*Sqrt(rho2))**2.5) -           &
    Sqrt(2.Q+0/3.Q+0)*lam31*Sqrt(rho2)
  hd40p22=lam50 - (5*Sqrt(3.Q+0)*ck)/(16.Q+0*(rho + Sqrt(6.Q+0)*Sqrt(rho2))**3.5)
  hd50p22=(35*Sqrt(3.Q+0)*ck)/(32.Q+0*(rho + Sqrt(6.Q+0)*Sqrt(rho2))**4.5)
  hd01p22=lam11 + (-2*Sqrt(3.Q+0)*lam01 +                                           &
       (3*ck)/Sqrt(rho + Sqrt(6.Q+0)*Sqrt(rho2)) - 4*Sqrt(3.Q+0)*lam02*rho2)/           &
     (6.Q+0*Sqrt(2.Q+0)*Sqrt(rho2))
  hd11p22=lam21 + (-4*Sqrt(3.Q+0)*lam11 - (3*ck)/(rho + Sqrt(6.Q+0)*Sqrt(rho2))**1.5 -  &
       8*Sqrt(3.Q+0)*lam12*rho2)/(12.Q+0*Sqrt(2.Q+0)*Sqrt(rho2))
  hd21p22=lam31 + (-8*Sqrt(3.Q+0)*lam21 +                                           &
       (9*ck)/(rho + Sqrt(6.Q+0)*Sqrt(rho2))**2.5)/(24.*Sqrt(2.Q+0)*Sqrt(rho2))
  hd31p22=(-16*Sqrt(3.Q+0)*lam31 - (45*ck)/(rho + Sqrt(6.Q+0)*Sqrt(rho2))**3.5)/        &
    (48.*Sqrt(2.Q+0)*Sqrt(rho2))
  hd02p22=(24*lam12 + (-3*ck*(Sqrt(2.Q+0)*rho + 3*Sqrt(3.Q+0)*Sqrt(rho2)) +             &
         2*(Sqrt(6.Q+0)*rho + 6*Sqrt(rho2))*Sqrt(rho + Sqrt(6.Q+0)*Sqrt(rho2))*         &
          (lam01 - 4*lam02*rho2))/                                              &
       ((rho + Sqrt(6.Q+0)*Sqrt(rho2))**1.5*rho2**1.5))/24.Q+0
  hd12p22=(4*Sqrt(6.Q+0)*lam11 +                                                    &
      (3*ck*(Sqrt(2.Q+0)*rho + 5*Sqrt(3.Q+0)*Sqrt(rho2)))/                              &
       (rho + Sqrt(6.Q+0)*Sqrt(rho2))**2.5 - 16*Sqrt(6.Q+0)*lam12*rho2)/                &
    (48.Q+0*rho2**1.5)
  hd00s19=(9*lam10 + 12*lam20*rho + lam01*(4*rho - 5*Sqrt(6.Q+0)*Sqrt(rho2)) -      &
      3*Sqrt(3.Q+0)*ck*Sqrt(rho + Sqrt(6.Q+0)*Sqrt(rho2)) -                             &
      6*Sqrt(6.Q+0)*lam20*Sqrt(rho2) - 8*Sqrt(6.Q+0)*lam11*rho*Sqrt(rho2) +             &
      24*lam11*rho2 + lam02*(8*rho*rho2 - 4*Sqrt(6.Q+0)*rho2**1.5))/9.
  hd10s19=(4*lam01 + 21*lam20 + 4*lam11*rho + 12*lam30*rho -                    &
      (3*Sqrt(3.Q+0)*ck)/(2.*Sqrt(rho + Sqrt(6.Q+0)*Sqrt(rho2))) -                      &
      13*Sqrt(6.Q+0)*lam11*Sqrt(rho2) - 6*Sqrt(6.Q+0)*lam30*Sqrt(rho2) -                &
      8*Sqrt(6.Q+0)*lam21*rho*Sqrt(rho2) + 8*lam02*rho2 + 24*lam21*rho2 +           &
      8*lam12*rho*rho2 - 4*Sqrt(6.Q+0)*lam12*rho2**1.5)/9.Q+0
  hd20s19=(32*lam11 + 132*lam30 + 16*lam21*rho + 48*lam40*rho +                 &
      (3*Sqrt(3.Q+0)*ck)/(rho + Sqrt(6.Q+0)*Sqrt(rho2))**1.5 -                          &
      84*Sqrt(6.Q+0)*lam21*Sqrt(rho2) - 24*Sqrt(6.Q+0)*lam40*Sqrt(rho2) -               &
      32*Sqrt(6.Q+0)*lam31*rho*Sqrt(rho2) + 64*lam12*rho2 + 96*lam31*rho2)/36.Q+0
  hd30s19=(96*lam21 + 360*lam40 + 32*lam31*rho + 96*lam50*rho -                 &
      (9*Sqrt(3.Q+0)*ck)/(rho + Sqrt(6.Q+0)*Sqrt(rho2))**2.5 -                          &
      232*Sqrt(6.Q+0)*lam31*Sqrt(rho2) - 48*Sqrt(6.Q+0)*lam50*Sqrt(rho2))/72.
  hd40s19=(16*lam31)/9.Q+0 + (19*lam50)/3.Q+0 +                                       &
    (5*Sqrt(3.Q+0)*ck)/(16.Q+0*(rho + Sqrt(6.Q+0)*Sqrt(rho2))**3.5)
  hd50s19=(-35*Sqrt(3.Q+0)*ck)/(32.Q+0*(rho + Sqrt(6.Q+0)*Sqrt(rho2))**4.5)
  hd01s19=(-9*Sqrt(2.Q+0)*ck - 2*Sqrt(rho + Sqrt(6.Q+0)*Sqrt(rho2))*                    &
       (5*Sqrt(6.Q+0)*lam01 + 6*Sqrt(6.Q+0)*lam20 + 8*Sqrt(6.Q+0)*lam11*rho -               &
         66*lam11*Sqrt(rho2) - 24*lam02*rho*Sqrt(rho2) -                        &
         24*lam21*rho*Sqrt(rho2) + 22*Sqrt(6.Q+0)*lam02*rho2 +                      &
         12*Sqrt(6.Q+0)*lam21*rho2 + 16*Sqrt(6.Q+0)*lam12*rho*rho2 -                    &
         48*lam12*rho2**1.5))/                                                  &
    (36.*Sqrt(rho + Sqrt(6.Q+0)*Sqrt(rho2))*Sqrt(rho2))
  hd11s19=(96*lam02 + 96*lam12*rho + 96*lam31*rho +                             &
      lam21*(360 - (32*Sqrt(6.Q+0)*rho)/Sqrt(rho2)) -                               &
      (52*Sqrt(6.Q+0)*lam11)/Sqrt(rho2) - (24*Sqrt(6.Q+0)*lam30)/Sqrt(rho2) +           &
      (9*Sqrt(2.Q+0)*ck)/((rho + Sqrt(6.Q+0)*Sqrt(rho2))**1.5*Sqrt(rho2)) -             &
      152*Sqrt(6.Q+0)*lam12*Sqrt(rho2) - 48*Sqrt(6.Q+0)*lam31*Sqrt(rho2))/72.Q+0
  hd21s19=(8*lam12)/3.Q+0 + (lam31*(57 - (4*Sqrt(6.Q+0)*rho)/Sqrt(rho2)))/9.Q+0 +         &
    (-56*Sqrt(3.Q+0)*lam21 - 16*Sqrt(3.Q+0)*lam40 -                                     &
       (9*ck)/(rho + Sqrt(6.Q+0)*Sqrt(rho2))**2.5)/(24.Q+0*Sqrt(2.Q+0)*Sqrt(rho2))
  hd31s19=(-464*Sqrt(3.Q+0)*lam31 - 96*Sqrt(3.Q+0)*lam50 +                              &
      (135*ck)/(rho + Sqrt(6.Q+0)*Sqrt(rho2))**3.5)/(144.Q+0*Sqrt(2.Q+0)*Sqrt(rho2))
  hd02s19=(lam12*(57 - (8*Sqrt(6.Q+0)*rho)/Sqrt(rho2)) +                            &
      (10*Sqrt(6.Q+0)*lam01 + 12*Sqrt(6.Q+0)*lam20 + 16*Sqrt(6.Q+0)*lam11*rho +             &
         (9*Sqrt(2.Q+0)*ck)/Sqrt(rho + Sqrt(6.Q+0)*Sqrt(rho2)) +                        &
         (9*Sqrt(3.Q+0)*ck*Sqrt(rho2))/(rho + Sqrt(6.Q+0)*Sqrt(rho2))**1.5 -            &
         64*Sqrt(6.Q+0)*lam02*rho2 - 48*Sqrt(6.Q+0)*lam21*rho2)/(8.Q+0*rho2**1.5))/9.Q+0
  hd12s19=(52*Sqrt(6.Q+0)*lam11 + 24*Sqrt(6.Q+0)*lam30 + 32*Sqrt(6.Q+0)*lam21*rho -         &
      (9*Sqrt(2.Q+0)*ck)/(rho + Sqrt(6.Q+0)*Sqrt(rho2))**1.5 -                          &
      (27*Sqrt(3.Q+0)*ck*Sqrt(rho2))/(rho + Sqrt(6.Q+0)*Sqrt(rho2))**2.5 -              &
      256*Sqrt(6.Q+0)*lam12*rho2 - 96*Sqrt(6.Q+0)*lam31*rho2)/(144.Q+0*rho2**1.5)

  dmbo2drho(:,1)=(/hd00p22,hd00s19/)
  dmbo2drho(:,2)=(/hd10p22,hd10s19/)
  dmbo2drho(:,3)=(/hd20p22,hd20s19/)
  dmbo2drho(:,4)=(/hd30p22,hd30s19/)
  dmbo2drho(:,5)=(/hd40p22,hd40s19/)
  dmbo2drho(:,6)=(/hd50p22,hd50s19/)
  dmbo2drho(:,7)=(/hd01p22,hd01s19/)
  dmbo2drho(:,8)=(/hd11p22,hd11s19/)
  dmbo2drho(:,9)=(/hd21p22,hd21s19/)
  dmbo2drho(:,10)=(/hd31p22,hd31s19/)
  dmbo2drho(:,11)=(/hd02p22,hd02s19/)
  dmbo2drho(:,12)=(/hd12p22,hd12s19/)

end
