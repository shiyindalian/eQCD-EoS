subroutine phypoint(Nflow,yflow,rho0,mPion,mSigma,mf,Vall)
!calculate the physical point rho0

  implicit none
  integer Nflow
  real(8) yflow(Nflow)
  real(8) rho0,rho
  real(8) mPion2,mPion,mSigma2,mSigma,mf
  real(8) Vall !total effective potential
  integer N_str(5) !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng
  real(8) lam0,lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(8) h
  real(8) Zphi,Zpsi,ZA
  real(8) c,kappa
  real(8) g
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  integer nn
  parameter(nn=1)
  real(8) x(nn)
  logical check
  external gapEq,gapEq_s
  real(8) rho_s,ms0
  real(8) drho1V,drho1V_s,mass_s
  real(8) c_bar_sigma_s
  real(8) lam1_s,lam2_s,lam3_s,lam4_s,lam5_s,lam6_s,lam7_s,c_sigma_s,kappa_s
  real(8) rho_s_i

  common /strucFun/ N_str
  common /gapPara/ lam1,lam2,lam3,lam4,lam5,lam6,lam7,c,kappa
  common /gapPara_s/ lam1_s,lam2_s,lam3_s,lam4_s,lam5_s,lam6_s,lam7_s,c_sigma_s,kappa_s
  common /drho1V_com/ drho1V,drho1V_s,mass_s
  common /rho_s_i_com/ rho_s_i


  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)


  lam1=yflow(1)
  lam2=yflow(2)
  lam3=yflow(3)
  lam4=yflow(4)
  lam5=yflow(5)
!  lam6=yflow(6)
!  lam7=yflow(7)
  lam6=0.
  lam7=0.
  lam0=yflow(Nv+1)
  h=yflow((Nv+1)+1)
  Zphi=yflow((Nv+1)+(Nh+1)+1)
  Zpsi=yflow((Nv+1)+(Nh+1)+2)
  ZA=yflow((Nv+1)+(Nh+1)+3)
  c=yflow((Nv+1)+(Nh+1)+Nz+1)
  kappa=yflow((Nv+1)+(Nh+1)+Nz+2)
  g=yflow((Nv+1)+(Nh+1)+Nz+Nck+1)


  x(1)=kappa*1.1       !avoiding numerical instability
  call newt(x, nn, check, gapEq)
  rho0=x(1)
!  rho0=kappa

!!!!!!!!! for reduced condesate, strange quark part
go to 100

!  ms0=120./hc
  ms0=190./hc
  rho_s=(sqrt(rho0/2.)+ms0/h)**2 !strange quark mass
  rho=2.*rho_s
  drho1V_s=lam1 + lam2*(-kappa + rho) + (lam3*(-kappa + rho)**2)/2. + &
           (lam4*(-kappa + rho)**3)/6. + (lam5*(-kappa + rho)**4)/24. + &
           (lam6*(-kappa + rho)**5)/120. + (lam7*(-kappa + rho)**6)/720.

  c_sigma_s=sqrt(2.*rho_s)*sqrt(2.)*drho1V_s    !explicit chiral symmetry breaking for strange quark
  c_bar_sigma_s=c_sigma_s*sqrt(Zphi)            !bare parameter

  write(*,"('c_bar_sigma_s=', e20.9)")c_bar_sigma_s
  write(*,"('rho_s=', e20.9)")rho_s

  stop

100 continue
!!!!!!!!!

!  c_bar_sigma_s=0.635253529E+04 ! ms0=120./hc
!  c_bar_sigma_s=0.100257386E+05 ! ms0=140./hc
!  c_bar_sigma_s=0.126905205E+05 ! ms0=150./hc
!  c_bar_sigma_s=0.140554E+05    ! ms0=155./hc
!  c_bar_sigma_s=0.161014503E+05 ! ms0=160./hc
!  c_bar_sigma_s=0.259224620E+05 ! ms0=180./hc
!  c_bar_sigma_s=0.328023504E+05 ! ms0=190./hc
!  c_bar_sigma_s=0.413777386E+05 ! ms0=200./hc
c_bar_sigma_s=100.798*(1.d3/hc)**3

  c_sigma_s=c_bar_sigma_s/sqrt(Zphi)

  lam1_s=lam1
  lam2_s=lam2
  lam3_s=lam3
  lam4_s=lam4
  lam5_s=lam5
  lam6_s=lam6
  lam7_s=lam7
  kappa_s=kappa

  x(1)=rho_s_i       !gap equation for strange quark
  call newt(x, nn, check, gapEq_s)
  rho_s=x(1)
  rho_s_i=rho_s

  drho1V_s=c_sigma_s/(sqrt(2.)*sqrt(2.*rho_s))

  mass_s=h*sqrt(2.*rho_s)/sqrt(2.)




!!!!!!!!!
  rho=rho0

  mPion2=lam1 + lam2*(-kappa + rho) + (lam3*(-kappa + rho)**2)/2. + &
        (lam4*(-kappa + rho)**3)/6. + (lam5*(-kappa + rho)**4)/24. + &
        (lam6*(-kappa + rho)**5)/120. + (lam7*(-kappa + rho)**6)/720.

  mPion=sqrt(mPion2)
  drho1V=mPion2

  mSigma2=lam1 + lam2*(-kappa + rho) + (lam3*(-kappa + rho)**2)/2. + &
        (lam4*(-kappa + rho)**3)/6. + (lam5*(-kappa + rho)**4)/24. + &
        (lam6*(-kappa + rho)**5)/120. + (lam7*(-kappa + rho)**6)/720. + &
        2*rho*(lam2 + lam3*(-kappa + rho) + (lam4*(-kappa + rho)**2)/2. + &
        (lam5*(-kappa + rho)**3)/6. + (lam6*(-kappa + rho)**4)/24. + &
        (lam7*(-kappa + rho)**5)/120.)

  mSigma=sqrt(mSigma2)

  mf=h*sqrt(rho/2.)

  Vall=lam0 - Sqrt(2.)*c*Sqrt(rho) + lam1*(-kappa + rho) + &
       (lam2*(-kappa + rho)**2)/2. + (lam3*(-kappa + rho)**3)/6. + &
       (lam4*(-kappa + rho)**4)/24. + (lam5*(-kappa + rho)**5)/120. + &
       (lam6*(-kappa + rho)**6)/720. + (lam7*(-kappa + rho)**7)/5040.

end

