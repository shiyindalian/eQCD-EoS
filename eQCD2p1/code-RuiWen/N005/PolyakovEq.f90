subroutine PolyakovEq(n, x, fvec)

  implicit none

  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  integer n
  real(16) x(n), fvec(n)
  real(16) l,lb
  integer i,j
  integer npoint1, ii
  parameter(npoint1=1024)
  real(16) w1(npoint1),y1(npoint1)
  !Guass integral

  real(16) yflowA(npoint1,50),yflowk(50)
  real(16) kA(npoint1)
  real(16) etaphiA(npoint1),etapsiA(npoint1),etaAA(npoint1)

  real(16) kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ZAk,ck,kappak,gk,etaphik,etapsik,etaAk,rhok

  real(16) rho0k0,Zphik0
  real(16) rho0
  real(16) dtVd1l,dtVd1lb
  real(16) Vd1l_1,Vd1lb_1,Vd1l_2,Vd1lb_2,Vd1l_4,Vd1lb_4,Vd1l,Vd1lb
  real(16) T,mu
  real(16) a1,a2,a3,a4,a5,b1,b2,b3,b4,c1,c2,c3,c4,c5,d1,d2,d3,d4,d5,T0,T_r,aPolya,cPolya,dPolya,bPolya

  common /y1w1/ y1, w1
  common /flowconfi/ kA,yflowA,etaphiA,etapsiA,etaAA
  common /rho0Zphi/ rho0k0,Zphik0
  common /Tmu/ T,mu

  l=x(1)
  lb=x(2)

  Vd1l_1=0.Q+0
  Vd1lb_1=0.Q+0
  do ii=1, npoint1
    kk=kA(ii)
    yflowk=yflowA(ii,:)
    etaphik=etaphiA(ii)
    etapsik=etapsiA(ii)

    Zphik=yflowA(ii,23)
    rhok=rho0k0*Zphik/Zphik0

    call dtVdiff1(kk,yflowk,etaphik,etapsik,rhok,l,lb,dtVd1l,dtVd1lb)

    Vd1l_1=Vd1l_1+w1(ii)/kk*dtVd1l
    Vd1lb_1=Vd1lb_1+w1(ii)/kk*dtVd1lb

  end do

  Vd1l_2=0.Q+0
  Vd1lb_2=0.Q+0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!glue potential
  a1=-44.14Q+0
  a2=151.4Q+0
  a3=-90.0677Q+0
  a4=2.77173Q+0
  a5=3.56403Q+0

  b1=-0.32665Q+0
  b2=-82.9823Q+0
  b3=3.Q+0
  b4=5.85559Q+0

  c1=-50.7961Q+0
  c2=114.038Q+0
  c3=-89.4596Q+0
  c4=3.08718Q+0
  c5=6.72812Q+0

  d1=27.0885Q+0
  d2=-56.0859Q+0
  d3=71.2225Q+0
  d4=2.9715Q+0
  d5=6.61433Q+0

  T0=225.Q+0/hc

  T_r=0.57Q+0*(T-T0)/T0+1.Q+0

  aPolya=(a1+a2/T_r+a3/T_r**2)/(1.Q+0+a4/T_r+a5/T_r**2)
  cPolya=(c1+c2/T_r+c3/T_r**2)/(1.Q+0+c4/T_r+c5/T_r**2)
  dPolya=(d1+d2/T_r+d3/T_r**2)/(1.Q+0+d4/T_r+d5/T_r**2)

  bPolya=b1*T_r**(-b4)*(1.Q+0-exp(b2/T_r**b3))

  Vd1l_4=((3*cPolya*l**2)/2.Q+0 - (aPolya*lb)/2.Q+0 + 2*dPolya*l*lb**2 +       &
    (bPolya*(12*l**2 - 6*lb - 6*l*lb**2))/                                     &
     (1 - 6*l*lb - 3*l**2*lb**2 + 4*(l**3 + lb**3)))*T**4

  Vd1lb_4=(-(aPolya*l)/2. + 2*dPolya*l**2*lb + (3*cPolya*lb**2)/2.Q+0 +        &
    (bPolya*(-6*l - 6*l**2*lb + 12*lb**2))/                                    &
     (1 - 6*l*lb - 3*l**2*lb**2 + 4*(l**3 + lb**3)))*T**4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Vd1l=Vd1l_1+Vd1l_2+Vd1l_4
  Vd1lb=Vd1lb_1+Vd1lb_2+Vd1lb_4

  fvec(1)=Vd1l
  fvec(2)=Vd1lb

end
