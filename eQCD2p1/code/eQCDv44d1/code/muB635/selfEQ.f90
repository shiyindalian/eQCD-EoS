subroutine selfEQ(kappa_UV_i,l_i,lb_i,kappa_UV,l,lb,rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,ZA,c,kappa,g)
!This subroutine solve the FRG equation and the equation for the Polyakov loop self-consistently,
!Inputing guess kappa_UV_i,l_i,lb_i, outputing kappa_UV,l,lb and other physical variables

  implicit none
  real(8) kappa_UV_i,l_i,lb_i,kappa_UV,l,lb
  real(8) rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,ZA,c,kappa,g
  real(8) l_new,lb_new,l_old,lb_old
  real(8) l_com,lb_com
  real(8) T,mu !temperature and chemical potential
  integer i,j
  integer jmax  !maximal number of loops
  parameter(jmax=600)
  INTEGER kmax,kount !variables in common block of subroutine odeint
  INTEGER KMAXX,NMAX
  PARAMETER (NMAX=50,KMAXX=2000)
  real(8) dxsav,xp(KMAXX),yp(NMAX,KMAXX) !variables in common block of subroutine odeint
  real(8) x_bef(KMAXX),y_bef(KMAXX) !used for interpolation
  real(8) k_UV,k_IR,t_UV,t_IR
  real(8) epsi_l,epsi_lb,epsi_err
  logical stopp
  real(8) Vall,Vinfi,Vglue
  external PolyakovEq
  integer nn
  parameter(nn=2)
  real(8) x(nn)
  integer npoint1, ii
  parameter(npoint1=512)
  real(8) w1(npoint1),y1(npoint1) !Guass integral
  integer npoint2
  parameter(npoint2=256)
  real(8) w2(npoint2),y2(npoint2) !Guass integral
  real(8) x_aft(npoint1),y_aft(npoint1) !used for interpolation

  real(8) kA(npoint1),lam1A(npoint1),lam2A(npoint1),lam3A(npoint1),lam4A(npoint1),lam5A(npoint1),lam0A(npoint1),&
          hA(npoint1),ZphiA(npoint1),ZpsiA(npoint1),ZAA(npoint1),cA(npoint1),kappaA(npoint1),gA(npoint1),&
          g3AA(npoint1),etaphiA(npoint1),etapsiA(npoint1),etaAA(npoint1)
  real(8) kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ZAk,ck,kappak,gk,g3Ak,etaphik,etapsik,etaAk
  real(8) rho0k0,Zphik0
  real(8) a1,a2,a3,a4,a5,b1,b2,b3,b4,c1,c2,c3,c4,c5,d1,d2,d3,d4,d5,T0,T_r,aPolya,cPolya,dPolya,bPolya
  logical check
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) k_infi


  COMMON /path/ kmax,kount,dxsav,xp,yp
  common /polyakov_com/ l_com,lb_com
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /y1w1/ y1, w1
  common /y2w2/ y2, w2
  common /flowconfi/ kA,lam1A,lam2A,lam3A,lam4A,lam5A,lam0A,hA,ZphiA,ZpsiA,ZAA,cA,kappaA,gA,g3AA,etaphiA,etapsiA,etaAA
  common /rho0Zphi/ rho0k0,Zphik0
  common /Tmu/ T,mu


  epsi_err=1.e-10

  l_new=l_i
  lb_new=lb_i

  j=0                    !start of loops
  stopp=.false.
  do while((.not. stopp).and.(j < jmax))
    j=j+1

    l_old=l_new
    lb_old=lb_new

    l_com=l_old
    lb_com=lb_old

    call FRG(kappa_UV_i,kappa_UV,rho0,mPion,mSigma,mf,Vall,fpi,h,Zphi,Zpsi,ZA,c,kappa,g)
    rho0k0=rho0
    Zphik0=Zphi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!interpolation
    if(j==1)then
      call gauleg(k_UV, k_IR, y1, w1, npoint1)

      k_infi=k_UV*2.
      call gauleg(k_UV, k_infi, y2, w2, npoint2)
    end if


    do ii=1, npoint1
      x_aft(ii)=y1(ii)
      kA(ii)=x_aft(ii)
    end do

    do i=1, kount
      x_bef(i)=k_UV*exp(xp(i))
    end do

    do i=1, kount
      y_bef(i)=yp(1,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      lam1A(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(2,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      lam2A(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(3,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      lam3A(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(4,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      lam4A(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(5,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      lam5A(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(8,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      lam0A(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(9,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      hA(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(10,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      ZphiA(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(11,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      ZpsiA(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(12,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      ZAA(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(13,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      cA(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(14,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      kappaA(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(15,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      gA(ii)=y_aft(ii)
    end do

    do i=1, kount
      y_bef(i)=yp(16,i)
    end do
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    do ii=1, npoint1
      g3AA(ii)=y_aft(ii)
    end do


    do ii=1, npoint1
      kk=kA(ii)
      lam1k=lam1A(ii)
      lam2k=lam2A(ii)
      lam3k=lam3A(ii)
      lam4k=lam4A(ii)
      lam5k=lam5A(ii)
      lam0k=lam0A(ii)
      hk=hA(ii)
      Zphik=ZphiA(ii)
      Zpsik=ZpsiA(ii)
      ZAk=ZAA(ii)
      ck=cA(ii)
      kappak=kappaA(ii)
      gk=gA(ii)
      g3Ak=g3AA(ii)

      call eta(kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ZAk,ck,kappak,gk,g3Ak,etaphik,etapsik,etaAk)

      etaphiA(ii)=etaphik
      etapsiA(ii)=etapsik
      etaAA(ii)=etaAk

!      etaphiA(ii)=0.
!      etapsiA(ii)=0.

    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    x(1)=l_com
    x(2)=lb_com
    call newt(x, nn, check, PolyakovEq)
    l_new=x(1)
    lb_new=x(2)
!    l_new=1.
!    lb_new=1.

    epsi_l=abs(l_new-l_old)/abs(l_new)
    epsi_lb=abs(lb_new-lb_old)/abs(lb_new)

    if(epsi_l<epsi_err.and.epsi_lb<epsi_err.and.j>=4)then
!    if(epsi_l<epsi_err.and.epsi_lb<epsi_err)then
      stopp=.true.                                           !循环控制
    end if

   write(*,"('polya loop number=', I4)")j
   write(*,"('epsi_l=', e20.9)")epsi_l
   write(*,"('epsi_lb=', e20.9)")epsi_lb
   write(*,"('l_old=', f15.12)")l_old
   write(*,"('l_new=', f15.12)")l_new

  end do

  l=l_old
  lb=lb_old

  l_com=l
  lb_com=lb

!  call vInf(k_UV,Vinfi)
  Vinfi=0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!glue potential

  a1=-44.14
  a2=151.4
  a3=-90.0677
  a4=2.77173
  a5=3.56403

  b1=-0.32665
  b2=-82.9823
  b3=3.
  b4=5.85559

  c1=-50.7961
  c2=114.038
  c3=-89.4596
  c4=3.08718
  c5=6.72812

  d1=27.0885
  d2=-56.0859
  d3=71.2225
  d4=2.9715
  d5=6.61433

!  T0=208./hc
!  T0=250./hc
  T0=225./hc


!  T_r=T/T0
  T_r=0.57*(T-T0)/T0+1.

  aPolya=(a1+a2/T_r+a3/T_r**2)/(1.+a4/T_r+a5/T_r**2)
  cPolya=(c1+c2/T_r+c3/T_r**2)/(1.+c4/T_r+c5/T_r**2)
  dPolya=(d1+d2/T_r+d3/T_r**2)/(1.+d4/T_r+d5/T_r**2)

  bPolya=b1*T_r**(-b4)*(1.-exp(b2/T_r**b3))

  Vglue=T**4*(-(aPolya*l*lb)/2. + dPolya*l**2*lb**2 + (cPolya*(l**3 + lb**3))/2. +  &
    bPolya*Log(1 - 6*l*lb - 3*l**2*lb**2 + 4*(l**3 + lb**3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Vtotal=Vall+Vinfi+Vglue

end








