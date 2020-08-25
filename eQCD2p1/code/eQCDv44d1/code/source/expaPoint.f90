subroutine expaPoint(Nflow,yflow,kappa_UV)

  implicit none

  integer Nflow
  real(8) yflow(Nflow) !sigma is the fixed expansion point at UV
  integer N_str(5) !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) kappa_UV, kappa_old, kappa_new
  real(8) T,mu
  real(8) k_UV,k_IR,t_UV,t_IR
  external derivs,rkqs
  real(8) eps_ode,h1,hmin !variables in subroutine odeint
  integer nok,nbad !variables in subroutine odeint
  INTEGER kmax,kount !variables in common block of subroutine odeint
  INTEGER KMAXX,NMAX
  PARAMETER (NMAX=50,KMAXX=2000)
  real(8) dxsav,xp(KMAXX),yp(NMAX,KMAXX) !variables in common block of subroutine odeint
  real(8) rho0,mPion,mSigma,mf,Vall
  real(8) fpi,h,Zphi,Zpsi,ZA,c,g
  real(8) epsi_rho0,epsi
  real(8) rho0_gaug
  real(8) kappa
  logical stopp
  integer i
  integer imax  !maximal number of loops
!  parameter(imax=600)
  parameter(imax=100)
  real(8) rescal,delta_IR
  integer k_num1
  real(8) k_value1


  common /strucFun/ N_str
  common /Tmu/ T,mu
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /odeContr/ eps_ode,h1,hmin
  COMMON /path/ kmax,kount,dxsav,xp,yp
  common /k_num1_com/k_num1,k_value1


  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)


  !  epsi_rho0=0.00001
  epsi_rho0=0.
  epsi=1.e-2
  rescal=20.

  kappa_old=kappa_UV
  kappa_new=kappa_UV

  i=0                    !start of loops
  stopp=.false.
  do while((.not. stopp).and.(i < imax))
     i=i+1

     kappa_old=kappa_new
     call initial(Nflow,yflow,kappa_old)
     k_num1=0
     k_value1=0.
     call odeint(yflow,Nflow,t_UV,t_IR,eps_ode,h1,hmin,nok,nbad,derivs,rkqs)
     call phypoint(Nflow,yflow,rho0,mPion,mSigma,mf,Vall)



     fpi=sqrt(2.*rho0)
     h=yflow((Nv+1)+1)
     Zphi=yflow((Nv+1)+(Nh+1)+1)
     Zpsi=yflow((Nv+1)+(Nh+1)+2)
     ZA=yflow((Nv+1)+(Nh+1)+3)
     c=yflow((Nv+1)+(Nh+1)+Nz+1)
     kappa=yflow((Nv+1)+(Nh+1)+Nz+2)
     g=yflow((Nv+1)+(Nh+1)+Nz+Nck+1)



     write(*,"('kappa/rho0=', f21.14)")kappa/rho0

!  goto 110

     write(*,"('kappa_UV=', f21.14)")kappa_old
     write(*,"('kappa=', f21.14)")kappa
     write(*,"('rho0=', f21.14)")rho0
     write(*,"('fpi=', f21.14)")fpi*hc
     write(*,"('mPion=', f21.14)")mPion*hc
     write(*,"('mSigma=', f21.14)")mSigma*hc
     write(*,"('mf=', f21.14)")mf*hc
     write(*,"('h=', f21.14)")h
     write(*,"('Zphi=', e20.9)")Zphi

   stop 



     open(unit=51,file='./buffer/yflow.dat')
     do i=1, Nflow
        write(51, "(e20.9)")yflow(i)
     end do
     close(51)

     open(unit=51,file='./buffer/t.dat')
     do i=1, kount
        write(51, "(e20.9)")xp(i)
     end do
     close(51)

     open(unit=51,file='./buffer/k.dat')
     do i=1, kount
        write(51, "(e20.9)")k_UV*exp(xp(i))*hc
     end do
     close(51)


     open(unit=51,file='./buffer/lam1.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(1,i)
     end do
     close(51)

     open(unit=51,file='./buffer/lam2.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(2,i)
     end do
     close(51)

     open(unit=51,file='./buffer/lam3.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(3,i)
     end do
     close(51)

     open(unit=51,file='./buffer/lam4.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(4,i)
     end do
     close(51)

     open(unit=51,file='./buffer/lam5.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(5,i)
     end do
     close(51)

     open(unit=51,file='./buffer/lam6.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(6,i)
     end do
     close(51)

     open(unit=51,file='./buffer/lam7.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(7,i)
     end do
     close(51)

     open(unit=51,file='./buffer/lam0.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(8,i)
     end do
     close(51)

     open(unit=51,file='./buffer/h.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(9,i)
     end do
     close(51)

     open(unit=51,file='./buffer/Zphi.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(10,i)
     end do
     close(51)

     open(unit=51,file='./buffer/Zpsi.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(11,i)
     end do
     close(51)

     open(unit=51,file='./buffer/ZA.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(12,i)
     end do
     close(51)

     open(unit=51,file='./buffer/invZA.dat')
     do i=1, kount
        write(51, "(e20.9)")1./yp(12,i)
     end do
     close(51)

     open(unit=51,file='./buffer/c.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(13,i)
     end do
     close(51)

     open(unit=51,file='./buffer/kappa.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(14,i)
     end do
     close(51)

     open(unit=51,file='./buffer/g.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(15,i)
     end do
     close(51)

     open(unit=51,file='./buffer/alphaS.dat')
     do i=1, kount
        write(51, "(e20.9)")yp(15,i)**2/(4.*pi)
     end do
     close(51)

     open(unit=51,file='./buffer/massPion.dat')
     do i=1, kount
        write(51, "(e20.9)")sqrt(yp(1,i))*hc
     end do
     close(51)

     open(unit=51,file='./buffer/massSigma.dat')
     do i=1, kount
        write(51, "(e20.9)")sqrt(yp(1,i) + 2*yp(2,i)*yp(14,i))*hc
     end do
     close(51)

     open(unit=51,file='./buffer/massQuark.dat')
     do i=1, kount
        write(51, "(e20.9)")sqrt((yp(9,i)**2*yp(14,i))/2.)*hc
     end do
     close(51)




     open(unit=51,file='./buffer/nok.dat')
     write(51, "(I4)")nok
     close(51)

     open(unit=51,file='./buffer/nbad.dat')
     write(51, "(I4)")nbad
     close(51)

     stop
110  continue


     kappa=yflow((Nv+1)+(Nh+1)+Nz+2)
!    rho0_gaug=rho0*(1.+epsi_rho0)
     rho0_gaug=rho0+epsi_rho0

     delta_IR=kappa-rho0_gaug
     if(abs(delta_IR)/kappa<epsi)then
        stopp=.true.
     else
        if(sqrt(2.*kappa)*hc>1.)then
          kappa_new=kappa_old-delta_IR/(Zphi*rescal)
        else
          stopp=.true.
        end if
     end if


!     write(*,"('delta_IR=', f15.7)")delta_IR
!     write(*,"('abs(delta_IR)/kappa=', f15.7)")abs(delta_IR)/kappa

  end do

  kappa_UV=kappa_old

end subroutine expaPoint



