program eQCD

  implicit none

  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) T,mu 
  !temperature and chemical potential
  real(16) l,lb 
  !polyakov loop
  real(16) rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,ZA
  real(16) c,kappa,g
  real(16) Zphi_p0
  real(16) Vtotal0
  real(16) Ti,dT,mu_down,mu_up
  integer i,iTmax,j,jmumax,m,mm,i_mm,j1
  parameter(iTmax=300,jmumax=0)
  real(16) mu_bound(0:iTmax),mu_bound_low,mu_bound_high
  real(16) T_MeV,muB_MeV
  real(16) mui,muBi,muB
  integer mmax
  parameter(mmax=10)
  !order of chebyshev polynomial
  real(16) dcd0mu(mmax),dcd1mu(mmax),dcd2mu(mmax),dcd3mu(mmax),dcd4mu(mmax),chi(mmax,0:iTmax)
  real(16) factorial
  real(16) fpi_res(0:jmumax,0:iTmax),mPion_res(0:jmumax,0:iTmax),mSigma_res(0:jmumax,0:iTmax),mf_res(0:jmumax,0:iTmax)
  integer iT,iv
  real(16) l_i1,l_i2,lb_i1,lb_i2,l_i_mu,lb_i_mu
  real(16) nfl,nfl0,nfl1,lset

  common /Tmu/ T,mu
  common /Zphi_p0_com/ Zphi_p0

  dT=1.Q+0

  muBi=0.Q+00/hc

  open(unit=11,file='m1.dat')
  read(11,*) j1
  close (11)
  write(*,*)'load OK',j1

  l_i1=1.Q-10
  l_i2=l_i1
  lb_i1=1.Q-10
  lb_i2=lb_i1

  do i=1,iTmax
    T_MeV=dT*real(i,kind=16)
    T=T_MeV/hc

    muB=muBi
    mu=muB/3.Q+0
    muB_MeV=muB*hc

    call selfEQ(l_i1,lb_i1,l,lb,rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,ZA,c,kappa,g)
    l_i1=2.Q+0*l-l_i2
    l_i2=l
    lb_i1=2.Q+0*lb-lb_i2
    lb_i2=lb

    write(*,"('i=', I4)")i
    write(*,"('muB_MeV=', f15.7, t25, 'T_MeV=', f15.7)")muB_MeV,T_MeV
    write(*,"('fpi=', f15.7, t25, 'mPion=', f15.7)")fpi*hc,mPion*hc
    write(*,*)'mpi_phy:',mPion*hc*sqrt(Zphi/Zphi_p0)
    write(*,*)'msg_phy:',mSigma*hc*sqrt(Zphi/Zphi_p0)

    open(unit=51,file='./buffer/TMeV.dat',position='append')
    write(51,*)T_MeV
    close(51)

    open(unit=51,file='./buffer/l.dat',position='append')
    write(51,*)l
    close(51)

    open(unit=51,file='./buffer/lb.dat',position='append')
    write(51,*)lb
    close(51)

    open(unit=51,file='./buffer/rho0.dat',position='append')
    write(51,*)rho0
    close(51)

    open(unit=51,file='./buffer/mPion.dat',position='append')
    write(51, "(e21.14)")mPion*hc
    close(51)

    open(unit=51,file='./buffer/mSigma.dat',position='append')
    write(51,*)mSigma*hc
    close(51)

    open(unit=51,file='./buffer/mf.dat',position='append')
    write(51,*)mf*hc
    close(51)

    open(unit=51,file='./buffer/Vtotal.dat',position='append')
    write(51,*)Vtotal
    close(51)

    open(unit=51,file='./buffer/fpi.dat',position='append')
    write(51,*)fpi*hc
    close(51)

    open(unit=51,file='./buffer/h.dat',position='append')
    write(51,*)h
    close(51)

    open(unit=51,file='./buffer/Zphi.dat',position='append')
    write(51, "(e21.14)")Zphi
    close(51)

    open(unit=51,file='./buffer/Zpsi.dat',position='append')
    write(51, "(e21.14)")Zpsi
    close(51)

    open(unit=51,file='./buffer/ZA.dat',position='append')
    write(51, *)ZA
    close(51)

    open(unit=51,file='./buffer/c.dat',position='append')
    write(51, *)c
    close(51)

    open(unit=51,file='./buffer/kappa.dat',position='append')
    write(51, *)kappa
    close(51)

    open(unit=51,file='./buffer/g.dat',position='append')
    write(51, *)g
    close(51)

    open(unit=51,file='./buffer/muB.dat',position='append')
    write(51, "(e21.14)")mu*3.Q+0*hc
    close(51)

    open(unit=51,file='./buffer/mpion_phy.dat',position='append')
    write(51, "(e21.14)")mPion*hc*sqrt(Zphi/Zphi_p0)
    close(51)

    open(unit=51,file='./buffer/mSigma_phy.dat',position='append')
    write(51, "(e21.14)")mSigma*hc*sqrt(Zphi/Zphi_p0)
    close(51)

    open(unit=51,file='./buffer/Zphi_p0.dat',position='append')
    write(51, "(e21.14)")Zphi_p0
    close(51)
  end do

end
