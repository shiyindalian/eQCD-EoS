program EQCD

  implicit none

  real(16) pi,hc
  parameter(pi=3.141592653589793Q+0)
  parameter(hc=197.33Q+0)
  real(16) T,mu 
  !temperature and chemical potential
  real(16) l,lb 
  !polyakov loop
  real(16) rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,ZA,c,kappa,g
  real(16) sigma_UV,kappa_UV_i,kappa_UV_i_mu,kappa_UV
  real(16) Vtotal0
  real(16) Ti,dT,mu_down,mu_up
  integer i,iTmax,j,jmumax,m,mm,i_mm,j1,munum
  parameter(iTmax=200,jmumax=0,munum=99)
  real(16) pre_res(0:jmumax,0:iTmax),T_res(0:iTmax),mu_res(0:jmumax,0:iTmax),pre_com(jmumax)
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

  real(16) g_max,alphas_max,k_g_max,ZAm1_max,k_ZAm1_max,g3A_max,alphas3A_max,k_g3A_max
  real(16) Zphi_p0
  real(16) mudelta

  common /Tmu/ T,mu
  common /alphas_max_com/g_max,alphas_max,k_g_max,ZAm1_max,k_ZAm1_max,g3A_max,alphas3A_max,k_g3A_max


  dT=1.Q+0

  muBi=0.Q+0/hc

  open(unit=11,file='m1.dat')
  read(11,*) j1
  close (11)
  write(*,*)'load OK',j1
  
  j=j1

  sigma_UV=2.8Q-2/hc
  kappa_UV_i=sigma_UV**2/2.Q+0

  l_i1=1.Q-10
  l_i2=l_i1
  lb_i1=1.Q-10
  lb_i2=lb_i1

  mu_up=250.Q+0/hc
  mu_down=-mu_up

  mudelta=cos(pi*(j-0.5Q+00)/real(munum,kind=16))*(0.5Q+00*(mu_up-mu_down))+0.5Q+00*(mu_up+mu_down)

  do i=1,iTmax
    T_MeV=dT*real(i,kind=16)
    T=T_MeV/hc

    muB=muBi+mudelta
    mu=muB/3.Q+0
    muB_MeV=muB*hc

    call selfEQ(kappa_UV_i,l_i1,lb_i1,kappa_UV,l,lb,rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,ZA,c,kappa,g)
    kappa_UV_i=kappa_UV
    l_i1=2.Q+0*l-l_i2
    l_i2=l
    lb_i1=2.Q+0*lb-lb_i2
    lb_i2=lb

    write(*,"('i=', I4)")i
    write(*,"('muB_MeV=', f15.7, t25, 'T_MeV=', f15.7)")muB_MeV,T_MeV
    write(*,"('fpi=', f15.7, t25, 'mPion=', f15.7)")fpi*hc,mPion*hc
    write(*,"('kappa_UV=', e21.14)")kappa_UV

    open(unit=51,file='./buffer/TMeV.dat',position='append')
    write(51,*)T_MeV
    close(51)

    open(unit=51,file='./buffer/l.dat',position='append')
    write(51,*)l
    close(51)

    open(unit=51,file='./buffer/lb.dat',position='append')
    write(51,*)lb
    close(51)

    open(unit=51,file='./buffer/kappaUV.dat',position='append')
    write(51,*)kappa_UV
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

    open(unit=51,file='./buffer/alphasmax.dat',position='append')
    write(51, "(e21.14)")alphas_max
    close(51)

    open(unit=51,file='./buffer/kgmax.dat',position='append')
    write(51, "(e21.14)")k_g_max
    close(51)

    open(unit=51,file='./buffer/g3A_max.dat',position='append')
    write(51, "(e21.14)")g3A_max
    close(51)

    open(unit=51,file='./buffer/alphas3A_max.dat',position='append')
    write(51, "(e21.14)")alphas3A_max
    close(51)

    open(unit=51,file='./buffer/k_g3A_max.dat',position='append')
    write(51, "(e21.14)")k_g3A_max
    close(51)

    open(unit=51,file='./buffer/ZAm1max.dat',position='append')
    write(51, "(e21.14)")ZAm1_max
    close(51)

    open(unit=51,file='./buffer/kZAm1max.dat',position='append')
    write(51, "(e21.14)")k_ZAm1_max
    close(51)

    open(unit=51,file='./buffer/muB.dat',position='append')
    write(51, "(e21.14)")mu*3.d0*hc
    close(51)
  end do

end
