program EQCD

  implicit none

  real(8) pi,hc
  parameter(pi=3.141592653589793Q+0)
  parameter(hc=197.33Q+0)
  real(8) T,mu 
  !temperature and chemical potential
  real(8) l,lb 
  !polyakov loop
  real(8) rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,ZA,c,kappa,g
  real(8) sigma_UV,kappa_UV_i,kappa_UV_i_mu,kappa_UV
  real(8) Vtotal0,cUV,cUV0
  real(8) Ti,dT,mu_down,mu_up
  integer i,iTmax,j,jmumax,m,mm,i_mm,j1,munum
  parameter(iTmax=250,jmumax=0,munum=101)
  real(8) pre_res(0:jmumax,0:iTmax),T_res(0:iTmax),mu_res(0:jmumax,0:iTmax),pre_com(jmumax)
  real(8) mu_bound(0:iTmax),mu_bound_low,mu_bound_high
  real(8) T_MeV,muB_MeV
  real(8) mui,muBi,muB
  integer mmax
  parameter(mmax=10)
  !order of chebyshev polynomial
  real(8) dcd0mu(mmax),dcd1mu(mmax),dcd2mu(mmax),dcd3mu(mmax),dcd4mu(mmax),chi(mmax,0:iTmax)
  real(8) factorial
  real(8) fpi_res(0:jmumax,0:iTmax),mPion_res(0:jmumax,0:iTmax),mSigma_res(0:jmumax,0:iTmax),mf_res(0:jmumax,0:iTmax)
  integer iT,iv
  real(8) l_i1,l_i2,lb_i1,lb_i2,l_i_mu,lb_i_mu
  real(8) nfl,nfl0,nfl1,lset

  real(8) g_max,alphas_max,k_g_max,ZAm1_max,k_ZAm1_max,g3A_max,alphas3A_max,k_g3A_max
  real(8) Zphi_p0
  real(8) mudelta

  common /Tmu/ T,mu
  common /alphas_max_com/g_max,alphas_max,k_g_max,ZAm1_max,k_ZAm1_max,g3A_max,alphas3A_max,k_g3A_max
  common /cUV/ cUV
  common /Zphi_p0_com/ Zphi_p0


  dT=1.Q+0

  mu=0.Q+0/hc
  cUV0=0.1Q+0*(1.Q+3/hc)**3 - (3.8Q+0-0.1Q+0)/101.*(1.Q+3/hc)**3  

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

  !mu_up=0.2Q+0*(1.Q+3/hc)**3 
  !mu_down=-mu_up

  !mudelta=cos(pi*(j-0.5Q+00)/real(munum,kind=16))*(0.5Q+00*(mu_up-mu_down))+0.5Q+00*(mu_up+mu_down)
  mudelta=(3.8Q+0-0.1Q+0)/101.*(1.Q+3/hc)**3 *j

  do i=1,iTmax
    T_MeV=dT*real(i,kind=16)
    T=T_MeV/hc

    cUV=cUV0+mudelta

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

    open(unit=51,file='./buffer/sigmaEoM.dat',position='append')
    write(51,*)mf*hc
    close(51)

  if(i==1)then
        open(unit=51,file='./buffer/mpion_phy.dat',position='append')
        write(51, "(e21.14)")mPion*hc*sqrt(Zphi/Zphi_p0)
        close(51)
  end if

  end do

end
