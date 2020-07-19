subroutine FRG(kappa_UV_i,kappa_UV,rho0,mPion,mSigma,mf,Vall,fpi,h,Zphi,Zpsi,ZA,c,kappa,g)
!This subroutine solve FRG flow equations with fixed expansion point, search for phyiscal point, find
!appropriate expansion point. Inputing a guess kappa_UV_i, outputing kappa_UV and other physical variables

  implicit none

  integer Nv,Nh
  !Nv: order of Tylor expansion for effective potential V, Nh: order of Yukawa coupling h
  parameter(Nv=7)
  parameter(Nh=0)
  integer Nz 
  !number of wave function renormalizations
  parameter(Nz=3)
  integer Nck
  parameter(Nck=2)
  integer Ng
  parameter(Ng=3)
  integer Nflow 
  !number of flow equations
  parameter(Nflow=(Nv+1)+(Nh+1)+Nz+Nck+Ng)
  real(16) yflow(Nflow)
  !dependent variables in flow equations
  !integer N_str(4) !store the structure of functions of ODE
  integer N_str(5) 
  !store the structure of functions of ODE
  real(16) pi,hc
  parameter(pi=3.141592653589793Q+0)
  parameter(hc=197.33Q+0)
  real(16) kappa_UV_i,kappa_UV
  real(16) T,mu
  real(16) k_UV,k_IR,t_UV,t_IR
  external derivs,rkqs
  real(16) eps_ode,h1,hmin 
  !variables in subroutine odeint
  integer nok,nbad 
  !variables in subroutine odeint
  INTEGER kmax,kount 
  !variables in common block of subroutine odeint
  INTEGER KMAXX,NMAX
  PARAMETER (NMAX=50,KMAXX=2000)
  real(16) dxsav,xp(KMAXX),yp(NMAX,KMAXX) 
  !variables in common block of subroutine odeint
  real(16) rho0,mPion,mSigma,mf,Vall
  real(16) fpi,h,Zphi,Zpsi,ZA,c,kappa,g
  integer i
  real(16) l_com,lb_com
  integer k_num
  real(16) k_value
  real(16) g_max,alphas_max,k_g_max,ZAm1_max,k_ZAm1_max,g3A_max,alphas3A_max,k_g3A_max

  integer n_x
  parameter(n_x=64)
  real(16) w_x(n_x),y_x(n_x) 
  !Guass integral
  integer n_cth
  parameter(n_cth=64)
  real(16) w_cth(n_cth),y_cth(n_cth) 
  !Guass integral

  common /strucFun/ N_str
  common /Tmu/ T,mu
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /odeContr/ eps_ode,h1,hmin
  COMMON /path/ kmax,kount,dxsav,xp,yp
  common /polyakov_com/ l_com,lb_com
  common /k_num_com/k_num,k_value
  common /alphas_max_com/g_max,alphas_max,k_g_max,ZAm1_max,k_ZAm1_max,g3A_max,alphas3A_max,k_g3A_max

  common /gausslegFRGcom/w_x,y_x,w_cth,y_cth

  N_str(1)=Nv
  N_str(2)=Nh
  N_str(3)=Nz
  N_str(4)=Nck
  N_str(5)=Ng

  k_UV=20.Q+3/hc 
  !in unit of fm**(-1)
  k_IR=0.01Q+0/hc   
  !in unit of fm**(-1)
  t_UV=0.Q+0
  t_IR=log(k_IR/k_UV)

  If(T<90)then 
    eps_ode=5.Q-7
  else
    eps_ode=5.Q-9
  end if
  h1=t_IR/20000.Q+0
  hmin=0.Q+0
  kmax=KMAXX
  dxsav=t_IR/10000.Q+0

  kappa_UV=kappa_UV_i
  call initial(Nflow,yflow,kappa_UV)

  k_num=0
  k_value=0.Q+0
  g_max=0.Q+0
  alphas_max=0.Q+0
  k_g_max=0.Q+0
  g3A_max=0.Q+0
  alphas3A_max=0.Q+0
  k_g3A_max=0.Q+0
  ZAm1_max=0.Q+0
  k_ZAm1_max=0.Q+0

  call gauleg(0.Q+0,1.Q+0, y_x, w_x, n_x)
  call gauleg(-1.Q+0,1.Q+0, y_cth, w_cth, n_cth)

  call odeint(yflow,Nflow,t_UV,t_IR,eps_ode,h1,hmin,nok,nbad,derivs,rkqs)
  call phypoint(Nflow,yflow,rho0,mPion,mSigma,mf,Vall)

  fpi=sqrt(2.Q+0*rho0)
  h=yflow((Nv+1)+1)
  Zphi=yflow((Nv+1)+(Nh+1)+1)
  Zpsi=yflow((Nv+1)+(Nh+1)+2)
  ZA=yflow((Nv+1)+(Nh+1)+3)
  c=yflow((Nv+1)+(Nh+1)+Nz+1)
  kappa=yflow((Nv+1)+(Nh+1)+Nz+2)
  g=yflow((Nv+1)+(Nh+1)+Nz+Nck+1)

end
