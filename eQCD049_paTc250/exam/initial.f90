subroutine initial(Nflow,yflow,kappa)
!make the initialization

  implicit none
  integer Nflow
  real(16) yflow(Nflow) !sigma is the fixed expansion point at UV
  integer N_str(5) !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng
  real(16) lam0,lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(16) h
  real(16) Zphi,Zpsi,ZA
  real(16) c,kappa,kappa0
  real(16) alphaS,g,b0,MU,alphaSMU,alphaS1L
  real(16) pi,hc
  parameter(pi=3.141592653589793Q+0)
  parameter(hc=197.33Q+0)
  real(16) lambda,nu
  real(16) k_UV,k_IR,t_UV,t_IR
  real(16) Nc,Nf
  parameter(Nc=3.Q+0,Nf=2.Q+0)

  common /strucFun/ N_str
  common /kRange/k_UV,k_IR,t_UV,t_IR

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)

  b0=11.Q+0-2.Q+0/3.Q+0*Nf
  MU=91.Q+0                       !GeV
  MU=MU*1.Q+3/hc                !fm**(0=-1)
  alphaSMU = 0.118Q+0

  alphaS1L=alphaSMU/(1 + (b0*alphaSMU*Log(k_UV/MU))/(2.Q+0*Pi))
  alphaS=alphaS1L
  alphaS=0.235Q+0


  g=sqrt(4.Q+0*pi*alphaS)

  lambda=0.Q+0
  nu=1.Q+4*k_UV**2

  h=1.Q+0
  !expansion coefficients of Yukawa coupling

  c=3.6Q+0*(1.Q+3/hc)**3  
  !explicit chiral symmetry breaking term, in unit of fm**(-3)




  kappa0=(c/nu)**2/2.Q+0
  !minimal point
!  kappa0=kappa



  lam0=0.Q+0
  lam1=kappa0*lambda+nu
  lam2=lambda
  lam3=0.Q+0
  lam4=0.Q+0
  lam5=0.Q+0
  lam6=0.Q+0
  lam7=0.Q+0
!expansion coefficients of effective potential V

  Zphi=1.Q+0 
  !meson wave function renormalization
  Zpsi=1.Q+0 
  !quark wave function renormalization
  ZA=1.Q+0   
  !gluon wave function renormalization

  yflow(1)=lam1
  yflow(2)=lam2
  yflow(3)=lam3
  yflow(4)=lam4
  yflow(5)=lam5
  yflow(6)=lam6
  yflow(7)=lam7
  yflow(Nv+1)=lam0
  yflow((Nv+1)+1)=h
  yflow((Nv+1)+(Nh+1)+1)=Zphi
  yflow((Nv+1)+(Nh+1)+2)=Zpsi
  yflow((Nv+1)+(Nh+1)+3)=ZA
  yflow((Nv+1)+(Nh+1)+Nz+1)=c
!  yflow((Nv+1)+(Nh+1)+Nz+2)=kappa
  yflow((Nv+1)+(Nh+1)+Nz+2)=kappa0
  yflow((Nv+1)+(Nh+1)+Nz+Nck+1)=g
  yflow((Nv+1)+(Nh+1)+Nz+Nck+2)=g   
  !Three gluon interaction
  yflow((Nv+1)+(Nh+1)+Nz+Nck+3)=g   
  !quark-gluon interaction for the strange quark

end
