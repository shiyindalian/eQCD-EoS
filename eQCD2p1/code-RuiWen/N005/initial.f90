subroutine initial(yflow)
!make the initialization

  implicit none
  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) Nc,Nf
  parameter(Nc=3.Q+0,Nf=2.Q+0)
  real(16) yflow(50) 
  integer N_str(5) 
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng
  real(16) lam00,lam10,lam20,lam30,lam40,lam50,lam60,lam70
  real(16) lam01,lam11,lam21,lam31,lam41,lam51
  real(16) lam02,lam12,lam22,lam32
  real(16) lam03,lam13
  real(16) h
  real(16) Z_pi,Z_k,Z_l,Z_s,Z_A,Zphi_p0
  real(16) jl,js,ck,Sl,Ss
  real(16) alphaS,g
  real(16) nu
  real(16) k_UV,k_IR,t_UV,t_IR

  common /strucFun/ N_str
  common /kRange/k_UV,k_IR,t_UV,t_IR

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)

  alphaS=0.2Q+0

  g=sqrt(4.Q+0*pi*alphaS)

  nu=1.Q+4*k_UV**2

  h=1.Q+0
  !expansion coefficients of Yukawa coupling

  ck=0.Q+0
  jl=3.6Q+0*(1.Q+3/hc)**3  
  js=27.Q+0*Jl 
  !explicit chiral symmetry breaking term, in unit of fm**(-3)

  Sl=jl/nu
  Ss=js/nu
  !minimal point

  lam00=0.Q+0
  lam10=nu
  lam20=0.Q+0
  lam30=0.Q+0
  lam40=0.Q+0
  lam50=0.Q+0
  lam60=0.Q+0
  lam70=0.Q+0
  lam01=0.Q+0
  lam11=0.Q+0
  lam21=0.Q+0
  lam31=0.Q+0
  lam41=0.Q+0
  lam51=0.Q+0
  lam02=0.Q+0
  lam12=0.Q+0
  lam22=0.Q+0
  lam32=0.Q+0
  lam03=0.Q+0
  lam13=0.Q+0
  !expansion coefficients of effective potential V

  Z_pi=1.Q+0 
  Z_K=1.Q+0 
  !meson wave function renormalization
  Z_l=1.Q+0 
  Z_s=1.Q+0 
  !quark wave function renormalization
  Z_A=1.Q+0   
  !gluon wave function renormalization
  Zphi_p0=1.Q+0

  yflow=0.Q+0

  yflow(1)=lam10
  yflow(2)=lam20
  yflow(3)=lam30
  yflow(4)=lam40
  yflow(5)=lam50
  yflow(6)=lam60
  yflow(7)=lam70
  yflow(8)=lam01
  yflow(9)=lam11
  yflow(10)=lam21
  yflow(11)=lam31
  yflow(12)=lam41
  yflow(13)=lam51
  yflow(14)=lam02
  yflow(15)=lam12
  yflow(16)=lam22
  yflow(17)=lam32
  yflow(18)=lam03
  yflow(19)=lam13
  yflow(Nv+1)=lam00
  yflow((Nv+1)+1)=h
  yflow((Nv+1)+2)=h
  yflow((Nv+1)+(Nh+2)+1)=Z_pi
  yflow((Nv+1)+(Nh+2)+2)=Z_K
  yflow((Nv+1)+(Nh+2)+3)=Z_l
  yflow((Nv+1)+(Nh+2)+4)=Z_s
  yflow((Nv+1)+(Nh+2)+5)=Z_A
  yflow((Nv+1)+(Nh+2)+6)=Zphi_p0
  yflow((Nv+1)+(Nh+2)+Nz+1)=jl
  yflow((Nv+1)+(Nh+2)+Nz+2)=js
  yflow((Nv+1)+(Nh+2)+Nz+3)=Sl
  yflow((Nv+1)+(Nh+2)+Nz+4)=Ss
  yflow((Nv+1)+(Nh+2)+Nz+5)=ck
  yflow((Nv+1)+(Nh+2)+Nz+Nck+1)=g
  yflow((Nv+1)+(Nh+2)+Nz+Nck+2)=g   
  !Three gluon interaction
  yflow((Nv+1)+(Nh+2)+Nz+Nck+3)=g   
  !quark-gluon interaction for the strange quark

end
