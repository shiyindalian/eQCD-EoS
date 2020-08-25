subroutine phypoint(yflow,rho0,mPion,mSigma,mf,Vall)
!calculate the physical point rho0

  implicit none
  real(16) yflow(50)
  real(16) rho0,rho
  real(16) mPion2,mPion,mSigma2,mSigma,mf
  real(16) Vall 
  !total effective potential
  integer N_str(5) 
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng
  real(16) lam00,lam10,lam20,lam30,lam40,lam50,lam60,lam70
  real(16) lam01,lam11,lam21,lam31,lam41,lam51
  real(16) lam02,lam12,lam22,lam32
  real(16) lam03,lam13
  real(16) h
  real(16) Zphi,Zpsi,ZA,Zphi_p0
  real(16) kappa
  real(16) g
  integer nn
  parameter(nn=1)
  real(16) x(nn)
  logical check
  external gapEq

  real(16) Z_pi,Z_K,Z_l,Z_s
  real(16) jl,js,Sl,Ss,ck

  common /strucFun/ N_str
  common /Zphi_p0_com/ Zphi_p0

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)

  lam10=yflow(1)
  lam20=yflow(2)
  lam30=yflow(3)
  lam40=yflow(4)
  lam50=yflow(5)
  lam60=yflow(6)
  lam70=yflow(7)
  lam01=yflow(8)
  lam11=yflow(9)
  lam21=yflow(10)
  lam31=yflow(11)
  lam41=yflow(12)
  lam51=yflow(13)
  lam02=yflow(14)
  lam12=yflow(15)
  lam22=yflow(16)
  lam32=yflow(17)
  lam03=yflow(18)
  lam13=yflow(19)
  h=yflow((Nv+1)+1)
  Zphi   =yflow((Nv+1)+(Nh+2)+1)
  Z_K    =yflow((Nv+1)+(Nh+2)+2)
  Zpsi   =yflow((Nv+1)+(Nh+2)+3)
  Z_s    =yflow((Nv+1)+(Nh+2)+4)
  ZA     =yflow((Nv+1)+(Nh+2)+5)
  Zphi_p0=yflow((Nv+1)+(Nh+2)+6)
  jl      =yflow((Nv+1)+(Nh+2)+Nz+1)
  js     =yflow((Nv+1)+(Nh+2)+Nz+2)
  Sl     =yflow((Nv+1)+(Nh+2)+Nz+3)
  Ss     =yflow((Nv+1)+(Nh+2)+Nz+4)
  ck     =yflow((Nv+1)+(Nh+2)+Nz+5)
  g      =yflow((Nv+1)+(Nh+2)+Nz+Nck+1)

  rho0=(Sl**2 + Ss**2)/2.Q+00

  mPion2=lam10 + (ck*Ss)/Sqrt(2.Q+0) + (lam01*(Sl**2 - 2*Ss**2))/6.Q+0

  mPion=sqrt(mPion2)

  mSigma2=(6*lam01*(3*Sl**2 - 2*Ss**2) + lam02*(Sl**3 - 2*Sl*Ss**2)**2 - &
      6*(-6*lam10 - 6*lam20*Sl**2 + 3*Sqrt(2.Q+0)*ck*Ss -               &
         2*lam11*(Sl**4 - 2*Sl**2*Ss**2)))/36.Q+0

  mSigma=sqrt(mSigma2)

  mf=(h**2.Q+00*Sl**2.Q+00)/4.Q+00

  Vall= lam00 -ck*(Sl**2.Q+00*Ss)/2.Q+00/Sqrt(2.Q+00)-jl*Sl-js*Ss

end

