subroutine dtVdiff1(kk,yflowk,etaphik,etapsik,rhok,l,lb,dtVd1l,dtVd1lb)
!This subroutine calculate various differential kernal corresponding 0<k<k_UV

  implicit none

  real(16) kk,yflowk(50),etaphik,etapsik,rhok
  real(16) l,lb,l_con,lb_con
  real(16) k,h,Zphi,Zpsi,jl,js,kappa,etaphi,etapsi,rho
  real(16) Sl,Ss
  real(16) mf2
  real(16) T,mu
  real(16) xff,xfa
  real(16) zb,zf 
  !distinguish the transverse and longituidanl wave function renormalization
  real(16) nf0,nf1,nf2
  real(16) Nc,Nf
  parameter(Nc=3.Q+0,Nf=2.Q+0)
  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) v3
  parameter(v3=1.Q+0/(2.Q+0*pi**2))
  real(16) nfd1l,nfd1lb
  real(16) nfd1lf,nfd1lbf,nfd1la,nfd1lba
  real(16) dtVd1l,dtVd1lb
  integer N_str(5)
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng

  common /Tmu/ T,mu
  common /strucFun/ N_str

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)

  l_con=l
  lb_con=lb

  k=kk
  h=yflowk((Nv+1)+1)
  Zphi=yflowk((Nv+1)+(Nh+2)+1)
  Zpsi=yflowk((Nv+1)+(Nh+2)+3)
  jl=yflowk((Nv+1)+(Nh+2)+Nz+1)
  js=yflowk((Nv+1)+(Nh+2)+Nz+2)
  Sl=yflowk((Nv+1)+(Nh+2)+Nz+3)

  kappa=Sl**2/2.Q+0

  etaphi=etaphik
  etapsi=etapsik

  rho=rhok

  zb=1.Q+0
  zf=1.Q+0

  mf2=(h**2*rho)/(k**2*Nf)

  xff=-mu + (k*Sqrt(1 + mf2))/zf
  xfa=mu + (k*Sqrt(1 + mf2))/zf

  l=l_con
  lb=lb_con

  call Fnf012(xff,T,l,lb,nf0,nf1,nf2)

  nfd1l=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
  nfd1lb=-(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))/2.

  nfd1lf=nfd1l
  nfd1lbf=nfd1lb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l=lb_con
  lb=l_con

  call Fnf012(xfa,T,l,lb,nf0,nf1,nf2)

  nfd1l=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
  nfd1lb=-(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))/2.

  nfd1la=nfd1l
  nfd1lba=nfd1lb

  l=l_con
  lb=lb_con
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtVd1l=(-2*(1 - etapsi/4.Q+0)*k**4*Nc*Nf*(-nfd1lba - nfd1lf)*v3)/(3.Q+0*Sqrt(1 + mf2)*zf)
  dtVd1lb=(-2*(1 - etapsi/4.Q+0)*k**4*Nc*Nf*(-nfd1la - nfd1lbf)*v3)/(3.Q+0*Sqrt(1 + mf2)*zf)

end
