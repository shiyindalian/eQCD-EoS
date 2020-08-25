subroutine dtVdiff1(kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ck,kappak,etaphik,etapsik,rhok,l,lb,dtVd1l,dtVd1lb)
!This subroutine calculate various differential kernal corresponding 0<k<k_UV

  implicit none

  real(8) kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ck,kappak,etaphik,etapsik,rhok
  real(8) l,lb,l_con,lb_con
  real(8) k,lam1,lam2,lam3,lam4,lam5,lam0,h,Zphi,Zpsi,c,kappa,etaphi,etapsi,rho
  real(8) mf2
  real(8) T,mu
  real(8) xff,xfa
  real(8) zb,zf !distinguish the transverse and longituidanl wave function renormalization
  real(8) nf0,nf1,nf2
  real(8) Fnf0,Fnf1,Fnf2
  external Fnf0,Fnf1,Fnf2
  real(8) Nc,Nf
  parameter(Nc=3.,Nf=2.)
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) v3
  parameter(v3=1./(2.*pi**2))
  real(8) nfd1l,nfd1lb
  real(8) nfd1lf,nfd1lbf,nfd1la,nfd1lba
  real(8) dtVd1l,dtVd1lb

  common /Tmu/ T,mu


  l_con=l
  lb_con=lb

  k=kk
  lam1=lam1k
  lam2=lam2k
  lam3=lam3k
  lam4=lam4k
  lam5=lam5k
  lam0=lam0k
  h=hk
  Zphi=Zphik
  Zpsi=Zpsik
  c=ck
  kappa=kappak

  etaphi=etaphik
  etapsi=etapsik

  rho=rhok

  zb=1.
  zf=1.

  mf2=(h**2*rho)/(k**2*Nf)

  xff=-mu + (k*Sqrt(1 + mf2))/zf
  xfa=mu + (k*Sqrt(1 + mf2))/zf

  l=l_con
  lb=lb_con

  nf0=Fnf0(xff,T,l,lb)
  nf1=Fnf1(xff,T,l,lb)
  nf2=Fnf2(xff,T,l,lb)

  nfd1l=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
  nfd1lb=-(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))/2.

  nfd1lf=nfd1l
  nfd1lbf=nfd1lb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l=lb_con
  lb=l_con

  nf0=Fnf0(xfa,T,l,lb)
  nf1=Fnf1(xfa,T,l,lb)
  nf2=Fnf2(xfa,T,l,lb)

  nfd1l=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
  nfd1lb=-(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))/2.

  nfd1la=nfd1l
  nfd1lba=nfd1lb

  l=l_con
  lb=lb_con
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtVd1l=(-2*(1 - etapsi/4.)*k**4*Nc*Nf*(-nfd1lba - nfd1lf)*v3)/(3.*Sqrt(1 + mf2)*zf)
  dtVd1lb=(-2*(1 - etapsi/4.)*k**4*Nc*Nf*(-nfd1la - nfd1lbf)*v3)/(3.*Sqrt(1 + mf2)*zf)

end







