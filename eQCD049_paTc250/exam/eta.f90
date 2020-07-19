subroutine eta(kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ZAk,ck,kappak,gk,g3Ak,etaphik,etapsik,etaAk)

!Calculating anomonous dimension eta

  implicit none

  real(16) kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ZAk,ck,kappak,gk,g3Ak,etaphik,etapsik,etaAk
  real(16) k ! IR cutoff in flow equations
  real(16) lam0,lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(16) h
  real(16) Zphi,Zpsi,ZA
  real(16) c,kappa
  real(16) pi,hc
  parameter(pi=3.141592653589793Q+0)
  parameter(hc=197.33Q+0)
  real(16) etaphi,etapsi,etaA
!meson and quark anomanous dimension
  real(16) Nc,Nf
  parameter(Nc=3.Q+0,Nf=2.Q+0)
  real(16) v3
  parameter(v3=1.Q+0/(2.Q+0*pi**2))
  real(16) rho !phi_a**2/2
  real(16) etaAYM_T_fun,etaA_unqQCD_T0_fun,etac_YM_T0_fun
  external etaAYM_T_fun,etaA_unqQCD_T0_fun,etac_YM_T0_fun
  real(16) zb,zf !distinguish the transverse and longituidanl wave function renormalization
  real(16) p0,p0c !temporal compontent of external momentum
  real(16) T,mu
  real(16) mu0
  real(16) l,lb !polyakov loop
  real(16) k_UV,k_IR,t_UV,t_IR
  real(16) mp2,ms2,mf2,mp2d1rho,mp2d2rho,mp2d3rho,mp2d4rho,mp2d5rho,ms2d1rho,ms2d2rho,ms2d3rho,ms2d4rho,ms2d5rho,mf2d1rho
  real(16) dr0dtV,dr1dtV,dr2dtV,dr3dtV,dr4dtV,dr5dtV
  real(16) dlam0dt,dlam1dt,dlam2dt,dlam3dt,dlam4dt,dlam5dt,dlam6dt,dlam7dt
  real(16) dth,dhdt
  real(16) dZphidt,dZpsidt,dZAdt,dcdt,dkappadt
  real(16) l_com,lb_com

  real(16) nffFd0x,nffFd1x,nffFd2x,nffFd3x,nffFd4x,nffFd5x,nfaFd0x,nfaFd1x,nfaFd2x,nfaFd3x,nfaFd4x,nfaFd5x
  real(16) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  real(16) nffF,nfaF,nffFd1,nfaFd1,nffFd2,nfaFd2,nffFd3,nfaFd3,nffFd4,nfaFd4,nffFd5,nfaFd5
  real(16) nbBd0x,nbBd1x,nbBd2x,nbBd3x,nbBd4x,nbBd5x
  real(16) nbB,nbBd1,nbBd2,nbBd3,nbBd4,nbBd5
  real(16) finvEB,finvEBd1,finvEBd2
  real(16) nbBa,nbBb,nbBad1,nbBbd1
  real(16) finvEBa,finvEBb,finvEBad1,finvEBbd1


  real(16) nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,nbSigma,nbd1xSigma,nbd2xSigma,nbd3xSigma,&
          nbd4xSigma,nbd5xSigma
  real(16) nbGluon,nbd1xGluon,nbd2xGluon,nbd3xGluon,nbd4xGluon,nbd5xGluon

  real(16) b2f1aP,b1f2P,b2f2aP,b1f3P,b3f1aP
  real(16) b2f1aS,b1f2S,b2f2aS,b1f3S
  real(16) b2f1aA,b1f2A,b2f2aA,b1f3A,b3f1aA,b3f2aA,b2f3aA
  real(16) b1f1A

  real(16) b2b2PS
  real(16) mb2
  real(16) mb2a,mb2b
  real(16) f2a,f3a,f3aEtaphi
  real(16) f2aT,f3aT
  real(16) b2f1a,b1f2,b2f2a,b1f3,b3f1a
  real(16) b2b2

  real(16) g,g3A,gccA,gAAA
  real(16) etaAYM_T,etaAYM_T0
  real(16) v4,C2Nc
  real(16) detaA,detaA_mass,detaA_mass_T,etaAM
  integer n_x,i_x
!  parameter(n_x=64)
  parameter(n_x=128)
  real(16) w_x(n_x),y_x(n_x) !Guass integral
  integer n_cth,i_cth
!  parameter(n_cth=128)
  parameter(n_cth=128)
  real(16) w_cth(n_cth),y_cth(n_cth) !Guass integral
  real(16) etaA_QL,cthinte_QL

  real(16) etaA_unq_T0
  real(16) etaAT0,etacT0
  real(16) fGL0,fGL1,fGhL0,fGhL1,fQL0,fQL1
  real(16) q_ti,p,p_ti,costhe0,costhe1

  real(16) r_4d_to_3d,k_4d

  common /Tmu/ T,mu
  common /mu0p0_com/ mu0,p0,p0c
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /polyakov_com/ l_com,lb_com
  common /nffFdx_com/ nffFd0x,nffFd1x,nffFd2x,nffFd3x,nffFd4x,nffFd5x,nfaFd0x,nfaFd1x,nfaFd2x,nfaFd3x,nfaFd4x,nfaFd5x
  common /nbBdx_com/ nbBd0x,nbBd1x,nbBd2x,nbBd3x,nbBd4x,nbBd5x
  common /nbBd_com/ nbB,nbBd1,nbBd2,nbBd3,nbBd4,nbBd5
  common /finvEBd_com/ finvEB,finvEBd1,finvEBd2
  common /nbBad_com/ nbBa,nbBb,nbBad1,nbBbd1
  common /finvEBad_com/ finvEBa,finvEBb,finvEBad1,finvEBbd1
  common /nbGluon_com/ nbGluon,nbd1xGluon,nbd2xGluon,nbd3xGluon,nbd4xGluon,nbd5xGluon
  common /nbPS_com/ nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,&
                    nbSigma,nbd1xSigma,nbd2xSigma,nbd3xSigma,nbd4xSigma,nbd5xSigma
  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa

  common /BF_mes_thr_com/ b1f2P,b1f2S,b1f3P,b1f3S,b2f1aP,b2f1aS,b2f2aP,b2f2aS,b3f1aP
  common /BF_thr_com/ b2f1a,b1f2,b2f2a,b1f3,b3f1a
  common /BF_A_thr_com/ b1f1A,b2f1aA,b1f2A,b2f2aA,b1f3A,b3f1aA,b3f2aA,b2f3aA
  common /r_4d_to_3d_com/r_4d_to_3d

  k=kk
  lam1=lam1k
  lam2=lam2k
  lam3=lam3k
  lam4=lam4k
  lam5=lam5k
  lam6=0.Q+0
  lam7=0.Q+0
  lam0=lam0k
  h=hk
  Zphi=Zphik
  Zpsi=Zpsik
  ZA=ZAk
  c=ck
  kappa=kappak
  g=gk
  g3A=g3Ak

  rho=kappa 
  !calculations are performed at expansion point kappa

  zb=1.
  zf=1.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!mass and their derivatives
  mp2=lam1/k**2
  ms2=(lam1 + 2*lam2*rho)/k**2
  mf2=(h**2*rho)/(k**2*Nf)
  mp2d1rho=lam2/k**2
  mp2d2rho=lam3/k**2
  mp2d3rho=lam4/k**2
  mp2d4rho=lam5/k**2
  mp2d5rho=lam6/k**2
  ms2d1rho=(3*lam2 + 2*lam3*rho)/k**2
  ms2d2rho=(5*lam3 + 2*lam4*rho)/k**2
  ms2d3rho=(7*lam4 + 2*lam5*rho)/k**2
  ms2d4rho=(9*lam5 + 2*lam6*rho)/k**2
  ms2d5rho=(11*lam6 + 2*lam7*rho)/k**2
  mf2d1rho=h**2/(k**2*Nf)

  if(abs(ms2-mp2)<1.Q-7*(ms2+mp2)/2.Q+0)then
    ms2=mp2+1.Q-4*(ms2+mp2)/2.Q+0
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  l=l_com
  lb=lb_com

  call nfdx_com(mf2,T,mu,l,lb,k)
  nff=nffFd0x
  nfd1xf=nffFd1x
  nfd2xf=nffFd2x
  nfd3xf=nffFd3x
  nfd4xf=nffFd4x
  nfd5xf=nffFd5x

  nfa=nfaFd0x
  nfd1xa=nfaFd1x
  nfd2xa=nfaFd2x
  nfd3xa=nfaFd3x
  nfd4xa=nfaFd4x
  nfd5xa=nfaFd5x

  call F_thr(mf2,k,f2a,f3a)
  call FT_thr(mf2,k,f2aT,f3aT)

  mb2=mp2
  call nbdx_com(mb2,T,k)
  nbPion=nbBd0x
  nbd1xPion=nbBd1x
  nbd2xPion=nbBd2x
  nbd3xPion=nbBd3x
  nbd4xPion=nbBd4x
  nbd5xPion=nbBd5x

  nbBa=nbB
  nbBad1=nbBd1

  finvEBa=finvEB
  finvEBad1=finvEBd1

  mb2=ms2
  call nbdx_com(mb2,T,k)
  nbSigma=nbBd0x
  nbd1xSigma=nbBd1x
  nbd2xSigma=nbBd2x
  nbd3xSigma=nbBd3x
  nbd4xSigma=nbBd4x
  nbd5xSigma=nbBd5x

  nbBb=nbB
  nbBbd1=nbBd1

  finvEBb=finvEB
  finvEBbd1=finvEBd1

  call BF_mes_thr(mp2,ms2,mf2,T,k)

  mb2a=mp2
  mb2b=ms2
  call BB_thr(mb2a,mb2b,k,b2b2)
  b2b2PS=b2b2

  mb2=0.
  call nbdx_com(mb2,T,k)
  nbGluon=nbBd0x
  nbd1xGluon=nbBd1x
  nbd2xGluon=nbBd2x
  nbd3xGluon=nbBd3x
  nbd4xGluon=nbBd4x
  nbd5xGluon=nbBd5x

  call BF_A_thr(mf2,T,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  v4=1./32./pi**2
  C2Nc=(Nc**2-1.)/(2.*Nc)

  f3aEtaphi=f3a

  k_4d=r_4d_to_3d*k

  etaAT0=etaA_unqQCD_T0_fun(k_4d)
  etacT0=etac_YM_T0_fun(k_4d)

  etaAYM_T=0.Q+0
  etaAYM_T0=0.Q+0

  gccA=g

  mb2=0.Q+0
  call nbdx_com(mb2,T,k)
  nbGluon=nbBd0x
  nbd1xGluon=nbBd1x
  nbd2xGluon=nbBd2x
  nbd3xGluon=nbBd3x
  nbd4xGluon=nbBd4x
  nbd5xGluon=nbBd5x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  etaphi=(v3*(((f2a - 2*f3a)*h**2*Nc*(4*b2f1aS*h**2 +                            &
           C2Nc*(9*b1f1A - 18*b1f2A + 8*b2f1aA - 2*b2f1aA*(etaAT0+(etaAYM_T-etaAYM_T0)))*g**2*  &
            Nf + 4*b2f1aP*h**2*(-1 + Nf**2))*v3)/Nf -                            &
      6*(3*f2a*h**2*Nc - 8*f3a*h**2*Nc - (4*b2b2PS*lam2**2*rho)/k**2)*           &
       (1 + ((b1f1A - 2*b1f2A)*C2Nc*g**2*v3)/4.)))/                              &
  (18.*(1 + ((b1f1A - 2*b1f2A)*C2Nc*g**2*v3)/4. +                                &
      ((f2a - 2*f3a)*h**4*Nc*(b2f1aS + b2f1aP*(-1 + Nf**2))*v3**2)/(18.*Nf)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  etapsi=-((b2f1aS*(-4 + etaphi)*h**2 +                                          &
       C2Nc*(-9*b1f1A + 18*b1f2A + 2*b2f1aA*(-4 + (etaAT0+(etaAYM_T-etaAYM_T0))))*g**2*Nf +    &
       b2f1aP*(-4 + etaphi)*h**2*(-1 + Nf**2))*v3)/                              &
  (3.*Nf*(4 + b1f1A*C2Nc*g**2*v3 - 2*b1f2A*C2Nc*g**2*v3))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  etaphik=etaphi
  etapsik=etapsi
  etaAk=0.Q+0

end
