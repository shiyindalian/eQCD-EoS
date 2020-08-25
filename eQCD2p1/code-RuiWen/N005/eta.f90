subroutine eta(kk,yflowk,etaphik,etapsik,etaAk)

!Calculating anomonous dimension eta

  implicit none

  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) Nc,Nf
  parameter(Nc=3.Q+0,Nf=2.Q+0)
  real(16) v3
  parameter(v3=1.Q+0/(2.Q+0*pi**2))
  real(16) kk,yflowk(50),etaphik,etapsik,etaAk
  real(16) k 
  ! IR cutoff in flow equations
  integer N_str(5)
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng
  real(16) lam00,lam10,lam20,lam30,lam40,lam50,lam60,lam70
  real(16) lam01,lam11,lam21,lam31,lam41,lam51
  real(16) lam02,lam12,lam22,lam32
  real(16) lam03,lam13
  real(16) h
  real(16) Zphi,Zpsi,ZA
  real(16) c,kappa1,kappa2
  real(16) Sl,Ss,jl,js
  real(16) etaphi,etapsi,etaA
  !meson and quark anomanous dimension
  real(16) rho,rho2
  real(16) etaA_unqQCD_T0_fun,etac_YM_T0_fun
  external etaA_unqQCD_T0_fun,etac_YM_T0_fun
  real(16) zb,zf 
  !distinguish the transverse and longituidanl wave function renormalization
  real(16) p0,p0c 
  !temporal compontent of external momentum
  real(16) T,mu
  real(16) mu0
  real(16) l,lb 
  !polyakov loop
  real(16) mp2,ms2,mf2
  real(16) dmb2drho(2,12)
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
  real(16) f2a,f3a
  real(16) f2aT,f3aT
  real(16) b2f1a,b1f2,b2f2a,b1f3,b3f1a
  real(16) b2b2

  real(16) g
  real(16) etaAYM_T,etaAYM_T0
  real(16) v4,C2Nc
  real(16) detaA,detaA_mass,detaA_mass_T,etaAM
  real(16) etaA_QL,cthinte_QL

  real(16) etaAT0,etacT0
  real(16) fGL0,fGL1,fGhL0,fGhL1,fQL0,fQL1

  real(16) r_4d_to_3d,k_4d

  common /Tmu/ T,mu
  common /mu0p0_com/ mu0,p0,p0c
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

  common /strucFun/ N_str

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)

  k=kk

  lam10=yflowk(1)
  lam20=yflowk(2)
  lam30=yflowk(3)
  lam40=yflowk(4)
  lam50=yflowk(5)
  lam60=yflowk(6)
  lam70=yflowk(7)
  lam01=yflowk(8)
  lam11=yflowk(9)
  lam21=yflowk(10)
  lam31=yflowk(11)
  lam41=yflowk(12)
  lam51=yflowk(13)
  lam02=yflowk(14)
  lam12=yflowk(15)
  lam22=yflowk(16)
  lam32=yflowk(17)
  lam03=yflowk(18)
  lam13=yflowk(19)
  lam00=yflowk(Nv+1)
  h=yflowk((Nv+1)+1)
  Zphi=yflowk((Nv+1)+(Nh+2)+1)
  Zpsi=yflowk((Nv+1)+(Nh+2)+3)
  ZA  =yflowk((Nv+1)+(Nh+2)+5)
  c   =yflowk((Nv+1)+(Nh+2)+Nz+1)
  Sl  =yflowk((Nv+1)+(Nh+2)+Nz+3)
  g   =yflowk((Nv+1)+(Nh+2)+Nz+Nck+1)

  kappa1=(Sl**2 + Ss**2)/2.Q+00
  kappa2=(Sl**2 - 2*Ss**2)**2/24.Q+00
  rho=kappa1
  rho2=kappa2
  !calculations are performed at expansion point kappa

  zb=1.Q+0
  zf=1.Q+0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!mass and their derivatives
  call mass2(yflowk,dmb2drho) 
  dmb2drho=dmb2drho/k**2
  mp2=dmb2drho(1,1)
  ms2=dmb2drho(2,1)
  mf2=(h**2.Q+00*Sl**2.Q+00)/(4.Q+00*k**2)

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

  mb2=0.Q+0
  call nbdx_com(mb2,T,k)
  nbGluon=nbBd0x
  nbd1xGluon=nbBd1x
  nbd2xGluon=nbBd2x
  nbd3xGluon=nbBd3x
  nbd4xGluon=nbBd4x
  nbd5xGluon=nbBd5x

  call BF_A_thr(mf2,T,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  v4=1.Q+0/32.Q+0/pi**2
  C2Nc=(Nc**2-1.Q+0)/(2.Q+0*Nc)

  k_4d=r_4d_to_3d*k

  etaAT0=etaA_unqQCD_T0_fun(k_4d)
  etacT0=etac_YM_T0_fun(k_4d)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  etaphi=(v3*(((f2a - 2*f3a)*h**2*Nc*(4*b2f1aS*h**2 +                            &
  C2Nc*(9*b1f1A - 18*b1f2A + 8*b2f1aA - 2*b2f1aA*(etaAT0))*g**2*  &
            Nf + 4*b2f1aP*h**2*(-1 + Nf**2))*v3)/Nf -                            &
      6*(3*f2a*h**2*Nc - 8*f3a*h**2*Nc - (4*b2b2PS*lam20**2*rho)/k**2)*           &
       (1 + ((b1f1A - 2*b1f2A)*C2Nc*g**2*v3)/4.)))/                              &
  (18.*(1 + ((b1f1A - 2*b1f2A)*C2Nc*g**2*v3)/4. +                                &
      ((f2a - 2*f3a)*h**4*Nc*(b2f1aS + b2f1aP*(-1 + Nf**2))*v3**2)/(18.*Nf)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  etapsi=-((b2f1aS*(-4 + etaphi)*h**2 +                                          &
       C2Nc*(-9*b1f1A + 18*b1f2A + 2*b2f1aA*(-4 + (etaAT0)))*g**2*Nf +    &
       b2f1aP*(-4 + etaphi)*h**2*(-1 + Nf**2))*v3)/                              &
  (3.*Nf*(4 + b1f1A*C2Nc*g**2*v3 - 2*b1f2A*C2Nc*g**2*v3))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  etaphik=etaphi
  etapsik=etapsi
  etaAk=0.Q+0

end
