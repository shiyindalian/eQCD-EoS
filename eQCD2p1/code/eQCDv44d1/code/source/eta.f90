subroutine eta(kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ZAk,ck,kappak,gk,g3Ak,etaphik,etapsik,etaAk)

!Calculating anomonous dimension eta

  implicit none

  real(8) kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ZAk,ck,kappak,gk,g3Ak,etaphik,etapsik,etaAk
  real(8) k ! IR cutoff in flow equations
  real(8) lam0,lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(8) h
  real(8) Zphi,Zpsi,ZA
  real(8) c,kappa
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) etaphi,etapsi,etaA
!meson and quark anomanous dimension
  real(8) Nc,Nf
  parameter(Nc=3.,Nf=2.)
  real(8) v3
  parameter(v3=1./(2.*pi**2))
  real(8) rho !phi_a**2/2
  real(8) Fnb,Fnf0,Fnf1,Fnf2,etaAYM_T_fun,etaA_unqQCD_T0_fun,etac_YM_T0_fun,g3A_T0_fun
  external Fnb,Fnf0,Fnf1,Fnf2,etaAYM_T_fun,etaA_unqQCD_T0_fun,etac_YM_T0_fun,g3A_T0_fun
  real(8) zb,zf !distinguish the transverse and longituidanl wave function renormalization
  real(8) p0,p0c !temporal compontent of external momentum
  real(8) T,mu
  real(8) mu0
  real(8) l,lb !polyakov loop
  real(8) k_UV,k_IR,t_UV,t_IR
  real(8) mp2,ms2,mf2,mp2d1rho,mp2d2rho,mp2d3rho,mp2d4rho,mp2d5rho,ms2d1rho,ms2d2rho,ms2d3rho,ms2d4rho,ms2d5rho,mf2d1rho
  real(8) dr0dtV,dr1dtV,dr2dtV,dr3dtV,dr4dtV,dr5dtV
  real(8) dlam0dt,dlam1dt,dlam2dt,dlam3dt,dlam4dt,dlam5dt,dlam6dt,dlam7dt
  real(8) dth,dhdt
  real(8) dZphidt,dZpsidt,dZAdt,dcdt,dkappadt
  real(8) l_com,lb_com
  complex(8) B1B1F3pPSC,B1B1F3aPSC,B2F3aPC,B2F3aSC,B2F3pPC,B2F3pSC,B2B1F2aPSC,B2B1F2aSPC,B2B1F2pPSC, &
             B2B1F2pSPC,B3F2aPC,B3F2aSC,B3F2pPC,B3F2pSC
  real(8) B1B1F3pPSI,B1B1F3aPSI,B2F3aPI,B2F3aSI,B2F3pPI,B2F3pSI,B2B1F2aPSI,B2B1F2aSPI,B2B1F2pPSI, &
          B2B1F2pSPI,B3F2aPI,B3F2aSI,B3F2pPI,B3F2pSI
  real(8) B1B1F3pPS,B1B1F3aPS,B2F3aP,B2F3aS,B2F3pP,B2F3pS,B2B1F2aPS,B2B1F2aSP,B2B1F2pPS, &
          B2B1F2pSP,B3F2aP,B3F2aS,B3F2pP,B3F2pS

  real(8) nffFd0x,nffFd1x,nffFd2x,nffFd3x,nffFd4x,nffFd5x,nfaFd0x,nfaFd1x,nfaFd2x,nfaFd3x,nfaFd4x,nfaFd5x
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  real(8) nffF,nfaF,nffFd1,nfaFd1,nffFd2,nfaFd2,nffFd3,nfaFd3,nffFd4,nfaFd4,nffFd5,nfaFd5
  real(8) nbBd0x,nbBd1x,nbBd2x,nbBd3x,nbBd4x,nbBd5x
  real(8) nbB,nbBd1,nbBd2,nbBd3,nbBd4,nbBd5
  real(8) finvEB,finvEBd1,finvEBd2
  real(8) nbBa,nbBb,nbBad1,nbBbd1
  real(8) finvEBa,finvEBb,finvEBad1,finvEBbd1


  real(8) nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,nbSigma,nbd1xSigma,nbd2xSigma,nbd3xSigma,&
          nbd4xSigma,nbd5xSigma
  real(8) nbGluon,nbd1xGluon,nbd2xGluon,nbd3xGluon,nbd4xGluon,nbd5xGluon
  real(8) nfGhost,nfd1xGhost,nfd2xGhost,nfd3xGhost,nfd4xGhost

  real(8) b2f1aP,b1f2P,b2f2aP,b1f3P,b3f1aP
  real(8) b2f1aS,b1f2S,b2f2aS,b1f3S,b3f1aS
  real(8) b2f1aA,b1f2A,b2f2aA,b1f3A,b3f1aA,b3f2aA,b2f3aA
  real(8) b1f1A

  real(8) b2b2PS
  real(8) mb2
  real(8) mb2a,mb2b
  real(8) f2a,f3a,f4a,f3aEtaphi
  real(8) F1F1,F2F1
  real(8) B3A,B4A,B4cA,B5A,B2B1A,B2A,F3Gh
  real(8) f2aT,f3aT,f4aT
  real(8) b2f1a,b1f2,b2f2a,b1f3,b3f1a
  real(8) b2b2

  real(8) g,g3A,gccA,gAAA
  real(8) etaAYM_T,etaAYM_T0
  real(8) v4,C2Nc
  real(8) etaphi_cut,n_etaphi_cut
  real(8) detaA,detaA_mass,detaA_mass_T,etaAM

  integer n_x,i_x
!  parameter(n_x=64)
  parameter(n_x=64)
  real(8) w_x(n_x),y_x(n_x) !Guass integral
  integer n_cth,i_cth
!  parameter(n_cth=128)
  parameter(n_cth=64)
  real(8) w_cth(n_cth),y_cth(n_cth) !Guass integral
  real(8) etaA_QL,cthinte_QL
  real(8) x_q,xprime,q,costhe,rF_plus_1
  real(8) etaA_GL,etaA_GhL,cthinte_GL,cthinte_GhL,xcos

  real(8) etaA_unq_T0
  real(8) etaAT0,etacT0
  real(8) fGL0,fGL1,fGhL0,fGhL1,fQL0,fQL1
  real(8) err1,err2,err3,err4
  real(8) x_low,x_high
  real(8) q_ti,p,p_ti,costhe0,costhe1

  real(8) r_4d_to_3d,k_4d


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
  common /nfGhost_com/ nfGhost,nfd1xGhost,nfd2xGhost,nfd3xGhost,nfd4xGhost
  common /nbPS_com/ nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,&
                    nbSigma,nbd1xSigma,nbd2xSigma,nbd3xSigma,nbd4xSigma,nbd5xSigma
  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa

  common /F_thr_com/ f2a,f3a,f4a
  common /FF_thr_com/ F1F1,F2F1
  common /FT_thr_com/ f2aT,f3aT,f4aT
  common /BF_mes_thr_com/ b1f2P,b1f2S,b1f3P,b1f3S,b2f1aP,b2f1aS,b2f2aP,b2f2aS,b3f1aP,b3f1aS,b2f3aP,b2f3aS,b3f2aP,b3f2aS
  common /BF_thr_com/ b2f1a,b1f2,b2f2a,b1f3,b3f1a
  common /BB_thr_com/ b2b2
  common /BF_A_thr_com/ b1f1A,b2f1aA,b1f2A,b2f2aA,b1f3A,b3f1aA,b3f2aA,b2f3aA
  common /BA_thr_com/ B3A,B4A,B4cA,B5A,B2B1A,B2A,F3Gh
  common /etaphi_cut_com/etaphi_cut,n_etaphi_cut
  common /err1_com/err1,err2,err3,err4,x_low,x_high

  common /r_4d_to_3d_com/r_4d_to_3d


  k=kk
  lam1=lam1k
  lam2=lam2k
  lam3=lam3k
  lam4=lam4k
  lam5=lam5k
  lam6=0.
  lam7=0.
  lam0=lam0k
  h=hk
  Zphi=Zphik
  Zpsi=Zpsik
  ZA=ZAk
  c=ck
  kappa=kappak
  g=gk
  g3A=g3Ak



  rho=kappa !calculations are performed at expansion point kappa

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


  if(abs(ms2-mp2)<1.d-7*(ms2+mp2)/2.)then
    ms2=mp2+1.d-4*(ms2+mp2)/2.
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

  call F_thr(mf2,k)
  call FT_thr(mf2,k)

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
  call BB_thr(mb2a,mb2b,k)
  b2b2PS=b2b2

  mb2=0.
  call nbdx_com(mb2,T,k)
  nbGluon=nbBd0x
  nbd1xGluon=nbBd1x
  nbd2xGluon=nbBd2x
  nbd3xGluon=nbBd3x
  nbd4xGluon=nbBd4x
  nbd5xGluon=nbBd5x

  call nfGhostdx_com(T,k)

!  call fabmf_com(mb2,mf2,T,k)
!  call fm_com(mb2,mf2,T,k)
!  call thr_A(mb2,mf2,T,k)
  call BF_A_thr(mf2,T,k)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  v4=1./32./pi**2
  C2Nc=(Nc**2-1.)/(2.*Nc)

!  f2a=F1F1
!  f3a=F1F2
!  f3aEtaphi=3./(4.*(1.+mf2))*f2a
  f3aEtaphi=f3a


  k_4d=r_4d_to_3d*k

  etaAT0=etaA_unqQCD_T0_fun(k_4d)
  etacT0=etac_YM_T0_fun(k_4d)
!  etaAYM_T=etaAYM_T_fun(T,k_4d)
!  etaAYM_T0=etaAYM_T_fun(1.d0/hc,k_4d)

  etaAYM_T=0.
  etaAYM_T0=0.


!  gccA=g3A
  gccA=g


! quark loop is calculated in what follows
!  call gauleg(0.d0,1.d0, y_x, w_x, n_x)
!  call gauleg(-1.d0,1.d0, y_cth, w_cth, n_cth)


!  fQL0=0.
!  fQL1=0.
!goto 310
!  do i_x=1,n_x
!    x_q=y_x(i_x)
!    q=k*Sqrt(x_q)

!    cthinte_QL=0.
!    do i_cth=1,n_cth
!      costhe=y_cth(i_cth)
!      call FF_thr(mf2,T,mu,l,lb,k,q,costhe)
!      xprime=x_q+1.-2.*Sqrt(x_q)*costhe
!      if(xprime<1.d0)then
!        rF_plus_1=1./Sqrt(xprime)
!      else
!        rF_plus_1=1.
!      end if

!      cthinte_QL=cthinte_QL+w_cth(i_cth)*((F1F1-F2F1) &
!                       +(Sqrt(x_q)*costhe**2-costhe)*rF_plus_1*(F2F1-F1F1/2.))
!    end do
!    fQL0=fQL0+w_x(i_x)*Sqrt(x_q)*cthinte_QL
!    fQL1=fQL1+w_x(i_x)*x_q*cthinte_QL
!  end do
!310 continue

!  gluon and ghost loops are calculated in what follows

  mb2=0.
  call nbdx_com(mb2,T,k)
  nbGluon=nbBd0x
  nbd1xGluon=nbBd1x
  nbd2xGluon=nbBd2x
  nbd3xGluon=nbBd3x
  nbd4xGluon=nbBd4x
  nbd5xGluon=nbBd5x


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  p=k              !external momentum for the gluon propagator
!  p_ti=p/k

!  fGL0=0.
!  fGL1=0.
!  fGhL0=0.
!  fGhL1=0.
!  call gauleg(x_low,x_high,y_x,w_x,n_x)
!  call gauleg(-1.d0,1.d0,y_cth,w_cth,n_cth)
!  do i_x=1,n_x
!    x_q=y_x(i_x)
!    q=k*Sqrt(x_q)

!    q_ti=Sqrt(x_q)

!    costhe0=(q_ti**2+p_ti**2-1.)/(2.*q_ti*p_ti)    !qmp_ti=1
!    costhe1=p_ti/(2.*q_ti)                         !q_ti=qmp_ti

!    cthinte_GL=0.
!    cthinte_GhL=0.
!    do i_cth=1,n_cth
!      costhe=y_cth(i_cth)


!      if(abs(costhe-costhe0)<err1.or.abs(costhe-costhe1)<err1)then
!        cthinte_GL=cthinte_GL+w_cth(i_cth)*(0.)
!        cthinte_GhL=cthinte_GhL+w_cth(i_cth)*(0.)
!      else
!        call BA_thr(T,k,q,p,costhe)
!        xcos=x_q*(1.-costhe**2)

!        cthinte_GL=cthinte_GL+w_cth(i_cth)*(xcos*p_ti**4*B5A                      &
!                   +(xcos*(-p_ti**2-2.*Sqrt(x_q)*p_ti*costhe)-8.*x_q*p_ti**2*costhe**2)*B4A   &
!                   +(xcos*(-3.*p_ti**2+2.*Sqrt(x_q)*p_ti*costhe)      &
!                       -8.*(Sqrt(x_q)*p_ti*costhe-p_ti**2)**2)*B4cA   &
!                   +(12.*xcos+16.*p_ti**2)*B3A-(xcos*B2B1A+4.*B2A))

!        cthinte_GhL=cthinte_GhL+w_cth(i_cth)*xcos*F3Gh
!      end if
!    end do

!    fGL0=fGL0+w_x(i_x)*Sqrt(x_q)*cthinte_GL
!    fGL1=fGL1+w_x(i_x)*(Sqrt(x_q))**3*cthinte_GL

!    fGhL0=fGhL0+w_x(i_x)*Sqrt(x_q)*cthinte_GhL
!    fGhL1=fGhL1+w_x(i_x)*(Sqrt(x_q))**3*cthinte_GhL

!  end do
!  fGL0=fGL0/p_ti**2
!  fGL1=fGL1/p_ti**2
!  fGhL0=fGhL0/p_ti**2
!  fGhL1=fGhL1/p_ti**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  etaphi=(v3*(((f2a - 2*f3a)*h**2*Nc*(4*b2f1aS*h**2 +                            &
           C2Nc*(9*b1f1A - 18*b1f2A + 8*b2f1aA - 2*b2f1aA*(etaAT0+(etaAYM_T-etaAYM_T0)))*g**2*  &
            Nf + 4*b2f1aP*h**2*(-1 + Nf**2))*v3)/Nf -                            &
      6*(3*f2a*h**2*Nc - 8*f3a*h**2*Nc - (4*b2b2PS*lam2**2*rho)/k**2)*           &
       (1 + ((b1f1A - 2*b1f2A)*C2Nc*g**2*v3)/4.)))/                              &
  (18.*(1 + ((b1f1A - 2*b1f2A)*C2Nc*g**2*v3)/4. +                                &
      ((f2a - 2*f3a)*h**4*Nc*(b2f1aS + b2f1aP*(-1 + Nf**2))*v3**2)/(18.*Nf)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  etaphi_cut=1.
!  n_etaphi_cut=1.
!  if(etaphi<0.)then
!    etaphi=etaphi*exp(-(abs(etaphi/etaphi_cut))**n_etaphi_cut)
!  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  etapsi=-((b2f1aS*(-4 + etaphi)*h**2 +                                          &
       C2Nc*(-9*b1f1A + 18*b1f2A + 2*b2f1aA*(-4 + (etaAT0+(etaAYM_T-etaAYM_T0))))*g**2*Nf +    &
       b2f1aP*(-4 + etaphi)*h**2*(-1 + Nf**2))*v3)/                              &
  (3.*Nf*(4 + b1f1A*C2Nc*g**2*v3 - 2*b1f2A*C2Nc*g**2*v3))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  etaphik=etaphi
  etapsik=etapsi
  etaAk=0.

end
