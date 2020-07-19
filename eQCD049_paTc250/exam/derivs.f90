subroutine derivs(x,y,dydx)
!Calculating the right hand side of differential equations

  implicit none

  integer NMAX 
  !maximal number of differential equations
  parameter(NMAX=50)
  real(16) x,y(NMAX),dydx(NMAX)
  integer N_str(5) 
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng
  real(16) k 
  ! IR cutoff in flow equations
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
  parameter(v3=1.Q+00/(2.Q+0*pi**2))
  real(16) rho 
  !phi_a**2/2
  real(16) etaAYM_T_fun,etaA_unqQCD_T0_fun,etac_YM_T0_fun
  external etaAYM_T_fun,etaA_unqQCD_T0_fun,etac_YM_T0_fun
  real(16) zb,zf 
  !distinguish the transverse and longituidanl wave function renormalization
  real(16) p0,p0c 
  !temporal compontent of external momentum
  real(16) T,mu
  real(16) mu0
  real(16) l,lb 
  !polyakov loop
  real(16) k_UV,k_IR,t_UV,t_IR
  real(16) mp2,ms2,mf2,mp2d1rho,mp2d2rho,mp2d3rho,mp2d4rho,mp2d5rho,ms2d1rho,ms2d2rho,ms2d3rho,ms2d4rho,ms2d5rho,mf2d1rho
  real(16) dr0dtV,dr1dtV,dr2dtV,dr3dtV,dr4dtV,dr5dtV
  real(16) dlam0dt,dlam1dt,dlam2dt,dlam3dt,dlam4dt,dlam5dt,dlam6dt,dlam7dt
  real(16) dth,dhdt
  real(16) dZphidt,dZpsidt,dZAdt,dcdt,dkappadt
  real(16) l_com,lb_com
  complex(16) B2B1F1aPSC,B2B1F1aSPC,B1B1F2aPSC,B1B1F3aPSC,B2B1F2aPSC,B2B1F2aSPC,B2F3aPC,B2F3aSC,B3F2aPC,B3F2aSC
  real(16) B2B1F1aPS,B2B1F1aSP,B1B1F2aPS,B1B1F3aPS,B2B1F2aPS,B2B1F2aSP,B2F3aP,B2F3aS,B3F2aP,B3F2aS

  real(16) dtA,dtlambdaSigmaPion
  real(16) nffFd0x,nffFd1x,nffFd2x,nffFd3x,nffFd4x,nffFd5x,nfaFd0x,nfaFd1x,nfaFd2x,nfaFd3x,nfaFd4x,nfaFd5x
  real(16) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
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
  real(16) b1f1AT,b2f1aAT,b1f2AT,b2f2aAT,b1f3AT,b3f1aAT
  real(16) b1f1A

  real(16) L11P,L11S,L11A
  real(16) b2b2PS
  real(16) mb2
  real(16) mb2a,mb2b
  real(16) f2a,f3a,f3aEtaphi
  real(16) F1F1,F2F1,F1F1phi,F2F1phi
  real(16) f2aT,f3aT
  real(16) b2f1a,b1f2,b2f2a,b1f3,b3f1a
  real(16) b2b2

  real(16) dtg,dtlam4,dtg3A

  real(16) g,gAAA
  real(16) etaAYM_T,etaAYM_T0
  real(16) ag,bg,dg,fginf,dtfginf
  real(16) v4,C2Nc,ashift,dashift
  real(16) N21mP,N21mS,N21gA,N12A,L12A,L111PS,L11gA,dthA
  real(16) dtlam4_gluon
  real(16) dtg_gluon,dtg_meson,dtg_gluon_T
  integer k_num,k_num1
  real(16) k_value,k_value1
  real(16) dtg_gluon_4d,dtg_meson_4d,dtg_non,dtlam4_gluon_4d,dtlam4_meson_4d,dtlam4_4d,dtlam4_meson_false,&
          dth_yukawa,dth_4fermi,etaphi_mesPqua,etaphi_qua,etaphi_f2a,etaphi_f3a,dtg_etaA,dtg_etapsi
  real(16) g_max,alphas_max,k_g_max,ZAm1_max,k_ZAm1_max,g3A_max,alphas3A_max,k_g3A_max
  real(16) etaphi_quark,etaphi_meson
  real(16) test1

  integer n_x,i_x
  parameter(n_x=64)
  real(16) w_x(n_x),y_x(n_x) 
  !Guass integral
  integer n_cth,i_cth
  parameter(n_cth=64)
  real(16) w_cth(n_cth),y_cth(n_cth) 
  !Guass integral
  real(16) etaA_QL,cthinte_QL
  real(16) x_q,xprime,q,costhe,rF_plus_1

  integer k_same_common

  real(16) etaA_unq_T0
  real(16) l_bar
  real(16) etaAT0,etacT0
  real(16) g3A,gccA
  real(16) fGL0,fGL1,fGhL0,fGhL1,fQL0,fQL1
  real(16) q_ti,p,p_ti,costhe0,costhe1

  real(16) r_4d_to_3d,k_4d
  real(16) dmass2S,ptdmass2S

  real(16) mfs2,g_qbAq_s
  real(16) b1f1As,b2f1aAs,b1f2As,b2f2aAs,b1f3As,b3f1aAs,b3f2aAs,b2f3aAs
  real(16) fQL0_s,fQL1_s
  real(16) etaA_QL_s
  real(16) dtg_qbAq_s_gluon,dtg_qbAq_s
  real(16) dtlambdaSigmaPionB
  real(16) etaA_QL_s_T0
  real(16) f2a_s,f3a_s,f4a_s
  real(16) amfs,bmfs,cmfs,dmfs,zmfs,mfs_MeV,mfs2_T0
  real(16) fQLphi0,fQLphi1,cthinte_QL_phi

  common /strucFun/ N_str
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
  common /nbPS_com/ nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,&
                    nbSigma,nbd1xSigma,nbd2xSigma,nbd3xSigma,nbd4xSigma,nbd5xSigma
  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  common /nbGluon_com/ nbGluon,nbd1xGluon,nbd2xGluon,nbd3xGluon,nbd4xGluon,nbd5xGluon
  common /FF_thr_com/ F1F1phi,F2F1phi

  common /BF_mes_thr_com/ b1f2P,b1f2S,b1f3P,b1f3S,b2f1aP,b2f1aS,b2f2aP,b2f2aS,b3f1aP
  common /BF_thr_com/ b2f1a,b1f2,b2f2a,b1f3,b3f1a
  common /BBF_thr_com/ B2B1F1aPS,B2B1F1aSP,B1B1F2aPS,B1B1F3aPS,B2B1F2aPS,B2B1F2aSP
  common /BF_A_thr_com/ b1f1A,b2f1aA,b1f2A,b2f2aA,b1f3A,b3f1aA,b3f2aA,b2f3aA
  common /BF_AT_thr_com/ b1f1AT,b2f1aAT,b1f2AT,b2f2aAT,b1f3AT,b3f1aAT

  common /k_num_com/k_num,k_value
  common /alphas_max_com/g_max,alphas_max,k_g_max,ZAm1_max,k_ZAm1_max,g3A_max,alphas3A_max,k_g3A_max
  common /k_num1_com/k_num1,k_value1
  common /r_4d_to_3d_com/r_4d_to_3d
  common /gausslegFRGcom/w_x,y_x,w_cth,y_cth

  common /k_same_com/k_same_common

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)

  k=k_UV*exp(x)
  lam1=y(1)
  lam2=y(2)
  lam3=y(3)
  lam4=y(4)
  lam5=y(5)
  lam6=0.Q+0
  lam7=0.Q+0
  lam0=y(Nv+1)
  h=y((Nv+1)+1)
  Zphi=y((Nv+1)+(Nh+1)+1)
  Zpsi=y((Nv+1)+(Nh+1)+2)
  ZA=y((Nv+1)+(Nh+1)+3)
  c=y((Nv+1)+(Nh+1)+Nz+1)
  kappa=y((Nv+1)+(Nh+1)+Nz+2)
  g=y((Nv+1)+(Nh+1)+Nz+Nck+1)          
  !quark-gluon interaction
  g3A=y((Nv+1)+(Nh+1)+Nz+Nck+2)        
  !Three gluon interaction
  g_qbAq_s=y((Nv+1)+(Nh+1)+Nz+Nck+3)   
  !quark-gluon interaction for the strange quark


  rho=kappa 
  !calculations are performed at expansion point kappa

  zb=1.Q+0
  zf=1.Q+0

  mu0=0.Q+0
  p0=pi*T
  p0c=pi*(T-1.Q+0/hc)*exp(-k/(pi*T))+pi*1.Q+0/hc

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

  mfs2=(sqrt(mf2)+120.Q+0/(hc*k))**2  !strange quark mass

  if(abs(ms2-mp2)<1.Q-7*(ms2+mp2)/2.)then
    ms2=mp2+1.Q-4*(ms2+mp2)/2.Q+0
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l=l_com
  lb=lb_com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  mb2=0.Q+0
  call nbdx_com(mb2,T,k)
  nbGluon=nbBd0x
  nbd1xGluon=nbBd1x
  nbd2xGluon=nbBd2x
  nbd3xGluon=nbBd3x
  nbd4xGluon=nbBd4x
  nbd5xGluon=nbBd5x

  call nfdx_com(mfs2,T,mu,l,lb,k)
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

  call BF_A_thr(mfs2,T,k)  
  !for strange quark

  b1f1As=b1f1A
  b2f1aAs=b2f1aA
  b1f2As=b1f2A
  b2f2aAs=b2f2aA
  b1f3As=b1f3A
  b3f1aAs=b3f1aA
  b3f2aAs=b3f2aA
  b2f3aAs=b2f3aA

! quark loop is calculated in what follows

  fQL0=0.Q+0
  fQL1=0.Q+0
  k_same_common=0
  do i_x=1,n_x
    x_q=y_x(i_x)
    q=k*Sqrt(x_q)

    cthinte_QL=0.Q+0
    do i_cth=1,n_cth
      costhe=y_cth(i_cth)
      call FF_thr(mfs2,T,mu,l,lb,k,q,costhe,F1F1,F2F1)
      xprime=x_q+1.Q+0-2.Q+0*Sqrt(x_q)*costhe
      if(xprime<1.Q+0)then
        rF_plus_1=1.Q+0/Sqrt(xprime)
      else
        rF_plus_1=1.Q+0
      end if

      cthinte_QL=cthinte_QL+w_cth(i_cth)*((F1F1-F2F1) &
                       +(Sqrt(x_q)*costhe**2-costhe)*rF_plus_1*(F2F1-F1F1/2.))
    end do
    fQL0=fQL0+w_x(i_x)*Sqrt(x_q)*cthinte_QL
    fQL1=fQL1+w_x(i_x)*x_q*cthinte_QL
  end do
  fQL0_s=fQL0
  fQL1_s=fQL1
! this part is for strange quark loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!for strange quarks above
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

  call BF_AT_thr(mf2,T,k)
  call BF_A_thr(mf2,T,k)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  v4=1.Q+0/32.Q+0/pi**2
  C2Nc=(Nc**2-1.Q+0)/(2.Q+0*Nc)

  f3aEtaphi=f3a

  r_4d_to_3d=1.Q+0
  k_4d=r_4d_to_3d*k

  etaAT0=etaA_unqQCD_T0_fun(k_4d)
  etacT0=etac_YM_T0_fun(k_4d)

  etaAYM_T=0.Q+0
  etaAYM_T0=0.Q+0

  gccA=g

! quark loop is calculated in what follows

  fQL0=0.Q+0
  fQL1=0.Q+0
  fQLphi0=0.Q+0
  fQLphi1=0.Q+0
  k_same_common=0
  do i_x=1,n_x
    x_q=y_x(i_x)
    q=k*Sqrt(x_q)

    cthinte_QL=0.Q+0
    cthinte_QL_phi=0.Q+0
    do i_cth=1,n_cth
      costhe=y_cth(i_cth)
      call FF_thr(mf2,T,mu,l,lb,k,q,costhe,F1F1,F2F1)
      xprime=x_q+1.Q+0-2.Q+0*Sqrt(x_q)*costhe
      if(xprime<1.Q+0)then
        rF_plus_1=1.Q+0/Sqrt(xprime)
      else
        rF_plus_1=1.Q+0
      end if

      cthinte_QL=cthinte_QL+w_cth(i_cth)*((F1F1-F2F1) &
                       +(Sqrt(x_q)*costhe**2-costhe)*rF_plus_1*(F2F1-F1F1/2.))

      cthinte_QL_phi=cthinte_QL_phi+w_cth(i_cth)*(((F1F1phi-f2a)-(F2F1phi-f3a)) &
           +((Sqrt(x_q)-costhe)*rF_plus_1*F2F1-f3a)-1.Q+0/2.Q+0*((Sqrt(x_q)-costhe)*rF_plus_1*F1F1-f2a) )

    end do
    fQL0=fQL0+w_x(i_x)*Sqrt(x_q)*cthinte_QL
    fQL1=fQL1+w_x(i_x)*x_q*cthinte_QL
    fQLphi0=fQLphi0+w_x(i_x)*Sqrt(x_q)*cthinte_QL_phi
    fQLphi1=fQLphi1+w_x(i_x)*x_q*cthinte_QL_phi
  end do

  mb2=0.Q+0
  call nbdx_com(mb2,T,k)
  nbGluon=nbBd0x
  nbd1xGluon=nbBd1x
  nbd2xGluon=nbBd2x
  nbd3xGluon=nbBd3x
  nbd4xGluon=nbBd4x
  nbd5xGluon=nbBd5x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  etaphi=(((fQLphi0 - fQLphi1)*h**2*Nc*                                          &
       (4*b2f1aS*h**2 + C2Nc*                                                    &
          (9*b1f1A - 18*b1f2A + 8*b2f1aA - 2*b2f1aA*etaAT0)*g**2*Nf +            &
         4*b2f1aP*h**2*(-1 + Nf**2))*v3**2)/(6.*Nf) -                            &
    (1 + ((b1f1A - 2*b1f2A)*C2Nc*g**2*v3)/4.)*                                   &
     (2*fQLphi0*h**2*Nc*v3 - (4*b2b2PS*lam2**2*rho*v3)/(3.*k**2)))/              &
  (1 + ((b1f1A - 2*b1f2A)*C2Nc*g**2*v3)/4. +                                     &
    ((fQLphi0 - fQLphi1)*h**4*Nc*(b2f1aS + b2f1aP*(-1 + Nf**2))*v3**2)/          &
     (6.*Nf))

  etapsi=-((b2f1aS*(-4 + etaphi)*h**2 +                                          &
       C2Nc*(-9*b1f1A + 18*b1f2A + 2*b2f1aA*(-4 + (etaAT0+(etaAYM_T-etaAYM_T0))))*g**2*Nf +    &
       b2f1aP*(-4 + etaphi)*h**2*(-1 + Nf**2))*v3)/                              &
  (3.*Nf*(4 + b1f1A*C2Nc*g**2*v3 - 2*b1f2A*C2Nc*g**2*v3))

  etaA_QL=2*((-1 + etapsi)*fQL0 - etapsi*fQL1)*g**2*Nf*v3

  etaA_QL_s=2*((-1 + etapsi)*fQL0_s - etapsi*fQL1_s)*g_qbAq_s**2*1.*v3


  amfs=5.94662682612336Q+0
  bmfs=9.27021872765222Q+0
  cmfs=340.55322960490327Q+0
  dmfs=120.60128700488963Q+0

  zmfs=bmfs*(log(k*hc)-amfs)
  if(zmfs>80.Q+0)then
    mfs_MeV=dmfs
  else
    mfs_MeV=cmfs/(exp(zmfs)+1.Q+0)+dmfs
  end if

  mfs2_T0=(mfs_MeV/(k*hc))**2


  f2a_s=1.Q+0/(4.Q+0*(1 + mfs2_T0)**1.5)
  f3a_s=3.Q+0/(16.Q+0*(1 + mfs2_T0)**2.5)
  f4a_s=15.Q+0/(96.Q+0*(1 + mfs2_T0)**3.5)
  etaA_QL_s_T0=1.Q+0/(30.Q+0*pi**2)*g_qbAq_s**2*(3.Q+0*(-3.Q+0+2.Q+0*etapsi)*f2a_s+4.Q+0*(8.Q+0-3.Q+0*etapsi)*f3a_s-8.Q+0*f4a_s)

  etaA_QL_s=etaA_QL_s+etaA_QL_s_T0

  etaA=etaAT0+etaA_QL+etaA_QL_s

  call massAScreen(k,T,dmass2S,ptdmass2S)

  etaA=etaA+(2.Q+0-etaA)*dmass2S/(ZA*k**2)-ptdmass2S/(ZA*k**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  gAAA=g3A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtg_gluon=3./(4.*Nc)*v3*g**3*mf2*(b2f2aA*((2.-etaA)/3.+etaA/5.)               &
       +2.*b1f3A*(2./3.*(1.-etapsi)+etapsi/2.))                                 &
       +3./4.*Nc*v3*g**2*gAAA*(b2f1aA*((1.-etapsi)/4.+etapsi/5.)                &
       -b1f2A*(2./3.*(1.-etapsi)+etapsi/2.)-b2f2aA*(-(1.-etapsi)/6.-etapsi/10.) &
       -2.*b2f1aA*((2.-etaA)/3.+etaA/5.)-2.*b3f1aA*(-(2.-etaA)/12.-etaA/30.) )

  dtg_meson= -1./(4.*Nf)*g*h**2*v3*( (b1f2S+2.*mf2*b1f3S)*(2./3.*(1.-etapsi)    &
             +etapsi/2.)+(b2f1aS+mf2*b2f2aS)*((2.-etaphi)/3.+etaphi/5.) )       &
             -(Nf**2-1.)/(4.*Nf)*g*h**2*v3*( (b1f2P+2.*mf2*b1f3P)               &
             *(2./3.*(1.-etapsi)+etapsi/2.)                                     &
             +(b2f1aP+mf2*b2f2aP)*((2.-etaphi)/3.+etaphi/5.) )

  dtg=((etaA + 2*etapsi)*g)/2.+dtg_gluon+dtg_meson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtg_gluon_T=3./(4.*Nc)*v3*g**3*mf2*(b2f2aAT*((2.-etaA)/3.+etaA/5.)            &
       +2.*b1f3AT*(2./3.*(1.-etapsi)+etapsi/2.))                                &
       +3./4.*Nc*v3*g**2*gAAA*(b2f1aAT*((1.-etapsi)/4.+etapsi/5.)               &
       -b1f2AT*(2./3.*(1.-etapsi)+etapsi/2.)-b2f2aAT*(-(1.-etapsi)/6.-etapsi/10.) &
       -2.*b2f1aAT*((2.-etaA)/3.+etaA/5.)-2.*b3f1aAT*(-(2.-etaA)/12.-etaA/30.) )

  dtg3A=((3.*etaA)*g3A)/2.-1./(6.*pi**2)*g**3*(1.-etapsi/4.)*(1.+2.*mf2)/(1.+mf2)**4    &
        +3./(64.*pi**2)*g3A**3*(11.-2.*etaA)                                   &
        +1./(64.*pi**2)*gccA**3*(1.-etacT0/8.)                                 &
        -(1./2.)*1./(6.*pi**2)*g_qbAq_s**3*(1.-etapsi/4.)*(1.+2.*mfs2)/(1.+mfs2)**4!Note zero temperature!!!

  dtg3A=dtg3A+dtg_gluon_T

  call IRenha(k,ashift,dashift)

  dtg=dashift*g + ashift*dtg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtg_qbAq_s_gluon=3./(4.*Nc)*v3*g_qbAq_s**3*mfs2*(b2f2aAs*((2.-etaA)/3.+etaA/5.)       &
       +2.*b1f3As*(2./3.*(1.-etapsi)+etapsi/2.))                                 &
       +3./4.*Nc*v3*g_qbAq_s**2*gAAA*(b2f1aAs*((1.-etapsi)/4.+etapsi/5.)                &
       -b1f2As*(2./3.*(1.-etapsi)+etapsi/2.)-b2f2aAs*(-(1.-etapsi)/6.-etapsi/10.) &
       -2.*b2f1aAs*((2.-etaA)/3.+etaA/5.)-2.*b3f1aAs*(-(2.-etaA)/12.-etaA/30.) )

  dtg_qbAq_s=((etaA + 2*etapsi)*g_qbAq_s)/2.+dtg_qbAq_s_gluon

  dtg_qbAq_s=dashift*g_qbAq_s + ashift*dtg_qbAq_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtlam4_gluon=-(3./2.)*C2Nc*(3./4.-1./Nc**2)*(g**4)*(2./k**2)*v3*              &
       (  (b3f1aA-mf2*b3f2aA)*((2.-etaA)/3.+etaA/5.)+(b2f2aA-2.*mf2*b2f3aA)*    &
           ((1.-etapsi)/3.+etapsi/4.)  )

  call BBF_thr(mp2,ms2,mf2,T,k)

  dtlambdaSigmaPion=h**4*v3/k**2*(Nf**2-2.)/(16.*Nf*Nc)*(                        &
    (1./3.*(2.-etaphi)+1./5.*etaphi)*(B2B1F1aSP+B2B1F1aPS                        &
    -mf2*(B2B1F2aSP+B2B1F2aPS))+(2./3.*(1.-etapsi)+1./2.*etapsi)*(B1B1F2aPS      &
    -2.*mf2*B1B1F3aPS)-(2./3.*(2.-etaphi)+2./5.*etaphi)*(b3f1aP-mf2*B3F2aP)      &
    -(2./3.*(1.-etapsi)+1./2.*etapsi)*(b2f2aP-2.*mf2*B2F3aP))

  dtlambdaSigmaPionB=h**4*v3/k**2*(Nf**2-2.)/(16.*Nf*Nc)*(                        &
    (1./3.*(2.-etaphi)+1./5.*etaphi)*(B2B1F1aSP+B2B1F1aPS                        &
    -0.*(B2B1F2aSP+B2B1F2aPS))+(2./3.*(1.-etapsi)+1./2.*etapsi)*(B1B1F2aPS      &
    -2.*0.*B1B1F3aPS)-(2./3.*(2.-etaphi)+2./5.*etaphi)*(b3f1aP-0.*B3F2aP)      &
    -(2./3.*(1.-etapsi)+1./2.*etapsi)*(b2f2aP-2.*0.*B2F3aP))

  dtlam4=dtlambdaSigmaPion+dtlam4_gluon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dtA=-k**2*dtlam4/h
  
  L11P=(2*b2f1aP*(1 - etaphi/5.Q+0))/3. + (2*b1f2P*(1 - etapsi/4.))/3.
  L11S=(2*b2f1aS*(1 - etaphi/5.Q+0))/3. + (2*b1f2S*(1 - etapsi/4.))/3.
  L11A=(2*b2f1aA*(1 - etaA/5.Q+0))/3. + (2*b1f2A*(1 - etapsi/4.))/3.

  dth=(etaphi/2. + etapsi)*h - (h*(3*g**2*L11A*(-1 + Nc**2)*Nf +                &
       h**2*Nc*(-L11S + L11P*(-1 + Nf**2)))*v3)/(2.*Nc*Nf)

  dth_yukawa=dth
  dth_4fermi=-mp2*dtA

  dth=dth-mp2*dtA
!dynamical hadronization
  dhdt=dth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dZphidt=-etaphi*Zphi
  dZpsidt=-etapsi*Zpsi
  dZAdt=-etaA*ZA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dr0dtV=(k**4*v3*((2*(1 - etaphi/5.Q+0)*(0.5Q+0 + nbSigma))/(3.Q+0*Sqrt(1 + ms2)*Sqrt(zb)) +   &
      (2*(1 - etaphi/5.Q+0)*(0.5 + nbPion)*(-1 + Nf**2))/                          &
       (3.Q+0*Sqrt(1 + mp2)*Sqrt(zb)) -                                            &
      (4*(1 - etapsi/4.Q+0)*Nc*Nf*(1 - nfa - nff))/(3.Q+0*Sqrt(1 + mf2)*zf)))/2.Q+0

  dr1dtV=(k**4*v3*(((1 - etaphi/5.)*k*ms2d1rho*nbd1xSigma)/(3.*(1 + ms2)*zb) +    &
      ((1 - etaphi/5.)*k*mp2d1rho*nbd1xPion*(-1 + Nf**2))/                      &
       (3.*(1 + mp2)*zb) - ((1 - etaphi/5.)*ms2d1rho*(0.5 + nbSigma))/          &
       (3.*(1 + ms2)**1.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*mp2d1rho*(0.5 + nbPion)*(-1 + Nf**2))/                   &
       (3.*(1 + mp2)**1.5*Sqrt(zb)) +                                           &
      (2*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*(1 - nfa - nff))/                       &
       (3.*(1 + mf2)**1.5*zf) -                                                 &
      (4*(1 - etapsi/4.)*Nc*Nf*                                                 &
         (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                          &
           (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/(3.*Sqrt(1 + mf2)*zf)    &
))/2.

  dr2dtV=(k**4*v3*(((1 - etaphi/5.)*k**2*ms2d1rho**2*nbd2xSigma)/               &
       (6.*(1 + ms2)**1.5*zb**1.5) +                                            &
      ((1 - etaphi/5.)*k**2*mp2d1rho**2*nbd2xPion*(-1 + Nf**2))/                &
       (6.*(1 + mp2)**1.5*zb**1.5) -                                            &
      ((1 - etaphi/5.)*k*ms2d1rho**2*nbd1xSigma)/(2.*(1 + ms2)**2*zb) +         &
      ((1 - etaphi/5.)*k*ms2d2rho*nbd1xSigma)/(3.*(1 + ms2)*zb) -               &
      ((1 - etaphi/5.)*k*mp2d1rho**2*nbd1xPion*(-1 + Nf**2))/                   &
       (2.*(1 + mp2)**2*zb) + ((1 - etaphi/5.)*k*mp2d2rho*nbd1xPion*            &
         (-1 + Nf**2))/(3.*(1 + mp2)*zb) +                                      &
      ((1 - etaphi/5.)*ms2d1rho**2*(0.5 + nbSigma))/                            &
       (2.*(1 + ms2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*ms2d2rho*(0.5 + nbSigma))/                               &
       (3.*(1 + ms2)**1.5*Sqrt(zb)) +                                           &
      ((1 - etaphi/5.)*mp2d1rho**2*(0.5 + nbPion)*(-1 + Nf**2))/                &
       (2.*(1 + mp2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*mp2d2rho*(0.5 + nbPion)*(-1 + Nf**2))/                   &
       (3.*(1 + mp2)**1.5*Sqrt(zb)) -                                           &
      ((1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*(1 - nfa - nff))/                      &
       ((1 + mf2)**2.5*zf) + (4*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                 &
         (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                          &
           (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/                         &
       (3.*(1 + mf2)**1.5*zf) -                                                 &
      (4*(1 - etapsi/4.)*Nc*Nf*                                                 &
         (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                     &
           (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                     &
           (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                      &
           (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                     &
       (3.*Sqrt(1 + mf2)*zf)))/2.

  dr3dtV=(k**4*v3*(((1 - etaphi/5.)*k**3*ms2d1rho**3*nbd3xSigma)/               &
       (12.*(1 + ms2)**2*zb**2) +                                               &
      ((1 - etaphi/5.)*k**3*mp2d1rho**3*nbd3xPion*(-1 + Nf**2))/                &
       (12.*(1 + mp2)**2*zb**2) -                                               &
      ((1 - etaphi/5.)*k**2*ms2d1rho**3*nbd2xSigma)/                            &
       (2.*(1 + ms2)**2.5*zb**1.5) +                                            &
      ((1 - etaphi/5.)*k**2*ms2d1rho*ms2d2rho*nbd2xSigma)/                      &
       (2.*(1 + ms2)**1.5*zb**1.5) -                                            &
      ((1 - etaphi/5.)*k**2*mp2d1rho**3*nbd2xPion*(-1 + Nf**2))/                &
       (2.*(1 + mp2)**2.5*zb**1.5) +                                            &
      ((1 - etaphi/5.)*k**2*mp2d1rho*mp2d2rho*nbd2xPion*(-1 + Nf**2))/          &
       (2.*(1 + mp2)**1.5*zb**1.5) +                                            &
      (5*(1 - etaphi/5.)*k*ms2d1rho**3*nbd1xSigma)/(4.*(1 + ms2)**3*zb) -       &
      (3*(1 - etaphi/5.)*k*ms2d1rho*ms2d2rho*nbd1xSigma)/                       &
       (2.*(1 + ms2)**2*zb) + ((1 - etaphi/5.)*k*ms2d3rho*nbd1xSigma)/          &
       (3.*(1 + ms2)*zb) + (5*(1 - etaphi/5.)*k*mp2d1rho**3*nbd1xPion*          &
         (-1 + Nf**2))/(4.*(1 + mp2)**3*zb) -                                   &
      (3*(1 - etaphi/5.)*k*mp2d1rho*mp2d2rho*nbd1xPion*(-1 + Nf**2))/           &
       (2.*(1 + mp2)**2*zb) + ((1 - etaphi/5.)*k*mp2d3rho*nbd1xPion*            &
         (-1 + Nf**2))/(3.*(1 + mp2)*zb) -                                      &
      (5*(1 - etaphi/5.)*ms2d1rho**3*(0.5 + nbSigma))/                          &
       (4.*(1 + ms2)**3.5*Sqrt(zb)) +                                           &
      (3*(1 - etaphi/5.)*ms2d1rho*ms2d2rho*(0.5 + nbSigma))/                    &
       (2.*(1 + ms2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*ms2d3rho*(0.5 + nbSigma))/                               &
       (3.*(1 + ms2)**1.5*Sqrt(zb)) -                                           &
      (5*(1 - etaphi/5.)*mp2d1rho**3*(0.5 + nbPion)*(-1 + Nf**2))/              &
       (4.*(1 + mp2)**3.5*Sqrt(zb)) +                                           &
      (3*(1 - etaphi/5.)*mp2d1rho*mp2d2rho*(0.5 + nbPion)*(-1 + Nf**2))/        &
       (2.*(1 + mp2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*mp2d3rho*(0.5 + nbPion)*(-1 + Nf**2))/                   &
       (3.*(1 + mp2)**1.5*Sqrt(zb)) +                                           &
      (5*(1 - etapsi/4.)*mf2d1rho**3*Nc*Nf*(1 - nfa - nff))/                    &
       (2.*(1 + mf2)**3.5*zf) -                                                 &
      (3*(1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*                                     &
         (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                          &
           (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/((1 + mf2)**2.5*zf)      &
+ (2*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                                            &
         (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                     &
           (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                     &
           (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                      &
           (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                     &
       ((1 + mf2)**1.5*zf) - (4*(1 - etapsi/4.)*Nc*Nf*                          &
         (-(k**3*mf2d1rho**3*nfd3xa)/(8.*(1 + mf2)**1.5*zf**3) -                &
           (k**3*mf2d1rho**3*nfd3xf)/(8.*(1 + mf2)**1.5*zf**3) +                &
           (3*k**2*mf2d1rho**3*nfd2xa)/(8.*(1 + mf2)**2*zf**2) +                &
           (3*k**2*mf2d1rho**3*nfd2xf)/(8.*(1 + mf2)**2*zf**2) -                &
           (3*k*mf2d1rho**3*nfd1xa)/(8.*(1 + mf2)**2.5*zf) -                    &
           (3*k*mf2d1rho**3*nfd1xf)/(8.*(1 + mf2)**2.5*zf)))/                   &
       (3.*Sqrt(1 + mf2)*zf)))/2.

  dr4dtV=(k**4*v3*(((1 - etaphi/5.)*k**4*ms2d1rho**4*nbd4xSigma)/               &
       (24.*(1 + ms2)**2.5*zb**2.5) +                                           &
      ((1 - etaphi/5.)*k**4*mp2d1rho**4*nbd4xPion*(-1 + Nf**2))/                &
       (24.*(1 + mp2)**2.5*zb**2.5) -                                           &
      (5*(1 - etaphi/5.)*k**3*ms2d1rho**4*nbd3xSigma)/                          &
       (12.*(1 + ms2)**3*zb**2) +                                               &
      ((1 - etaphi/5.)*k**3*ms2d1rho**2*ms2d2rho*nbd3xSigma)/                   &
       (2.*(1 + ms2)**2*zb**2) -                                                &
      (5*(1 - etaphi/5.)*k**3*mp2d1rho**4*nbd3xPion*(-1 + Nf**2))/              &
       (12.*(1 + mp2)**3*zb**2) +                                               &
      ((1 - etaphi/5.)*k**3*mp2d1rho**2*mp2d2rho*nbd3xPion*(-1 + Nf**2))/       &
       (2.*(1 + mp2)**2*zb**2) +                                                &
      (15*(1 - etaphi/5.)*k**2*ms2d1rho**4*nbd2xSigma)/                         &
       (8.*(1 + ms2)**3.5*zb**1.5) -                                            &
      (3*(1 - etaphi/5.)*k**2*ms2d1rho**2*ms2d2rho*nbd2xSigma)/                 &
       ((1 + ms2)**2.5*zb**1.5) +                                               &
      ((1 - etaphi/5.)*k**2*ms2d2rho**2*nbd2xSigma)/                            &
       (2.*(1 + ms2)**1.5*zb**1.5) +                                            &
      (2*(1 - etaphi/5.)*k**2*ms2d1rho*ms2d3rho*nbd2xSigma)/                    &
       (3.*(1 + ms2)**1.5*zb**1.5) +                                            &
      (15*(1 - etaphi/5.)*k**2*mp2d1rho**4*nbd2xPion*(-1 + Nf**2))/             &
       (8.*(1 + mp2)**3.5*zb**1.5) -                                            &
      (3*(1 - etaphi/5.)*k**2*mp2d1rho**2*mp2d2rho*nbd2xPion*(-1 + Nf**2))/     &
       ((1 + mp2)**2.5*zb**1.5) +                                               &
      ((1 - etaphi/5.)*k**2*mp2d2rho**2*nbd2xPion*(-1 + Nf**2))/                &
       (2.*(1 + mp2)**1.5*zb**1.5) +                                            &
      (2*(1 - etaphi/5.)*k**2*mp2d1rho*mp2d3rho*nbd2xPion*(-1 + Nf**2))/        &
       (3.*(1 + mp2)**1.5*zb**1.5) -                                            &
      (35*(1 - etaphi/5.)*k*ms2d1rho**4*nbd1xSigma)/(8.*(1 + ms2)**4*zb) +      &
      (15*(1 - etaphi/5.)*k*ms2d1rho**2*ms2d2rho*nbd1xSigma)/                   &
       (2.*(1 + ms2)**3*zb) - (3*(1 - etaphi/5.)*k*ms2d2rho**2*nbd1xSigma)/     &
       (2.*(1 + ms2)**2*zb) - (2*(1 - etaphi/5.)*k*ms2d1rho*ms2d3rho*           &
         nbd1xSigma)/((1 + ms2)**2*zb) +                                        &
      ((1 - etaphi/5.)*k*ms2d4rho*nbd1xSigma)/(3.*(1 + ms2)*zb) -               &
      (35*(1 - etaphi/5.)*k*mp2d1rho**4*nbd1xPion*(-1 + Nf**2))/                &
       (8.*(1 + mp2)**4*zb) + (15*(1 - etaphi/5.)*k*mp2d1rho**2*mp2d2rho*       &
         nbd1xPion*(-1 + Nf**2))/(2.*(1 + mp2)**3*zb) -                         &
      (3*(1 - etaphi/5.)*k*mp2d2rho**2*nbd1xPion*(-1 + Nf**2))/                 &
       (2.*(1 + mp2)**2*zb) - (2*(1 - etaphi/5.)*k*mp2d1rho*mp2d3rho*           &
         nbd1xPion*(-1 + Nf**2))/((1 + mp2)**2*zb) +                            &
      ((1 - etaphi/5.)*k*mp2d4rho*nbd1xPion*(-1 + Nf**2))/                      &
       (3.*(1 + mp2)*zb) + (35*(1 - etaphi/5.)*ms2d1rho**4*                     &
         (0.5 + nbSigma))/(8.*(1 + ms2)**4.5*Sqrt(zb)) -                        &
      (15*(1 - etaphi/5.)*ms2d1rho**2*ms2d2rho*(0.5 + nbSigma))/                &
       (2.*(1 + ms2)**3.5*Sqrt(zb)) +                                           &
      (3*(1 - etaphi/5.)*ms2d2rho**2*(0.5 + nbSigma))/                          &
       (2.*(1 + ms2)**2.5*Sqrt(zb)) +                                           &
      (2*(1 - etaphi/5.)*ms2d1rho*ms2d3rho*(0.5 + nbSigma))/                    &
       ((1 + ms2)**2.5*Sqrt(zb)) -                                              &
      ((1 - etaphi/5.)*ms2d4rho*(0.5 + nbSigma))/                               &
       (3.*(1 + ms2)**1.5*Sqrt(zb)) +                                           &
      (35*(1 - etaphi/5.)*mp2d1rho**4*(0.5 + nbPion)*(-1 + Nf**2))/             &
       (8.*(1 + mp2)**4.5*Sqrt(zb)) -                                           &
      (15*(1 - etaphi/5.)*mp2d1rho**2*mp2d2rho*(0.5 + nbPion)*                  &
         (-1 + Nf**2))/(2.*(1 + mp2)**3.5*Sqrt(zb)) +                           &
      (3*(1 - etaphi/5.)*mp2d2rho**2*(0.5 + nbPion)*(-1 + Nf**2))/              &
       (2.*(1 + mp2)**2.5*Sqrt(zb)) +                                           &
      (2*(1 - etaphi/5.)*mp2d1rho*mp2d3rho*(0.5 + nbPion)*(-1 + Nf**2))/        &
       ((1 + mp2)**2.5*Sqrt(zb)) -                                              &
      ((1 - etaphi/5.)*mp2d4rho*(0.5 + nbPion)*(-1 + Nf**2))/                   &
       (3.*(1 + mp2)**1.5*Sqrt(zb)) -                                           &
      (35*(1 - etapsi/4.)*mf2d1rho**4*Nc*Nf*(1 - nfa - nff))/                   &
       (4.*(1 + mf2)**4.5*zf) +                                                 &
      (10*(1 - etapsi/4.)*mf2d1rho**3*Nc*Nf*                                    &
         (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                          &
           (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/((1 + mf2)**3.5*zf)      &
- (6*(1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*                                         &
         (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                     &
           (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                     &
           (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                      &
           (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                     &
       ((1 + mf2)**2.5*zf) + (8*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                 &
         (-(k**3*mf2d1rho**3*nfd3xa)/(8.*(1 + mf2)**1.5*zf**3) -                &
           (k**3*mf2d1rho**3*nfd3xf)/(8.*(1 + mf2)**1.5*zf**3) +                &
           (3*k**2*mf2d1rho**3*nfd2xa)/(8.*(1 + mf2)**2*zf**2) +                &
           (3*k**2*mf2d1rho**3*nfd2xf)/(8.*(1 + mf2)**2*zf**2) -                &
           (3*k*mf2d1rho**3*nfd1xa)/(8.*(1 + mf2)**2.5*zf) -                    &
           (3*k*mf2d1rho**3*nfd1xf)/(8.*(1 + mf2)**2.5*zf)))/                   &
       (3.*(1 + mf2)**1.5*zf) -                                                 &
      (4*(1 - etapsi/4.)*Nc*Nf*                                                 &
         (-(k**4*mf2d1rho**4*nfd4xa)/(16.*(1 + mf2)**2*zf**4) -                 &
           (k**4*mf2d1rho**4*nfd4xf)/(16.*(1 + mf2)**2*zf**4) +                 &
           (3*k**3*mf2d1rho**4*nfd3xa)/(8.*(1 + mf2)**2.5*zf**3) +              &
           (3*k**3*mf2d1rho**4*nfd3xf)/(8.*(1 + mf2)**2.5*zf**3) -              &
           (15*k**2*mf2d1rho**4*nfd2xa)/(16.*(1 + mf2)**3*zf**2) -              &
           (15*k**2*mf2d1rho**4*nfd2xf)/(16.*(1 + mf2)**3*zf**2) +              &
           (15*k*mf2d1rho**4*nfd1xa)/(16.*(1 + mf2)**3.5*zf) +                  &
           (15*k*mf2d1rho**4*nfd1xf)/(16.*(1 + mf2)**3.5*zf)))/                 &
       (3.Q+0*Sqrt(1 + mf2)*zf)))/2.Q+0

  dr5dtV=(k**4*v3*(((1 - etaphi/5.)*k**5*ms2d1rho**5*nbd5xSigma)/               &
       (48.Q+0*(1 + ms2)**3*zb**3) +                                               &
      ((1 - etaphi/5.)*k**5*mp2d1rho**5*nbd5xPion*(-1 + Nf**2))/                &
       (48.*(1 + mp2)**3*zb**3) -                                               &
      (5*(1 - etaphi/5.)*k**4*ms2d1rho**5*nbd4xSigma)/                          &
       (16.*(1 + ms2)**3.5*zb**2.5) +                                           &
      (5*(1 - etaphi/5.)*k**4*ms2d1rho**3*ms2d2rho*nbd4xSigma)/                 &
       (12.*(1 + ms2)**2.5*zb**2.5) -                                           &
      (5*(1 - etaphi/5.)*k**4*mp2d1rho**5*nbd4xPion*(-1 + Nf**2))/              &
       (16.*(1 + mp2)**3.5*zb**2.5) +                                           &
      (5*(1 - etaphi/5.)*k**4*mp2d1rho**3*mp2d2rho*nbd4xPion*(-1 + Nf**2))/     &
       (12.*(1 + mp2)**2.5*zb**2.5) +                                           &
      (35*(1 - etaphi/5.)*k**3*ms2d1rho**5*nbd3xSigma)/                         &
       (16.*(1 + ms2)**4*zb**2) -                                               &
      (25*(1 - etaphi/5.)*k**3*ms2d1rho**3*ms2d2rho*nbd3xSigma)/                &
       (6.*(1 + ms2)**3*zb**2) +                                                &
      (5*(1 - etaphi/5.)*k**3*ms2d1rho*ms2d2rho**2*nbd3xSigma)/                 &
       (4.*(1 + ms2)**2*zb**2) +                                                &
      (5*(1 - etaphi/5.)*k**3*ms2d1rho**2*ms2d3rho*nbd3xSigma)/                 &
       (6.*(1 + ms2)**2*zb**2) +                                                &
      (35*(1 - etaphi/5.)*k**3*mp2d1rho**5*nbd3xPion*(-1 + Nf**2))/             &
       (16.*(1 + mp2)**4*zb**2) -                                               &
      (25*(1 - etaphi/5.)*k**3*mp2d1rho**3*mp2d2rho*nbd3xPion*                  &
         (-1 + Nf**2))/(6.*(1 + mp2)**3*zb**2) +                                &
      (5*(1 - etaphi/5.)*k**3*mp2d1rho*mp2d2rho**2*nbd3xPion*(-1 + Nf**2))/     &
       (4.*(1 + mp2)**2*zb**2) +                                                &
      (5*(1 - etaphi/5.)*k**3*mp2d1rho**2*mp2d3rho*nbd3xPion*(-1 + Nf**2))/     &
       (6.*(1 + mp2)**2*zb**2) -                                                &
      (35*(1 - etaphi/5.)*k**2*ms2d1rho**5*nbd2xSigma)/                         &
       (4.*(1 + ms2)**4.5*zb**1.5) +                                            &
      (75*(1 - etaphi/5.)*k**2*ms2d1rho**3*ms2d2rho*nbd2xSigma)/                &
       (4.*(1 + ms2)**3.5*zb**1.5) -                                            &
      (15*(1 - etaphi/5.)*k**2*ms2d1rho*ms2d2rho**2*nbd2xSigma)/                &
       (2.*(1 + ms2)**2.5*zb**1.5) -                                            &
      (5*(1 - etaphi/5.)*k**2*ms2d1rho**2*ms2d3rho*nbd2xSigma)/                 &
       ((1 + ms2)**2.5*zb**1.5) +                                               &
      (5*(1 - etaphi/5.)*k**2*ms2d2rho*ms2d3rho*nbd2xSigma)/                    &
       (3.*(1 + ms2)**1.5*zb**1.5) +                                            &
      (5*(1 - etaphi/5.)*k**2*ms2d1rho*ms2d4rho*nbd2xSigma)/                    &
       (6.*(1 + ms2)**1.5*zb**1.5) -                                            &
      (35*(1 - etaphi/5.)*k**2*mp2d1rho**5*nbd2xPion*(-1 + Nf**2))/             &
       (4.*(1 + mp2)**4.5*zb**1.5) +                                            &
      (75*(1 - etaphi/5.)*k**2*mp2d1rho**3*mp2d2rho*nbd2xPion*                  &
         (-1 + Nf**2))/(4.*(1 + mp2)**3.5*zb**1.5) -                            &
      (15*(1 - etaphi/5.)*k**2*mp2d1rho*mp2d2rho**2*nbd2xPion*                  &
         (-1 + Nf**2))/(2.*(1 + mp2)**2.5*zb**1.5) -                            &
      (5*(1 - etaphi/5.)*k**2*mp2d1rho**2*mp2d3rho*nbd2xPion*(-1 + Nf**2))/     &
       ((1 + mp2)**2.5*zb**1.5) +                                               &
      (5*(1 - etaphi/5.)*k**2*mp2d2rho*mp2d3rho*nbd2xPion*(-1 + Nf**2))/        &
       (3.*(1 + mp2)**1.5*zb**1.5) +                                            &
      (5*(1 - etaphi/5.)*k**2*mp2d1rho*mp2d4rho*nbd2xPion*(-1 + Nf**2))/        &
       (6.*(1 + mp2)**1.5*zb**1.5) +                                            &
      (315*(1 - etaphi/5.)*k*ms2d1rho**5*nbd1xSigma)/                           &
       (16.*(1 + ms2)**5*zb) -                                                  &
      (175*(1 - etaphi/5.)*k*ms2d1rho**3*ms2d2rho*nbd1xSigma)/                  &
       (4.*(1 + ms2)**4*zb) + (75*(1 - etaphi/5.)*k*ms2d1rho*ms2d2rho**2*       &
         nbd1xSigma)/(4.*(1 + ms2)**3*zb) +                                     &
      (25*(1 - etaphi/5.)*k*ms2d1rho**2*ms2d3rho*nbd1xSigma)/                   &
       (2.*(1 + ms2)**3*zb) - (5*(1 - etaphi/5.)*k*ms2d2rho*ms2d3rho*           &
         nbd1xSigma)/((1 + ms2)**2*zb) -                                        &
      (5*(1 - etaphi/5.)*k*ms2d1rho*ms2d4rho*nbd1xSigma)/                       &
       (2.*(1 + ms2)**2*zb) + ((1 - etaphi/5.)*k*ms2d5rho*nbd1xSigma)/          &
       (3.*(1 + ms2)*zb) + (315*(1 - etaphi/5.)*k*mp2d1rho**5*nbd1xPion*        &
         (-1 + Nf**2))/(16.*(1 + mp2)**5*zb) -                                  &
      (175*(1 - etaphi/5.)*k*mp2d1rho**3*mp2d2rho*nbd1xPion*(-1 + Nf**2))/      &
       (4.*(1 + mp2)**4*zb) + (75*(1 - etaphi/5.)*k*mp2d1rho*mp2d2rho**2*       &
         nbd1xPion*(-1 + Nf**2))/(4.*(1 + mp2)**3*zb) +                         &
      (25*(1 - etaphi/5.)*k*mp2d1rho**2*mp2d3rho*nbd1xPion*(-1 + Nf**2))/       &
       (2.*(1 + mp2)**3*zb) - (5*(1 - etaphi/5.)*k*mp2d2rho*mp2d3rho*           &
         nbd1xPion*(-1 + Nf**2))/((1 + mp2)**2*zb) -                            &
      (5*(1 - etaphi/5.)*k*mp2d1rho*mp2d4rho*nbd1xPion*(-1 + Nf**2))/           &
       (2.*(1 + mp2)**2*zb) + ((1 - etaphi/5.)*k*mp2d5rho*nbd1xPion*            &
         (-1 + Nf**2))/(3.*(1 + mp2)*zb) -                                      &
      (315*(1 - etaphi/5.)*ms2d1rho**5*(0.5 + nbSigma))/                        &
       (16.*(1 + ms2)**5.5*Sqrt(zb)) +                                          &
      (175*(1 - etaphi/5.)*ms2d1rho**3*ms2d2rho*(0.5 + nbSigma))/               &
       (4.*(1 + ms2)**4.5*Sqrt(zb)) -                                           &
      (75*(1 - etaphi/5.)*ms2d1rho*ms2d2rho**2*(0.5 + nbSigma))/                &
       (4.*(1 + ms2)**3.5*Sqrt(zb)) -                                           &
      (25*(1 - etaphi/5.)*ms2d1rho**2*ms2d3rho*(0.5 + nbSigma))/                &
       (2.*(1 + ms2)**3.5*Sqrt(zb)) +                                           &
      (5*(1 - etaphi/5.)*ms2d2rho*ms2d3rho*(0.5 + nbSigma))/                    &
       ((1 + ms2)**2.5*Sqrt(zb)) +                                              &
      (5*(1 - etaphi/5.)*ms2d1rho*ms2d4rho*(0.5 + nbSigma))/                    &
       (2.*(1 + ms2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*ms2d5rho*(0.5 + nbSigma))/                               &
       (3.*(1 + ms2)**1.5*Sqrt(zb)) -                                           &
      (315*(1 - etaphi/5.)*mp2d1rho**5*(0.5 + nbPion)*(-1 + Nf**2))/            &
       (16.*(1 + mp2)**5.5*Sqrt(zb)) +                                          &
      (175*(1 - etaphi/5.)*mp2d1rho**3*mp2d2rho*(0.5 + nbPion)*                 &
         (-1 + Nf**2))/(4.*(1 + mp2)**4.5*Sqrt(zb)) -                           &
      (75*(1 - etaphi/5.)*mp2d1rho*mp2d2rho**2*(0.5 + nbPion)*                  &
         (-1 + Nf**2))/(4.*(1 + mp2)**3.5*Sqrt(zb)) -                           &
      (25*(1 - etaphi/5.)*mp2d1rho**2*mp2d3rho*(0.5 + nbPion)*                  &
         (-1 + Nf**2))/(2.*(1 + mp2)**3.5*Sqrt(zb)) +                           &
      (5*(1 - etaphi/5.)*mp2d2rho*mp2d3rho*(0.5 + nbPion)*(-1 + Nf**2))/        &
       ((1 + mp2)**2.5*Sqrt(zb)) +                                              &
      (5*(1 - etaphi/5.)*mp2d1rho*mp2d4rho*(0.5 + nbPion)*(-1 + Nf**2))/        &
       (2.*(1 + mp2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*mp2d5rho*(0.5 + nbPion)*(-1 + Nf**2))/                   &
       (3.*(1 + mp2)**1.5*Sqrt(zb)) +                                           &
      (315*(1 - etapsi/4.)*mf2d1rho**5*Nc*Nf*(1 - nfa - nff))/                  &
       (8.*(1 + mf2)**5.5*zf) -                                                 &
      (175*(1 - etapsi/4.)*mf2d1rho**4*Nc*Nf*                                   &
         (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                          &
           (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/                         &
       (4.*(1 + mf2)**4.5*zf) +                                                 &
      (25*(1 - etapsi/4.)*mf2d1rho**3*Nc*Nf*                                    &
         (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                     &
           (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                     &
           (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                      &
           (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                     &
       ((1 + mf2)**3.5*zf) - (10*(1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*             &
         (-(k**3*mf2d1rho**3*nfd3xa)/(8.*(1 + mf2)**1.5*zf**3) -                &
           (k**3*mf2d1rho**3*nfd3xf)/(8.*(1 + mf2)**1.5*zf**3) +                &
           (3*k**2*mf2d1rho**3*nfd2xa)/(8.*(1 + mf2)**2*zf**2) +                &
           (3*k**2*mf2d1rho**3*nfd2xf)/(8.*(1 + mf2)**2*zf**2) -                &
           (3*k*mf2d1rho**3*nfd1xa)/(8.*(1 + mf2)**2.5*zf) -                    &
           (3*k*mf2d1rho**3*nfd1xf)/(8.*(1 + mf2)**2.5*zf)))/                   &
       ((1 + mf2)**2.5*zf) + (10*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                &
         (-(k**4*mf2d1rho**4*nfd4xa)/(16.*(1 + mf2)**2*zf**4) -                 &
           (k**4*mf2d1rho**4*nfd4xf)/(16.*(1 + mf2)**2*zf**4) +                 &
           (3*k**3*mf2d1rho**4*nfd3xa)/(8.*(1 + mf2)**2.5*zf**3) +              &
           (3*k**3*mf2d1rho**4*nfd3xf)/(8.*(1 + mf2)**2.5*zf**3) -              &
           (15*k**2*mf2d1rho**4*nfd2xa)/(16.*(1 + mf2)**3*zf**2) -              &
           (15*k**2*mf2d1rho**4*nfd2xf)/(16.*(1 + mf2)**3*zf**2) +              &
           (15*k*mf2d1rho**4*nfd1xa)/(16.*(1 + mf2)**3.5*zf) +                  &
           (15*k*mf2d1rho**4*nfd1xf)/(16.*(1 + mf2)**3.5*zf)))/                 &
       (3.Q+0*(1 + mf2)**1.5*zf) -                                                 &
      (4*(1 - etapsi/4.Q+0)*Nc*Nf*                                                 &
         (-(k**5*mf2d1rho**5*nfd5xa)/(32.Q+0*(1 + mf2)**2.5Q+0*zf**5) -               &
           (k**5*mf2d1rho**5*nfd5xf)/(32.Q+0*(1 + mf2)**2.5Q+0*zf**5) +               &
           (5*k**4*mf2d1rho**5*nfd4xa)/(16.Q+0*(1 + mf2)**3Q+0*zf**4) +               &
           (5*k**4*mf2d1rho**5*nfd4xf)/(16.Q+0*(1 + mf2)**3Q+0*zf**4) -               &
           (45*k**3*mf2d1rho**5*nfd3xa)/(32.Q+0*(1 + mf2)**3.5Q+0*zf**3) -            &
           (45*k**3*mf2d1rho**5*nfd3xf)/(32.Q+0*(1 + mf2)**3.5Q+0*zf**3) +            &
           (105*k**2*mf2d1rho**5*nfd2xa)/(32.Q+0*(1 + mf2)**4*zf**2) +             &
           (105*k**2*mf2d1rho**5*nfd2xf)/(32.Q+0*(1 + mf2)**4*zf**2) -             &
           (105*k*mf2d1rho**5*nfd1xa)/(32.Q+0*(1 + mf2)**4.5Q+0*zf) -                 &
           (105*k*mf2d1rho**5*nfd1xf)/(32.Q+0*(1 + mf2)**4.5Q+0*zf)))/                &
       (3.Q+0*Sqrt(1 + mf2)*zf)))/2.Q+0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dcdt=(1.Q+0/2.Q+0)*etaphi*c

  dkappadt=(c*dcdt*lam1-(dr1dtV+etaphi*(lam1+kappa*lam2))*c**2)/(lam1**3+lam2*c**2)

  dlam0dt=dr0dtV+(dkappadt+etaphi*kappa)*lam1
  dlam1dt=(etaphi*lam1*lam2/2.Q+0+lam1**3/c**2*(dr1dtV+etaphi*(lam1+kappa*lam2)))/(lam1**3/c**2+lam2)

  dlam2dt=2.Q+0*etaphi*lam2+dr2dtV+(dkappadt+etaphi*kappa)*lam3
  dlam3dt=3.Q+0*etaphi*lam3+dr3dtV+(dkappadt+etaphi*kappa)*lam4
  dlam4dt=4.Q+0*etaphi*lam4+dr4dtV+(dkappadt+etaphi*kappa)*lam5
  dlam5dt=5.Q+0*etaphi*lam5+dr5dtV+(dkappadt+etaphi*kappa)*lam6
  dlam6dt=0.Q+0
  dlam7dt=0.Q+0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dydx=0.Q+00

  dydx(1)=dlam1dt
  dydx(2)=dlam2dt
  dydx(3)=dlam3dt
  dydx(4)=dlam4dt
  dydx(5)=dlam5dt
  dydx(6)=dlam6dt
  dydx(7)=dlam7dt
  dydx(Nv+1)=dlam0dt
  dydx((Nv+1)+1)=dhdt
  dydx((Nv+1)+(Nh+1)+1)=dZphidt
  dydx((Nv+1)+(Nh+1)+2)=dZpsidt
  dydx((Nv+1)+(Nh+1)+3)=dZAdt
  dydx((Nv+1)+(Nh+1)+Nz+1)=dcdt
  dydx((Nv+1)+(Nh+1)+Nz+2)=dkappadt
  dydx((Nv+1)+(Nh+1)+Nz+Nck+1)=dtg
  dydx((Nv+1)+(Nh+1)+Nz+Nck+2)=dtg3A
  dydx((Nv+1)+(Nh+1)+Nz+Nck+3)=dtg_qbAq_s


!  goto 210


  if(k_num==0.or.k<k_value)then

    k_num=k_num+1
    k_value=k

    if(g_max<g)then
      g_max=g
      alphas_max=g_max**2/(4.Q+0*pi)
      k_g_max=k*hc
    end if

   if(g3A_max<g3A)then
     g3A_max=g3A
     alphas3A_max=g3A_max**2/(4.Q+0*pi)
     k_g3A_max=k*hc
   end if

    if(ZAm1_max<1.Q+0/ZA)then
      ZAm1_max=1.Q+0/ZA
      k_ZAm1_max=k*hc
    end if

  end if

210 continue

end
