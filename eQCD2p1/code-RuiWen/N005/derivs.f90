subroutine derivs(x,y,dydx)
!Calculating the right hand side of differential equations

  implicit none

  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) Nc,Nf
  parameter(Nc=3.Q+0,Nf=2.Q+0)
  real(16) v3
  parameter(v3=1.Q+00/(2.Q+0*pi**2))
  integer NMAX 
  !maximal number of differential equations
  parameter(NMAX=50)
  real(16) x,y(NMAX),dydx(NMAX)
  integer N_str(5) 
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng
  real(16) k 
  ! IR cutoff in flow equations
  real(16) lam00,lam10,lam20,lam30,lam40,lam50,lam60,lam70
  real(16) lam01,lam11,lam21,lam31,lam41,lam51
  real(16) lam02,lam12,lam22,lam32
  real(16) lam03,lam13
  real(16) h,hlk,hsk
  real(16) Zpsi,ZA,Zphi_p0
  real(16) Z_pi,Z_K,Z_l,Z_s
  real(16) kappa1,kappa2
  real(16) rho,rho2
  real(16) etaphi,etapsi,etaA,etaphi_p0
  !meson and quark anomanous dimension
  real(16) etaA_unqQCD_T0_fun,etac_YM_T0_fun
  external etaA_unqQCD_T0_fun,etac_YM_T0_fun
  real(16) p0,p0c 
  !temporal compontent of external momentum
  real(16) T,mu
  real(16) mu0
  real(16) l,lb 
  !polyakov loop
  real(16) k_UV,k_IR,t_UV,t_IR
  real(16) mp2,ms2,mf2
  real(16) dr00dtv,dr10dtV,dr20dtV,dr30dtV,dr40dtV,dr50dtV
  real(16) dr01dtV,dr11dtV,dr21dtV,dr31dtV
  real(16) dr02dtV,dr12dtV
  real(16) dmb2drho(2,12),dmf2drho(2,4)
  real(16) nbdnx(2,6),nffdnx(2,6),nfadnx(2,6)
  real(16) lbt(2,12),lft(2,12)
  real(16) lbt00(2),lbt10(2),lbt20(2),lbt30(2),lbt40(2),lbt50(2)
  real(16) lbt01(2),lbt11(2),lbt21(2),lbt31(2),lbt02(2),lbt12(2)
  real(16) lft00(2),lft10(2),lft20(2),lft30(2),lft40(2),lft50(2)
  real(16) lft01(2),lft11(2),lft21(2),lft31(2),lft02(2),lft12(2)
  !terms of drdtV ,terms of function 26
  real(16) k44pi2
  !calculate accelerate,no sence
  real(16) dlam00dt,dlam10dt,dlam20dt,dlam30dt,dlam40dt,dlam50dt,dlam60dt
  real(16) dlam70dt,dlam01dt,dlam11dt,dlam21dt,dlam31dt,dlam41dt,dlam51dt
  real(16) dlam02dt,dlam12dt,dlam22dt,dlam32dt
  real(16) dlam03dt,dlam13dt
  real(16) dth,dhdt
  real(16) dZ_pidt,dZpsidt,dZAdt,dZphi_p0dt
  real(16) dkappa1dt,dkappa2dt
  real(16) l_com,lb_com
  real(16) B2B1F1aPS,B2B1F1aSP,B1B1F2aPS,B1B1F3aPS,B2B1F2aPS,B2B1F2aSP,B2F3aP,B3F2aP

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
  real(16) f2a,f3a
  real(16) F1F1,F2F1
  real(16) f2aT,f3aT
  real(16) b2f1a,b1f2,b2f2a,b1f3,b3f1a
  real(16) b2b2

  real(16) dtg,dtlam4,dtg3A

  real(16) g,gAAA
  real(16) v4,C2Nc,ashift,dashift
  real(16) N21mP,N21mS,N21gA,N12A,L12A,L111PS,L11gA
  real(16) dtlam4_gluon
  real(16) dtg_gluon,dtg_meson,dtg_gluon_T

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

  real(16) etaAT0,etacT0
  real(16) g3A,gccA
  real(16) fGL0,fGL1,fGhL0,fGhL1,fQL0,fQL1

  real(16) r_4d_to_3d,k_4d
  real(16) dmass2S,ptdmass2S

  real(16) mfs2,g_qbAq_s
  real(16) b1f1As,b2f1aAs,b1f2As,b2f2aAs,b1f3As,b3f1aAs,b3f2aAs,b2f3aAs
  real(16) fQL0_s,fQL1_s
  real(16) etaA_QL_s
  real(16) dtg_qbAq_s_gluon,dtg_qbAq_s
  real(16) etaA_QL_s_T0
  real(16) f2a_s,f3a_s,f4a_s
  real(16) amfs,bmfs,cmfs,dmfs,zmfs,mfs_MeV,mfs2_T0
  real(16) fQLphi0,fQLphi1,cthinte_QL_phi

  real(16) jl,js,Sl,Ss,ck
  real(16) djldt,djsdt,dSldt,dSsdt,dckdt
  real(16) nff_s,nfd1xf_s,nfd2xf_s,nfd3xf_s,nfd4xf_s,nfd5xf_s
  real(16) nfa_s,nfd1xa_s,nfd2xa_s,nfd3xa_s,nfd4xa_s,nfd5xa_s

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

  common /BF_mes_thr_com/ b1f2P,b1f2S,b1f3P,b1f3S,b2f1aP,b2f1aS,b2f2aP,b2f2aS,b3f1aP,b2f3aP,b3f2aP
  common /BF_thr_com/ b2f1a,b1f2,b2f2a,b1f3,b3f1a
  common /BBF_thr_com/ B2B1F1aPS,B2B1F1aSP,B1B1F2aPS,B1B1F3aPS,B2B1F2aPS,B2B1F2aSP
  common /BF_A_thr_com/ b1f1A,b2f1aA,b1f2A,b2f2aA,b1f3A,b3f1aA,b3f2aA,b2f3aA
  common /BF_AT_thr_com/ b1f1AT,b2f1aAT,b1f2AT,b2f2aAT,b1f3AT,b3f1aAT

  common /r_4d_to_3d_com/r_4d_to_3d
  common /gausslegFRGcom/w_x,y_x,w_cth,y_cth

  common /k_same_com/k_same_common

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)

  k=k_UV*exp(x)
  
  lam10=y(1)
  lam20=y(2)
  lam30=y(3)
  lam40=y(4)
  lam50=y(5)
  lam60=y(6)
  lam70=y(7)
  lam01=y(8)
  lam11=y(9)
  lam21=y(10)
  lam31=y(11)
  lam41=y(12)
  lam51=y(13)
  lam02=y(14)
  lam12=y(15)
  lam22=y(16)
  lam32=y(17)
  lam03=y(18)
  lam13=y(19)
  lam00=y(Nv+1)
  h=y((Nv+1)+1)
  Z_pi =y((Nv+1)+(Nh+2)+1)
  Z_K  =y((Nv+1)+(Nh+2)+2)
  Zpsi =y((Nv+1)+(Nh+2)+3)
  Z_s  =y((Nv+1)+(Nh+2)+4)
  ZA   =y((Nv+1)+(Nh+2)+5)
  Zphi_p0=y((Nv+1)+(Nh+2)+6)
  jl=y((Nv+1)+(Nh+2)+Nz+1)
  js=y((Nv+1)+(Nh+2)+Nz+2)
  Sl=y((Nv+1)+(Nh+2)+Nz+3)
  Ss=y((Nv+1)+(Nh+2)+Nz+4)
  ck=y((Nv+1)+(Nh+2)+Nz+5)
  g=y((Nv+1)+(Nh+2)+Nz+Nck+1)          
  !quark-gluon interaction
  g3A=y((Nv+1)+(Nh+2)+Nz+Nck+2)        
  !Three gluon interaction
  g_qbAq_s=y((Nv+1)+(Nh+2)+Nz+Nck+3)   
  !quark-gluon interaction for the strange quark

  kappa1=(Sl**2 + Ss**2)/2.Q+00
  kappa2=(Sl**2 - 2*Ss**2)**2/24.Q+00
  rho=kappa1
  rho2=kappa2
  !calculations are performed at expansion point kappa

  hlk=h
  hsk=h

  mu0=0.Q+0
  p0=pi*T
  p0c=pi*(T-1.Q+0/hc)*exp(-k/(pi*T))+pi*1.Q+0/hc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call mass2(y,dmb2drho) 

  !mass and their derivatives
  dmf2drho(1,1)=(hlk**2.Q+00*Sl**2.Q+00)/4.Q+00
  dmf2drho(1,2)=hlk**2.Q+00/3.Q+00
  dmf2drho(1,3)=hlk**2.Q+00/(Sl**2.Q+00 - 2.Q+00*Ss**2.Q+00)
  dmf2drho(1,4)=(-12.Q+00*hlk**2.Q+00)/(Sl**2.Q+00 - 2.Q+00*Ss**2.Q+00)**3.Q+00
  dmf2drho(2,1)=(hsk**2.Q+00*Ss**2.Q+00)/2.Q+00
  dmf2drho(2,2)=hsk**2.Q+00/3.Q+00
  dmf2drho(2,3)=(-2.Q+00*hsk**2.Q+00)/(Sl**2.Q+00 - 2.Q+00*Ss**2.Q+00)
  dmf2drho(2,4)=(24.Q+00*hsk**2.Q+00)/(Sl**2.Q+00 - 2.Q+00*Ss**2.Q+00)**3.Q+00
 
  dmb2drho=dmb2drho/k**2
  dmf2drho=dmf2drho/k**2

  mp2=dmb2drho(1,1)
  ms2=dmb2drho(2,1)
  mf2=dmf2drho(1,1)
  mfs2=dmf2drho(2,1)

  if(abs(ms2-mp2)<1.Q-7*(ms2+mp2)/2.Q+0)then
    ms2=mp2+1.Q-4*(ms2+mp2)/2.Q+0
  end if

!  V_pipiSg=lam2
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

  nff_s=nff
  nfd1xf_s=nfd1xf
  nfd2xf_s=nfd2xf
  nfd3xf_s=nfd3xf
  nfd4xf_s=nfd4xf
  nfd5xf_s=nfd5xf

  nfa_s=nfa
  nfd1xa_s=nfd1xa
  nfd2xa_s=nfd2xa
  nfd3xa_s=nfd3xa
  nfd4xa_s=nfd4xa
  nfd5xa_s=nfd5xa

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
  v4=1.Q+0/32.Q+0/pi**2.Q+0
  C2Nc=(Nc**2-1.Q+0)/(2.Q+0*Nc)

  r_4d_to_3d=1.Q+0
  k_4d=r_4d_to_3d*k

  etaAT0=etaA_unqQCD_T0_fun(k_4d)
  etacT0=etac_YM_T0_fun(k_4d)

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

      cthinte_QL_phi=cthinte_QL_phi+w_cth(i_cth)*(((F1F1-f2a)-(F2F1-f3a)) &
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
         4*b2f1aP*h**2*(-1 + Nf**2))*v3**2)/(6.Q+0*Nf) -                            &
    (1 + ((b1f1A - 2*b1f2A)*C2Nc*g**2*v3)/4.)*                                   &
     (2*fQLphi0*h**2*Nc*v3 - (4*b2b2PS*lam20**2*rho*v3)/(3.Q+0*k**2)))/              &
  (1 + ((b1f1A - 2*b1f2A)*C2Nc*g**2*v3)/4. +                                     &
    ((fQLphi0 - fQLphi1)*h**4*Nc*(b2f1aS + b2f1aP*(-1 + Nf**2))*v3**2)/          &
     (6.Q+0*Nf))

  etapsi=-((b2f1aS*(-4 + etaphi)*h**2 +                                          &
       C2Nc*(-9*b1f1A + 18*b1f2A + 2*b2f1aA*(-4 + (etaAT0)))*g**2*Nf +    &
       b2f1aP*(-4 + etaphi)*h**2*(-1 + Nf**2))*v3)/                              &
  (3.Q+0*Nf*(4 + b1f1A*C2Nc*g**2*v3 - 2*b1f2A*C2Nc*g**2*v3))

  etaphi_p0=((((-3 + 2*etapsi)*f2a - 4*(-2 + etapsi)*f3a)*h**2*Nc +              &
            (4*b2b2PS*lam20**2*rho)/k**2)*v3)/3.

  etaA_QL=2*((-1 + etapsi)*fQL0 - etapsi*fQL1)*g**2*Nf*v3

  etaA_QL_s=2*((-1 + etapsi)*fQL0_s - etapsi*fQL1_s)*g_qbAq_s**2*1.Q+0*v3


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

  dtg_gluon=3./(4.Q+0*Nc)*v3*g**3*mf2*(b2f2aA*((2.-etaA)/3.+etaA/5.)               &
       +2.Q+0*b1f3A*(2./3.Q+0*(1.-etapsi)+etapsi/2.))                                 &
       +3./4.Q+0*Nc*v3*g**2*gAAA*(b2f1aA*((1.-etapsi)/4.+etapsi/5.)                &
       -b1f2A*(2./3.Q+0*(1.-etapsi)+etapsi/2.)-b2f2aA*(-(1.-etapsi)/6.-etapsi/10.) &
       -2.Q+0*b2f1aA*((2.-etaA)/3.+etaA/5.)-2.Q+0*b3f1aA*(-(2.-etaA)/12.-etaA/30.) )

  dtg_meson= -1./(4.Q+0*Nf)*g*h**2*v3*( (b1f2S+2.Q+0*mf2*b1f3S)*(2./3.Q+0*(1.-etapsi)    &
             +etapsi/2.)+(b2f1aS+mf2*b2f2aS)*((2.-etaphi)/3.+etaphi/5.) )       &
             -(Nf**2-1.)/(4.Q+0*Nf)*g*h**2*v3*( (b1f2P+2.Q+0*mf2*b1f3P)               &
             *(2./3.Q+0*(1.-etapsi)+etapsi/2.)                                     &
             +(b2f1aP+mf2*b2f2aP)*((2.-etaphi)/3.+etaphi/5.) )

  dtg=((etaA + 2*etapsi)*g)/2.+dtg_gluon+dtg_meson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtg_gluon_T=3./(4.Q+0*Nc)*v3*g**3*mf2*(b2f2aAT*((2.-etaA)/3.+etaA/5.)            &
       +2.Q+0*b1f3AT*(2./3.Q+0*(1.-etapsi)+etapsi/2.))                                &
       +3./4.Q+0*Nc*v3*g**2*gAAA*(b2f1aAT*((1.-etapsi)/4.+etapsi/5.)               &
       -b1f2AT*(2./3.Q+0*(1.-etapsi)+etapsi/2.)-b2f2aAT*(-(1.-etapsi)/6.-etapsi/10.) &
       -2.Q+0*b2f1aAT*((2.-etaA)/3.+etaA/5.)-2.Q+0*b3f1aAT*(-(2.-etaA)/12.-etaA/30.) )

  dtg3A=((3.Q+0*etaA)*g3A)/2.-1./(6.Q+0*pi**2)*g**3*(1.-etapsi/4.)*(1.+2.Q+0*mf2)/(1.+mf2)**4    &
        +3./(64.Q+0*pi**2)*g3A**3*(11.-2.Q+0*etaA)                                   &
        +1./(64.Q+0*pi**2)*gccA**3*(1.-etacT0/8.)                                 &
        -(1./2.)*1./(6.Q+0*pi**2)*g_qbAq_s**3*(1.-etapsi/4.)*(1.+2.Q+0*mfs2)/(1.+mfs2)**4!Note zero temperature!!!

  dtg3A=dtg3A+dtg_gluon_T

  call IRenha(k,ashift,dashift)

  dtg=dashift*g + ashift*dtg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtg_qbAq_s_gluon=3./(4.Q+0*Nc)*v3*g_qbAq_s**3*mfs2*(b2f2aAs*((2.-etaA)/3.+etaA/5.)       &
       +2.Q+0*b1f3As*(2./3.Q+0*(1.-etapsi)+etapsi/2.))                                 &
       +3./4.Q+0*Nc*v3*g_qbAq_s**2*gAAA*(b2f1aAs*((1.-etapsi)/4.+etapsi/5.)                &
       -b1f2As*(2./3.Q+0*(1.-etapsi)+etapsi/2.)-b2f2aAs*(-(1.-etapsi)/6.-etapsi/10.) &
       -2.Q+0*b2f1aAs*((2.-etaA)/3.+etaA/5.)-2.Q+0*b3f1aAs*(-(2.-etaA)/12.-etaA/30.) )

  dtg_qbAq_s=((etaA + 2*etapsi)*g_qbAq_s)/2.+dtg_qbAq_s_gluon

  dtg_qbAq_s=dashift*g_qbAq_s + ashift*dtg_qbAq_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtlam4_gluon=-(3./2.)*C2Nc*(3./4.-1./Nc**2)*(g**4)*(2./k**2)*v3*              &
       (  (b3f1aA-mf2*b3f2aA)*((2.-etaA)/3.+etaA/5.)+(b2f2aA-2.Q+0*mf2*b2f3aA)*    &
           ((1.-etapsi)/3.+etapsi/4.)  )

  call BBF_thr(mp2,ms2,mf2,T,k)

  dtlambdaSigmaPion=h**4*v3/k**2*(Nf**2-2.)/(16.Q+0*Nf*Nc)*(                        &
    (1./3.Q+0*(2.-etaphi)+1./5.Q+0*etaphi)*(B2B1F1aSP+B2B1F1aPS                        &
    -mf2*(B2B1F2aSP+B2B1F2aPS))+(2./3.Q+0*(1.-etapsi)+1./2.Q+0*etapsi)*(B1B1F2aPS      &
    -2.Q+0*mf2*B1B1F3aPS)-(2./3.Q+0*(2.-etaphi)+2./5.Q+0*etaphi)*(b3f1aP-mf2*B3F2aP)      &
    -(2./3.Q+0*(1.-etapsi)+1./2.Q+0*etapsi)*(b2f2aP-2.Q+0*mf2*B2F3aP))

  dtlam4=dtlambdaSigmaPion+dtlam4_gluon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dtA=-k**2*dtlam4/h
  
  L11P=(2*b2f1aP*(1 - etaphi/5.Q+0))/3. + (2*b1f2P*(1 - etapsi/4.))/3.
  L11S=(2*b2f1aS*(1 - etaphi/5.Q+0))/3. + (2*b1f2S*(1 - etapsi/4.))/3.
  L11A=(2*b2f1aA*(1 - etaA/5.Q+0))/3. + (2*b1f2A*(1 - etapsi/4.))/3.

  dth=(etaphi/2. + etapsi)*h - (h*(3*g**2*L11A*(-1 + Nc**2)*Nf +                &
       h**2*Nc*(-L11S + L11P*(-1 + Nf**2)))*v3)/(2.Q+0*Nc*Nf)
  
  dth=dth-mp2*dtA
!dynamical hadronization
  dhdt=dth
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dZ_pidt=-etaphi*Z_pi
  dZpsidt=-etapsi*Zpsi
  dZAdt=-etaA*ZA
  dZphi_p0dt=-etaphi_p0*Z_pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !WR

  nbdnx(:,1)=(/nbPion,nbSigma/)
  nbdnx(:,2)=(/nbd1xPion,nbd1xSigma/)
  nbdnx(:,3)=(/nbd2xPion,nbd2xSigma/)
  nbdnx(:,4)=(/nbd3xPion,nbd3xSigma/)
  nbdnx(:,5)=(/nbd4xPion,nbd4xSigma/)
  nbdnx(:,6)=(/nbd5xPion,nbd5xSigma/)

  nffdnx(:,1)=(/nff,nff_s/)
  nffdnx(:,2)=(/nfd1xf,nfd1xf_s/)
  nffdnx(:,3)=(/nfd2xf,nfd2xf_s/)
  nffdnx(:,4)=(/nfd3xf,nfd3xf_s/)
  nffdnx(:,5)=(/nfd4xf,nfd4xf_s/)
  nffdnx(:,6)=(/nfd5xf,nfd5xf_s/)
  
  nfadnx(:,1)=(/nfa,nfa_s/)
  nfadnx(:,2)=(/nfd1xa,nfd1xa_s/)
  nfadnx(:,3)=(/nfd2xa,nfd2xa_s/)
  nfadnx(:,4)=(/nfd3xa,nfd3xa_s/)
  nfadnx(:,5)=(/nfd4xa,nfd4xa_s/)
  nfadnx(:,6)=(/nfd5xa,nfd5xa_s/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !WR
  call lbtem(k,dmb2drho,nbdnx,etaphi,lbt)
  lbt00=lbt(:,1)
  lbt10=lbt(:,2)
  lbt20=lbt(:,3)
  lbt30=lbt(:,4)
  lbt40=lbt(:,5)
  lbt50=lbt(:,6)
  lbt01=lbt(:,7)
  lbt11=lbt(:,8)
  lbt21=lbt(:,9)
  lbt31=lbt(:,10)
  lbt02=lbt(:,11)
  lbt12=lbt(:,12)

  call lftem(k,dmf2drho,nffdnx,nfadnx,etapsi,lft)
  lft00=lft(:,1)
  lft10=lft(:,2)
  lft20=lft(:,3)
  lft30=lft(:,4)
  lft40=lft(:,5)
  lft50=lft(:,6)
  lft01=lft(:,7)
  lft11=lft(:,8)
  lft21=lft(:,9)
  lft31=lft(:,10)
  lft02=lft(:,11)
  lft12=lft(:,12)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !WR
  k44pi2=k**4.Q+00/4.Q+00/pi/pi

  dr00dtV=k44pi2*(3*lbt00(1)+lbt00(2)-4.Q+00*Nc*(2*(lft00(1))+lft00(2)))
  dr10dtV=k44pi2*(3*lbt10(1)+lbt10(2)-4.Q+00*Nc*(2*(lft10(1))+lft10(2)))
  dr20dtV=k44pi2*(3*lbt20(1)+lbt20(2)-4.Q+00*Nc*(2*(lft20(1))+lft20(2)))
  dr30dtV=k44pi2*(3*lbt30(1)+lbt30(2)-4.Q+00*Nc*(2*(lft30(1))+lft30(2)))
  dr40dtV=k44pi2*(3*lbt40(1)+lbt40(2)-4.Q+00*Nc*(2*(lft40(1))+lft40(2)))
  dr50dtV=k44pi2*(3*lbt50(1)+lbt50(2)-4.Q+00*Nc*(2*(lft50(1))+lft50(2)))
  dr01dtV=k44pi2*(3*lbt01(1)+lbt01(2)-4.Q+00*Nc*(2*(lft01(1))+lft01(2)))
  dr11dtV=k44pi2*(3*lbt11(1)+lbt11(2)-4.Q+00*Nc*(2*(lft11(1))+lft11(2)))
  dr21dtV=k44pi2*(3*lbt21(1)+lbt21(2)-4.Q+00*Nc*(2*(lft21(1))+lft21(2)))
  dr31dtV=k44pi2*(3*lbt31(1)+lbt31(2)-4.Q+00*Nc*(2*(lft31(1))+lft31(2)))
  dr02dtV=k44pi2*(3*lbt02(1)+lbt02(2)-4.Q+00*Nc*(2*(lft02(1))+lft02(2)))
  dr12dtV=k44pi2*(3*lbt12(1)+lbt12(2)-4.Q+00*Nc*(2*(lft12(1))+lft12(2)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  djldt=(1.Q+0/2.Q+0)*etaphi*jl
  djsdt=(1.Q+0/2.Q+0)*etaphi*js

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dSldt=-(((dr01dtV + 2.Q+0*etaphi*lam01 + 2.Q+0*etaphi*kappa2*lam02 +                   &
           etaphi*kappa1*lam11 +                                                &
           (etaphi*js*Sl - etaphi*jl*Ss)/(Sl**3*Ss - 2.Q+0*Sl*Ss**3))*              &
         (lam20*Ss - (lam11*Ss*(Sl**2 - 2.Q+0*Ss**2.Q+0))/3.Q+0 +                          &
           (4*js + Sqrt(2.Q+0)*ck*(Sl**2 - 4*Ss**2.Q+0))/(12.Q+0*Ss**2.Q+0)) -                 &
        (dr10dtV + (etaphi*(6*lam10 + 12.Q+0*kappa2*lam11 + 6*kappa1*lam20 -        &
                (2.Q+0*jl)/Sl - js/Ss))/6.Q+0)*                                        &
         (lam11*Ss - (lam02*Ss*(Sl**2 - 2.Q+0*Ss**2.Q+0))/3.Q+0 -                          &
           (16*jl*Ss**3 + Sqrt(2.Q+0)*ck*Sl*(Sl**2 - 2.Q+0*Ss**2.Q+0)**2 +                  &
              4*js*(Sl**3 - 6*Sl*Ss**2.Q+0))/(2.Q+0*Sl*Ss**2.Q+0*(Sl**2 - 2.Q+0*Ss**2.Q+0)**2.Q+0)     &
))/((lam11*Sl + (lam02*Sl*(Sl**2 - 2.Q+0*Ss**2.Q+0))/6.Q+0 +                               &
           (-4*js*Sl**3 + 6*jl*Sl**2.Q+0*Ss - 4*jl*Ss**3)/                          &
            (Ss*(Sl**3 - 2.Q+0*Sl*Ss**2.Q+0)**2.Q+0))*                                      &
         (lam20*Ss - (lam11*Ss*(Sl**2 - 2.Q+0*Ss**2.Q+0))/3.Q+0 +                          &
           (4*js + Sqrt(2.Q+0)*ck*(Sl**2 - 4*Ss**2.Q+0))/(12.Q+0*Ss**2.Q+0)) -                 &
        ((2.Q+0*jl)/(3.Q+0*Sl**2.Q+0) + (Sl*                                               &
              (-(Sqrt(2.Q+0)*ck) + 6*lam20*Ss + lam11*Sl**2.Q+0*Ss -                    &
                2.Q+0*lam11*Ss**3))/(6.Q+0*Ss))*                                       &
         (lam11*Ss - (lam02*Ss*(Sl**2 - 2.Q+0*Ss**2.Q+0))/3.Q+0 -                          &
           (16*jl*Ss**3 + Sqrt(2.Q+0)*ck*Sl*(Sl**2 - 2.Q+0*Ss**2.Q+0)**2 +                  &
              4*js*(Sl**3 - 6*Sl*Ss**2.Q+0))/(2.Q+0*Sl*Ss**2.Q+0*(Sl**2 - 2.Q+0*Ss**2.Q+0)**2.Q+0))    &
))
   dSsdt=(2.Q+0*Sl*Ss*(6*Ss*(dr10dtV*                                                &
            (24*js*Sl**3 + Ss*                                                  &
               (jl*(-36*Sl**2 + 24*Ss**2.Q+0) -                                     &
                 Sl**3*(Sl**2 - 2.Q+0*Ss**2.Q+0)**2.Q+0*                                    &
                  (6*lam11 + lam02*Sl**2 - 2.Q+0*lam02*Ss**2.Q+0))) +                   &
           dr01dtV*(Sl**2 - 2.Q+0*Ss**2.Q+0)**2.Q+0*                                        &
            (-(Sqrt(2.Q+0)*ck*Sl**3) +                                              &
              Ss*(4*jl + Sl**3*(6*lam20 + lam11*Sl**2 - 2.Q+0*lam11*Ss**2.Q+0)))) +     &
        etaphi*(-24*js**2.Q+0*Sl**3 +                                               &
           js*(-6*Sqrt(2.Q+0)*ck*Sl**3*(Sl**2 - 2.Q+0*Ss**2.Q+0) +                          &
              Ss*(12.Q+0*jl*(Sl**2 - 6*Ss**2.Q+0) +                                     &
                 Sl**3*(144*lam10 + 288*kappa2*lam11 + 144*kappa1*lam20 +       &
                    36*lam20*Sl**2 + 12.Q+0*lam11*Sl**4 + lam02*Sl**6 -             &
                    72.Q+0*lam20*Ss**2 - 48*lam11*Sl**2.Q+0*Ss**2 -                     &
                    6*lam02*Sl**4*Ss**2 + 48*lam11*Ss**4 +                      &
                    12.Q+0*lam02*Sl**2.Q+0*Ss**4 - 8*lam02*Ss**6))) +                   &
           2.Q+0*Ss*(-3*Sqrt(2.Q+0)*ck*Sl**2.Q+0*(Sl**2 - 2.Q+0*Ss**2.Q+0)*                         &
               (-jl + (2.Q+0*lam01 + 2.Q+0*kappa2*lam02 + kappa1*lam11)*Sl*             &
                  (Sl**2 - 2.Q+0*Ss**2.Q+0)) +                                          &
              Ss*(24*jl**2.Q+0*Sl -                                                 &
                 3*Sl**3*(Sl**2 - 2.Q+0*Ss**2.Q+0)**2.Q+0*                                  &
                  (-12.Q+0*lam01*lam20 + 12.Q+0*kappa2*(lam11**2 - lam02*lam20) -       &
                    2.Q+0*lam01*lam11*Sl**2 - kappa1*lam11**2.Q+0*Sl**2 +               &
                    kappa1*lam02*lam20*Sl**2 + 4*lam01*lam11*Ss**2 +            &
                    2.Q+0*kappa1*lam11**2.Q+0*Ss**2 - 2.Q+0*kappa1*lam02*lam20*Ss**2 +      &
                    lam10*(6*lam11 + lam02*(Sl**2 - 2.Q+0*Ss**2.Q+0))) +                &
                 jl*(-108*kappa1*lam20*Sl**2 + 24*lam01*Sl**4 +                 &
                    12.Q+0*kappa1*lam11*Sl**4 - 18*lam20*Sl**4 +                    &
                    3*lam11*Sl**6 + lam02*Sl**8 + 72.Q+0*kappa1*lam20*Ss**2 -       &
                    96*lam01*Sl**2.Q+0*Ss**2 - 48*kappa1*lam11*Sl**2.Q+0*Ss**2 +        &
                    36*lam20*Sl**2.Q+0*Ss**2 - 12.Q+0*lam11*Sl**4*Ss**2 -               &
                    6*lam02*Sl**6*Ss**2 + 96*lam01*Ss**4 +                      &
                    48*kappa1*lam11*Ss**4 + 12.Q+0*lam11*Sl**2.Q+0*Ss**4 +              &
                    12.Q+0*lam02*Sl**4*Ss**4 - 8*lam02*Sl**2.Q+0*Ss**6 -                &
                    36*lam10*(3*Sl**2 - 2.Q+0*Ss**2.Q+0) -                              &
                    24*kappa2*(lam11*(9*Sl**2 - 6*Ss**2.Q+0) -                      &
                       lam02*(Sl**2 - 2.Q+0*Ss**2.Q+0)**2.Q+0)))))))/                       &
    (-96*js**2.Q+0*Sl**4 - 12.Q+0*ck**2.Q+0*Sl**4*(Sl**2 - 2.Q+0*Ss**2.Q+0)**2 +                    &
      Sqrt(2.Q+0)*ck*Sl*Ss*(12.Q+0*jl*(5*Sl**4 - 30*Sl**2.Q+0*Ss**2 + 16*Ss**4) +           &
         Sl**3*(Sl**2 - 2.Q+0*Ss**2.Q+0)**2.Q+0*                                            &
          (36*lam20 + (Sl**2 - 2.Q+0*Ss**2.Q+0)*                                        &
             (12.Q+0*lam11 + lam02*Sl**2 - 8*lam02*Ss**2.Q+0))) +                       &
      4*Ss**4*(96*jl**2 - 9*(lam11**2 - lam02*lam20)*Sl**4*                     &
          (Sl**2 - 2.Q+0*Ss**2.Q+0)**3 +                                                &
         4*jl*Sl*(9*lam20*(7*Sl**2 - 2.Q+0*Ss**2.Q+0) +                                 &
            (Sl**2 - 2.Q+0*Ss**2.Q+0)**2.Q+0*(-6*lam11 + lam02*(Sl**2 - 2.Q+0*Ss**2.Q+0)))) -       &
      4*js*Sl*(12.Q+0*Sqrt(2.Q+0)*ck*Sl**3*(Sl**2 - 5*Ss**2.Q+0) +                          &
         Ss*(jl*(-60*Sl**2 + 168*Ss**2.Q+0) -                                       &
            Sl**3*(36*lam20*(Sl**2 - 8*Ss**2.Q+0) +                                 &
 (Sl**2 - 2.Q+0*Ss**2.Q+0)**2.Q+0*(12.Q+0*lam11 + lam02*Sl**2 - 2.Q+0*lam02*Ss**2.Q+0)))))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dkappa1dt=dSldt*Sl+dSsdt*Ss
  dkappa2dt=((dSldt*Sl - 2*dSsdt*Ss)*(Sl**2 - 2*Ss**2))/6.Q+00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dlam00dt=dr00dtV+(dkappa1dt+etaphi*kappa1)*lam10+(dkappa2dt+2*etaphi*kappa2)*lam01+0*etaphi*lam00
  dlam10dt=dr10dtV+(dkappa1dt+etaphi*kappa1)*lam20+(dkappa2dt+2*etaphi*kappa2)*lam11+1*etaphi*lam10
  dlam20dt=dr20dtV+(dkappa1dt+etaphi*kappa1)*lam30+(dkappa2dt+2*etaphi*kappa2)*lam21+2*etaphi*lam20
  dlam30dt=dr30dtV+(dkappa1dt+etaphi*kappa1)*lam40+(dkappa2dt+2*etaphi*kappa2)*lam31+3*etaphi*lam30
  dlam40dt=dr40dtV+(dkappa1dt+etaphi*kappa1)*lam50+(dkappa2dt+2*etaphi*kappa2)*lam41+4*etaphi*lam40
  dlam50dt=dr50dtV+(dkappa1dt+etaphi*kappa1)*lam60+(dkappa2dt+2*etaphi*kappa2)*lam51+5*etaphi*lam50
  dlam01dt=dr01dtV+(dkappa1dt+etaphi*kappa1)*lam11+(dkappa2dt+2*etaphi*kappa2)*lam02+2*etaphi*lam01
  dlam11dt=dr11dtV+(dkappa1dt+etaphi*kappa1)*lam21+(dkappa2dt+2*etaphi*kappa2)*lam12+3*etaphi*lam11
  dlam21dt=dr21dtV+(dkappa1dt+etaphi*kappa1)*lam31+(dkappa2dt+2*etaphi*kappa2)*lam22+4*etaphi*lam21
  dlam31dt=dr31dtV+(dkappa1dt+etaphi*kappa1)*lam41+(dkappa2dt+2*etaphi*kappa2)*lam32+5*etaphi*lam31
  dlam02dt=dr02dtV+(dkappa1dt+etaphi*kappa1)*lam12+(dkappa2dt+2*etaphi*kappa2)*lam03+4*etaphi*lam02
  dlam12dt=dr12dtV+(dkappa1dt+etaphi*kappa1)*lam22+(dkappa2dt+2*etaphi*kappa2)*lam13+5*etaphi*lam12

  dlam60dt=0.Q+00
  dlam70dt=0.Q+00
  dlam41dt=0.Q+00
  dlam51dt=0.Q+00
  dlam22dt=0.Q+00
  dlam32dt=0.Q+00
  dlam03dt=0.Q+00
  dlam13dt=0.Q+00

  dckdt=0.Q+0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dydx=0.Q+00

  dydx(1)=dlam10dt
  dydx(2)=dlam20dt
  dydx(3)=dlam30dt
  dydx(4)=dlam40dt
  dydx(5)=dlam50dt
  dydx(6)=dlam60dt
  dydx(7)=dlam70dt
  dydx(8)=dlam01dt
  dydx(9)=dlam11dt
  dydx(10)=dlam21dt
  dydx(11)=dlam31dt
  dydx(12)=dlam41dt
  dydx(13)=dlam51dt
  dydx(14)=dlam02dt
  dydx(15)=dlam12dt
  dydx(16)=dlam22dt
  dydx(17)=dlam32dt
  dydx(18)=dlam03dt
  dydx(19)=dlam13dt
  dydx(Nv+1)=dlam00dt
  dydx((Nv+1)+1)=dhdt
  dydx((Nv+1)+(Nh+2)+1)=dZ_pidt
  dydx((Nv+1)+(Nh+2)+2)=0.Q+0    !Z_K
  dydx((Nv+1)+(Nh+2)+3)=dZpsidt  !Z_l
  dydx((Nv+1)+(Nh+2)+4)=0.Q+0    !Z_s
  dydx((Nv+1)+(Nh+2)+5)=dZAdt
  dydx((Nv+1)+(Nh+2)+6)=dZphi_p0dt
  dydx((Nv+1)+(Nh+2)+Nz+1)=djldt
  dydx((Nv+1)+(Nh+2)+Nz+2)=djsdt
  dydx((Nv+1)+(Nh+2)+Nz+3)=dSldt
  dydx((Nv+1)+(Nh+2)+Nz+4)=dSsdt
  dydx((Nv+1)+(Nh+2)+Nz+5)=dckdt
  dydx((Nv+1)+(Nh+2)+Nz+Nck+1)=dtg
  dydx((Nv+1)+(Nh+2)+Nz+Nck+2)=dtg3A
  dydx((Nv+1)+(Nh+2)+Nz+Nck+3)=dtg_qbAq_s

  write(*,*)k*hc
  open(unit=51,file='./buffer/k_save.dat',position='append')
    write(51,*)k*hc
  close(51)
  open(unit=51,file='./buffer/g_save.dat',position='append')
    write(51,*)g
  close(51)
  open(unit=51,file='./buffer/lam_save.dat',position='append')
    write(51,*)lam10,lam20,lam30,lam40,lam50,lam01,lam11,lam21,lam31,lam02,lam12
  close(51)
  open(unit=51,file='./buffer/SlSs_save.dat',position='append')
    write(51,*)Sl,Ss
  close(51)
  open(unit=51,file='./buffer/h_save.dat',position='append')
    write(51,*)h
  close(51)
  open(unit=51,file='./buffer/etaphi.dat',position='append')
    write(51,*)etaphi
  close(51)
end
