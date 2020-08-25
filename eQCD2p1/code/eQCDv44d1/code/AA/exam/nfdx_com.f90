subroutine nfdx_com(mf2,T,mu,l,lb,k)
!Calculating the right hand side of differential equations

  implicit none


  real(8) mf2,T,mu,l,lb,k
  real(8) Fnf0,Fnf1,Fnf2
  external Fnf0,Fnf1,Fnf2
  real(8) x,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x
  real(8) nffFd0x,nffFd1x,nffFd2x,nffFd3x,nffFd4x,nffFd5x,nfaFd0x,nfaFd1x,nfaFd2x,nfaFd3x,nfaFd4x,nfaFd5x
  real(8) nffF,nfaF,nffFd1,nfaFd1,nffFd2,nfaFd2,nffFd3,nfaFd3,nffFd4,nfaFd4,nffFd5,nfaFd5
  real(8) finvEF,finvEFd1,finvEFd2

  common /nffFdx_com/ nffFd0x,nffFd1x,nffFd2x,nffFd3x,nffFd4x,nffFd5x,nfaFd0x,nfaFd1x,nfaFd2x,nfaFd3x,nfaFd4x,nfaFd5x
  common /nffFd_com/ nffF,nfaF,nffFd1,nfaFd1,nffFd2,nfaFd2,nffFd3,nfaFd3,nffFd4,nfaFd4,nffFd5,nfaFd5
  common /finvEFd_com/ finvEF,finvEFd1,finvEFd2

  x=k*Sqrt(1 + mf2) - mu
  nf0=Fnf0(x,T,l,lb)
  nf1=Fnf1(x,T,l,lb)
  nf2=Fnf2(x,T,l,lb)
  call nfdx(l,lb,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
  nffFd0x=nfd0x
  nffFd1x=nfd1x
  nffFd2x=nfd2x
  nffFd3x=nfd3x
  nffFd4x=nfd4x
  nffFd5x=nfd5x

  x=k*Sqrt(1 + mf2) + mu
  nf0=Fnf0(x,T,lb,l)
  nf1=Fnf1(x,T,lb,l)
  nf2=Fnf2(x,T,lb,l)
  call nfdx(lb,l,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
  nfaFd0x=nfd0x
  nfaFd1x=nfd1x
  nfaFd2x=nfd2x
  nfaFd3x=nfd3x
  nfaFd4x=nfd4x
  nfaFd5x=nfd5x


  nffF=nffFd0x
  nfaF=nfaFd0x


  nffFd1=(k*nffFd1x)/(2.*Sqrt(1 + mf2))
  nfaFd1=(k*nfaFd1x)/(2.*Sqrt(1 + mf2))

  nffFd2=(k*(-nffFd1x + k*Sqrt(1 + mf2)*nffFd2x))/(4.*(1 + mf2)**1.5)
  nfaFd2=(k*(-nfaFd1x + k*Sqrt(1 + mf2)*nfaFd2x))/(4.*(1 + mf2)**1.5)

  nffFd3=(k*(3*nffFd1x + k*(-3*Sqrt(1 + mf2)*nffFd2x + k*(1 + mf2)*nffFd3x)))/    &
  (8.*(1 + mf2)**2.5)
  nfaFd3=(k*(3*nfaFd1x + k*(-3*Sqrt(1 + mf2)*nfaFd2x + k*(1 + mf2)*nfaFd3x)))/    &
  (8.*(1 + mf2)**2.5)

  nffFd4=(k*(-15*nffFd1x + k*(15*Sqrt(1 + mf2)*nffFd2x +                          &
         k*(1 + mf2)*(-6*nffFd3x + k*Sqrt(1 + mf2)*nffFd4x))))/                   &
  (16.*(1 + mf2)**3.5)
  nfaFd4=(k*(-15*nfaFd1x + k*(15*Sqrt(1 + mf2)*nfaFd2x +                          &
         k*(1 + mf2)*(-6*nfaFd3x + k*Sqrt(1 + mf2)*nfaFd4x))))/                   &
  (16.*(1 + mf2)**3.5)

  nffFd5=(k*(105*nffFd1x + k*(-105*Sqrt(1 + mf2)*nffFd2x +                        &
         k*(1 + mf2)*(45*nffFd3x +                                                &
            k*(-10*Sqrt(1 + mf2)*nffFd4x + k*(1 + mf2)*nffFd5x)))))/              &
  (32.*(1 + mf2)**4.5)
  nfaFd5=(k*(105*nfaFd1x + k*(-105*Sqrt(1 + mf2)*nfaFd2x +                        &
         k*(1 + mf2)*(45*nfaFd3x +                                                &
            k*(-10*Sqrt(1 + mf2)*nfaFd4x + k*(1 + mf2)*nfaFd5x)))))/              &
  (32.*(1 + mf2)**4.5)


  finvEF=1/(k*Sqrt(1 + mf2))
  finvEFd1=-1/(2.*k*(1 + mf2)**1.5)
  finvEFd2=3/(4.*k*(1 + mf2)**2.5)


end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nfqmp_com(mf2,T,mu,l,lb,k,q,costhe)
!Calculating the right hand side of differential equations
!costhe is the cosine angle between the vectors q amd p

  implicit none

  real(8) mf2,T,mu,l,lb,k,q,costhe
  real(8) k2,q2_qmp
  real(8) Fnf0,Fnf1,Fnf2
  external Fnf0,Fnf1,Fnf2
  real(8) x,nf0,nf1,nf2
  real(8) nffqmp,nfaqmp,Eqmp

  common /nffqmp_com/ nffqmp,nfaqmp,Eqmp



  k2=k**2
  q2_qmp=q**2+k2-2.*q*k*costhe
  if(q2_qmp<k2)then
    q2_qmp=k2
  end if
  Eqmp=Sqrt(q2_qmp+k2*mf2)

  x=Eqmp-mu
  nf0=Fnf0(x,T,l,lb)
  nf1=Fnf1(x,T,l,lb)
  nf2=Fnf2(x,T,l,lb)

  nffqmp=nf0 + lb*nf1 + l*nf2


  x=Eqmp+mu
  nf0=Fnf0(x,T,lb,l)
  nf1=Fnf1(x,T,lb,l)
  nf2=Fnf2(x,T,lb,l)

  nfaqmp=nf0 + l*nf1 + lb*nf2

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nbqmp_com(T,q,p,costhe)
!Calculating the right hand side of differential equations
!costhe is the cosine angle between the vectors q amd p

  implicit none

  real(8) T,q,p,costhe
  real(8) qmp2
  real(8) Fnb
  external Fnb
  real(8) x,nb
  real(8) nbqGluon,nbd1xqGluon,nbqmpGluon,nbd1xqmpGluon,qmp

  common /nbqmp_comm/ nbqGluon,nbd1xqGluon,nbqmpGluon,nbd1xqmpGluon,qmp


  qmp2=q**2+p**2-2.*q*p*costhe
  qmp=Sqrt(qmp2)

  x=q
  nb=Fnb(x,T)
  nbqGluon=nb
  nbd1xqGluon=-((nb*(1 + nb))/T)

  x=qmp
  nb=Fnb(x,T)
  nbqmpGluon=nb
  nbd1xqmpGluon=-((nb*(1 + nb))/T)

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nbdx_com(mb2,T,k)
!Calculating the right hand side of differential equations

  implicit none

  real(8) mb2,T,k
  real(8) Fnb
  external Fnb
  real(8) nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x
  real(8) nbBd0x,nbBd1x,nbBd2x,nbBd3x,nbBd4x,nbBd5x
  real(8) nbB,nbBd1,nbBd2,nbBd3,nbBd4,nbBd5
  real(8) finvEB,finvEBd1,finvEBd2

  common /nbBdx_com/ nbBd0x,nbBd1x,nbBd2x,nbBd3x,nbBd4x,nbBd5x
  common /nbBd_com/ nbB,nbBd1,nbBd2,nbBd3,nbBd4,nbBd5
  common /finvEBd_com/ finvEB,finvEBd1,finvEBd2


  nb=Fnb(k*Sqrt(1.d0 + mb2),T)
  call nbdx(nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x)
  nbBd0x=nbd0x
  nbBd1x=nbd1x
  nbBd2x=nbd2x
  nbBd3x=nbd3x
  nbBd4x=nbd4x
  nbBd5x=nbd5x


  nbB=nbBd0x
  nbBd1=(k*nbBd1x)/(2.*Sqrt(1 + mb2))
  nbBd2=(k*(-nbBd1x + k*Sqrt(1 + mb2)*nbBd2x))/(4.*(1 + mb2)**1.5)
  nbBd3=(k*(3*nbBd1x + k*(-3*Sqrt(1 + mb2)*nbBd2x + k*(1 + mb2)*nbBd3x)))/        &
  (8.*(1 + mb2)**2.5)
  nbBd4=(k*(-15*nbBd1x + k*(15*Sqrt(1 + mb2)*nbBd2x +                             &
         k*(1 + mb2)*(-6*nbBd3x + k*Sqrt(1 + mb2)*nbBd4x))))/                     &
  (16.*(1 + mb2)**3.5)
  nbBd5=(k*(105*nbBd1x + k*(-105*Sqrt(1 + mb2)*nbBd2x +                           &
         k*(1 + mb2)*(45*nbBd3x +                                                 &
            k*(-10*Sqrt(1 + mb2)*nbBd4x + k*(1 + mb2)*nbBd5x)))))/                &
  (32.*(1 + mb2)**4.5)


  finvEB=1/(k*Sqrt(1 + mb2))
  finvEBd1=-1/(2.*k*(1 + mb2)**1.5)
  finvEBd2=3/(4.*k*(1 + mb2)**2.5)

end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine nfGhostdx_com(T,k)
!Calculating the right hand side of differential equations

  implicit none

  real(8) T,k
  real(8) Fnf
  real(8) nf,nfd1x,nfd2x,nfd3x,nfd4x
  real(8) nfGhost,nfd1xGhost,nfd2xGhost,nfd3xGhost,nfd4xGhost
  external Fnf

  common /nfGhost_com/ nfGhost,nfd1xGhost,nfd2xGhost,nfd3xGhost,nfd4xGhost

  nf=Fnf(k,T)
  nfd1x=((-1 + nf)*nf)/T
  nfd2x=((-1 + nf)*nf*(-1 + 2*nf))/T**2
  nfd3x=((-1 + nf)*nf*(1 + 6*(-1 + nf)*nf))/T**3
  nfd4x=((-1 + nf)*nf*(-1 + 2*nf)*(1 + 12*(-1 + nf)*nf))/T**4

  nfGhost=nf
  nfd1xGhost=nfd1x
  nfd2xGhost=nfd2x
  nfd3xGhost=nfd3x
  nfd4xGhost=nfd4x

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nfqmpGhost_com(T,q,p,costhe)
!Calculating the right hand side of differential equations
!costhe is the cosine angle between the vectors q amd p

  implicit none

  real(8) T,q,p,costhe
  real(8) qmp2,qmp
  real(8) Fnf
  external Fnf
  real(8) nf
  real(8) nfqmpGhost

  common /nfqmpGhost_comm/ nfqmpGhost


  qmp2=q**2+p**2-2.*q*p*costhe
  qmp=Sqrt(qmp2)


  nf=Fnf(qmp,T)

  nfqmpGhost=nf

end





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine fabmf_com(mb2,mf2,T,k)
!Calculating the right hand side of differential equations

  implicit none


  real(8) mb2,mf2,T,k
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) mu0,p0,p0c

  complex(8) fabmfBFkc,fabmfFBMkc,fabmfBFk,fabmfFBk,fabmfBFMk,fabmfFBMk,fabmfBFMkc,fabmfFBkc,&
             fabmfBFkcd0Bd1F,fabmfFBMkcd0Bd1F,fabmfBFkd0Bd1F,fabmfFBkd0Bd1F,fabmfBFMkd0Bd1F,&
             fabmfFBMkd0Bd1F,fabmfBFMkcd0Bd1F,fabmfFBkcd0Bd1F,fabmfBFkcd1Bd0F,fabmfFBMkcd1Bd0F,&
             fabmfBFkd1Bd0F,fabmfFBkd1Bd0F,fabmfBFMkd1Bd0F,fabmfFBMkd1Bd0F,fabmfBFMkcd1Bd0F,&
             fabmfFBkcd1Bd0F,fabmfBFkd1Bd1F,fabmfBFMkd1Bd1F,fabmfBFkcd1Bd1F,fabmfFBMkcd1Bd1F,&
             fabmfFBkd1Bd1F,fabmfFBMkd1Bd1F,fabmfBFMkcd1Bd1F,fabmfFBkcd1Bd1F,fabmfBFkd0Bd2F,&
             fabmfBFMkd0Bd2F,fabmfBFkcd0Bd2F,fabmfFBkd0Bd2F,fabmfFBMkd0Bd2F,fabmfFBMkcd0Bd2F,&
             fabmfBFMkcd0Bd2F,fabmfFBkcd0Bd2F,fabmfFBkd2Bd0F,fabmfFBMkd2Bd0F,fabmfFBMkcd2Bd0F,&
             fabmfBFkd2Bd0F,fabmfBFMkd2Bd0F,fabmfBFkcd2Bd0F,fabmfBFMkcd2Bd0F,fabmfFBkcd2Bd0F,&
             fabmfBFkcd2Bd1F,fabmfFBMkcd2Bd1F,fabmfBFkd2Bd1F,fabmfFBkd2Bd1F,fabmfBFMkd2Bd1F,fabmfFBMkd2Bd1F,&
             fabmfBFMkcd2Bd1F,fabmfFBkcd2Bd1F,fabmfBFkcd1Bd2F,fabmfFBMkcd1Bd2F,fabmfBFkd1Bd2F,fabmfFBkd1Bd2F,&
             fabmfBFMkd1Bd2F,fabmfFBMkd1Bd2F,fabmfBFMkcd1Bd2F,fabmfFBkcd1Bd2F,&
             fabmfFFk2,fabmfFFMk2,fabmfFFk2c,fabmfFFMk2c,fabmfFFk2d0Bd1F,fabmfFFMk2d0Bd1F,fabmfFFk2cd0Bd1F,&
             fabmfFFMk2cd0Bd1F,fabmfFFk2d1Bd0F,fabmfFFMk2d1Bd0F,fabmfFFk2cd1Bd0F,fabmfFFMk2cd1Bd0F

  common /fabmf_com_com/ fabmfBFkc,fabmfFBMkc,fabmfBFk,fabmfFBk,fabmfBFMk,fabmfFBMk,fabmfBFMkc,fabmfFBkc,&
                         fabmfBFkcd0Bd1F,fabmfFBMkcd0Bd1F,fabmfBFkd0Bd1F,fabmfFBkd0Bd1F,fabmfBFMkd0Bd1F,&
                         fabmfFBMkd0Bd1F,fabmfBFMkcd0Bd1F,fabmfFBkcd0Bd1F,fabmfBFkcd1Bd0F,fabmfFBMkcd1Bd0F,&
                         fabmfBFkd1Bd0F,fabmfFBkd1Bd0F,fabmfBFMkd1Bd0F,fabmfFBMkd1Bd0F,fabmfBFMkcd1Bd0F,&
                         fabmfFBkcd1Bd0F,fabmfBFkd1Bd1F,fabmfBFMkd1Bd1F,fabmfBFkcd1Bd1F,fabmfFBMkcd1Bd1F,&
                         fabmfFBkd1Bd1F,fabmfFBMkd1Bd1F,fabmfBFMkcd1Bd1F,fabmfFBkcd1Bd1F,fabmfBFkd0Bd2F,&
                         fabmfBFMkd0Bd2F,fabmfBFkcd0Bd2F,fabmfFBkd0Bd2F,fabmfFBMkd0Bd2F,fabmfFBMkcd0Bd2F,&
                         fabmfBFMkcd0Bd2F,fabmfFBkcd0Bd2F,fabmfFBkd2Bd0F,fabmfFBMkd2Bd0F,fabmfFBMkcd2Bd0F,&
                         fabmfBFkd2Bd0F,fabmfBFMkd2Bd0F,fabmfBFkcd2Bd0F,fabmfBFMkcd2Bd0F,fabmfFBkcd2Bd0F,&
                         fabmfBFkcd2Bd1F,fabmfFBMkcd2Bd1F,fabmfBFkd2Bd1F,fabmfFBkd2Bd1F,fabmfBFMkd2Bd1F,fabmfFBMkd2Bd1F,&
                         fabmfBFMkcd2Bd1F,fabmfFBkcd2Bd1F,fabmfBFkcd1Bd2F,fabmfFBMkcd1Bd2F,fabmfBFkd1Bd2F,fabmfFBkd1Bd2F,&
                         fabmfBFMkd1Bd2F,fabmfFBMkd1Bd2F,fabmfBFMkcd1Bd2F,fabmfFBkcd1Bd2F,&
                         fabmfFFk2,fabmfFFMk2,fabmfFFk2c,fabmfFFMk2c,fabmfFFk2d0Bd1F,fabmfFFMk2d0Bd1F,fabmfFFk2cd0Bd1F,&
                         fabmfFFMk2cd0Bd1F,fabmfFFk2d1Bd0F,fabmfFFMk2d1Bd0F,fabmfFFk2cd1Bd0F,fabmfFFMk2cd1Bd0F

  common /mu0p0_com/ mu0,p0,p0c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fabmfBFkc=1/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0c)**2)
  fabmfFBMkc=1/(-(k**2*(1 + mb2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)
  fabmfBFk=1/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0)**2)
  fabmfFBk=1/(-(k**2*(1 + mb2)) + (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0)**2)
  fabmfBFMk=1/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0)**2)
  fabmfFBMk=1/(-(k**2*(1 + mb2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)
  fabmfBFMkc=1/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0c)**2)
  fabmfFBkc=1/(-(k**2*(1 + mb2)) + (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c)**2)

  fabmfBFkcd0Bd1F=k**2/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0c)**2)**2
  fabmfFBMkcd0Bd1F=k**2/(-(k**2*(1 + mb2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2
  fabmfBFkd0Bd1F=k**2/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0)**2)**2
  fabmfFBkd0Bd1F=k**2/(-(k**2*(1 + mb2)) + (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0)**2)**2
  fabmfBFMkd0Bd1F=k**2/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0)**2)**2
  fabmfFBMkd0Bd1F=k**2/(-(k**2*(1 + mb2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2
  fabmfBFMkcd0Bd1F=k**2/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0c)**2)**2
  fabmfFBkcd0Bd1F=k**2/(-(k**2*(1 + mb2)) + (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c)**2)**2

  fabmfBFkcd1Bd0F=-((k*(k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0c))/             &
    (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                         &
         (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0c)**2)**2))
  fabmfFBMkcd1Bd0F=(k*(-(k*Sqrt(1 + mf2)) + mu0 - Complex(0,1)*p0c))/           &
  (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2)
  fabmfBFkd1Bd0F=-((k*(k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0))/               &
    (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                         &
         (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0)**2)**2))
  fabmfFBkd1Bd0F=-((k*(k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0))/               &
    (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                         &
         (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0)**2)**2))
  fabmfBFMkd1Bd0F=(k*(-(k*Sqrt(1 + mb2)) + mu0 - Complex(0,1)*p0))/             &
  (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0)**2)**2)
  fabmfFBMkd1Bd0F=(k*(-(k*Sqrt(1 + mf2)) + mu0 - Complex(0,1)*p0))/             &
  (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)
  fabmfBFMkcd1Bd0F=(k*(-(k*Sqrt(1 + mb2)) + mu0 - Complex(0,1)*p0c))/           &
  (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0c)**2)**2)
  fabmfFBkcd1Bd0F=-((k*(k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c))/             &
    (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                         &
         (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c)**2)**2))

  fabmfBFkd1Bd1F=(-2*k**3*(k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0))/           &
  (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0)**2)**3)
  fabmfBFMkd1Bd1F=(2*k**3*(-(k*Sqrt(1 + mb2)) + mu0 - Complex(0,1)*p0))/        &
  (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0)**2)**3)
  fabmfBFkcd1Bd1F=(-2*k**3*(k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0c))/         &
  (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0c)**2)**3)
  fabmfFBMkcd1Bd1F=(2*k**3*(-(k*Sqrt(1 + mf2)) + mu0 - Complex(0,1)*p0c))/      &
  (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3)
  fabmfFBkd1Bd1F=(-2*k**3*(k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0))/           &
  (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                           &
       (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0)**2)**3)
  fabmfFBMkd1Bd1F=(2*k**3*(-(k*Sqrt(1 + mf2)) + mu0 - Complex(0,1)*p0))/        &
  (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3)
  fabmfBFMkcd1Bd1F=(2*k**3*(-(k*Sqrt(1 + mb2)) + mu0 - Complex(0,1)*p0c))/      &
  (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0c)**2)**3)
  fabmfFBkcd1Bd1F=(-2*k**3*(k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c))/         &
  (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                           &
       (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c)**2)**3)

  fabmfBFkd0Bd2F=(2*k**4)/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0)**2)** &
   3
  fabmfBFMkd0Bd2F=(2*k**4)/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0)**2)** &
   3
  fabmfBFkcd0Bd2F=(2*k**4)/(-(k**2*(1 + mf2)) +                                 &
     (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0c)**2)**3
  fabmfFBkd0Bd2F=(2*k**4)/(-(k**2*(1 + mb2)) + (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0)**2)** &
   3
  fabmfFBMkd0Bd2F=(2*k**4)/(-(k**2*(1 + mb2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)** &
   3
  fabmfFBMkcd0Bd2F=(2*k**4)/(-(k**2*(1 + mb2)) +                                &
     (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3
  fabmfBFMkcd0Bd2F=(2*k**4)/(-(k**2*(1 + mf2)) +                                &
     (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0c)**2)**3
  fabmfFBkcd0Bd2F=(2*k**4)/(-(k**2*(1 + mb2)) +                                 &
     (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c)**2)**3


  fabmfFBkd2Bd0F=(k*(4*k**3*(1 + mf2)**1.5 + k**2*(8 - mb2 + 9*mf2)*            &
       (mu0 - Complex(0,1)*p0) +                                                &
      6*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0)**2 +                            &
      (mu0 - Complex(0,1)*p0)**3))/                                             &
  (2.*(1 + mf2)**1.5*(k**2*(-mb2 + mf2) +                                       &
       2*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0) +                              &
       (mu0 - Complex(0,1)*p0)**2)**3)
  fabmfFBMkd2Bd0F=-(k*(-4*k**3*(1 + mf2)**1.5 + k**2*(8 - mb2 + 9*mf2)*         &
        (mu0 - Complex(0,1)*p0) -                                               &
       6*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0)**2 +                           &
       (mu0 - Complex(0,1)*p0)**3))/                                            &
  (2.*(1 + mf2)**1.5*(k**2*(-mb2 + mf2) -                                       &
       2*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0) +                              &
       (mu0 - Complex(0,1)*p0)**2)**3)
  fabmfFBMkcd2Bd0F=-(k*(-4*k**3*(1 + mf2)**1.5 + k**2*(8 - mb2 + 9*mf2)*        &
        (mu0 - Complex(0,1)*p0c) -                                              &
       6*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0c)**2 +                          &
       (mu0 - Complex(0,1)*p0c)**3))/                                           &
  (2.*(1 + mf2)**1.5*(k**2*(-mb2 + mf2) -                                       &
       2*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0c) +                             &
       (mu0 - Complex(0,1)*p0c)**2)**3)
  fabmfBFkd2Bd0F=(k*(4*k**3*(1 + mb2)**1.5 + k**2*(8 + 9*mb2 - mf2)*            &
       (mu0 - Complex(0,1)*p0) +                                                &
      6*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0)**2 +                            &
      (mu0 - Complex(0,1)*p0)**3))/                                             &
  (2.*(1 + mb2)**1.5*(k**2*(mb2 - mf2) +                                        &
       2*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0) +                              &
       (mu0 - Complex(0,1)*p0)**2)**3)
  fabmfBFMkd2Bd0F=-(k*(-4*k**3*(1 + mb2)**1.5 + k**2*(8 + 9*mb2 - mf2)*         &
        (mu0 - Complex(0,1)*p0) -                                               &
       6*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0)**2 +                           &
       (mu0 - Complex(0,1)*p0)**3))/                                            &
  (2.*(1 + mb2)**1.5*(k**2*(mb2 - mf2) -                                        &
       2*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0) +                              &
       (mu0 - Complex(0,1)*p0)**2)**3)
  fabmfBFkcd2Bd0F=(k*(4*k**3*(1 + mb2)**1.5 + k**2*(8 + 9*mb2 - mf2)*           &
       (mu0 - Complex(0,1)*p0c) +                                               &
      6*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0c)**2 +                           &
      (mu0 - Complex(0,1)*p0c)**3))/                                            &
  (2.*(1 + mb2)**1.5*(k**2*(mb2 - mf2) +                                        &
       2*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0c) +                             &
       (mu0 - Complex(0,1)*p0c)**2)**3)
  fabmfBFMkcd2Bd0F=-(k*(-4*k**3*(1 + mb2)**1.5 + k**2*(8 + 9*mb2 - mf2)*        &
        (mu0 - Complex(0,1)*p0c) -                                              &
       6*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0c)**2 +                          &
       (mu0 - Complex(0,1)*p0c)**3))/                                           &
  (2.*(1 + mb2)**1.5*(k**2*(mb2 - mf2) -                                        &
       2*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0c) +                             &
       (mu0 - Complex(0,1)*p0c)**2)**3)
  fabmfFBkcd2Bd0F=(k*(4*k**3*(1 + mf2)**1.5 + k**2*(8 - mb2 + 9*mf2)*           &
       (mu0 - Complex(0,1)*p0c) +                                               &
      6*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0c)**2 +                           &
      (mu0 - Complex(0,1)*p0c)**3))/                                            &
  (2.*(1 + mf2)**1.5*(k**2*(-mb2 + mf2) +                                       &
       2*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0c) +                             &
       (mu0 - Complex(0,1)*p0c)**2)**3)

  fabmfBFkcd2Bd1F=(k**3*(6*k**3*(1 + mb2)**1.5 +                                &
      k**2*(12 + 13*mb2 - mf2)*(mu0 - Complex(0,1)*p0c) +                       &
      8*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0c)**2 +                           &
      (mu0 - Complex(0,1)*p0c)**3))/                                            &
  ((1 + mb2)**1.5*(k**2*(mb2 - mf2) +                                           &
       2*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0c) +                             &
       (mu0 - Complex(0,1)*p0c)**2)**4)
  fabmfFBMkcd2Bd1F=-((k**3*(-6*k**3*(1 + mf2)**1.5 +                            &
        k**2*(12 - mb2 + 13*mf2)*(mu0 - Complex(0,1)*p0c) -                     &
        8*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0c)**2 +                         &
        (mu0 - Complex(0,1)*p0c)**3))/                                          &
    ((1 + mf2)**1.5*(k**2*(-mb2 + mf2) -                                        &
         2*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0c) +                           &
         (mu0 - Complex(0,1)*p0c)**2)**4))
  fabmfBFkd2Bd1F=(k**3*(6*k**3*(1 + mb2)**1.5 +                                 &
      k**2*(12 + 13*mb2 - mf2)*(mu0 - Complex(0,1)*p0) +                        &
      8*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0)**2 +                            &
      (mu0 - Complex(0,1)*p0)**3))/                                             &
  ((1 + mb2)**1.5*(k**2*(mb2 - mf2) +                                           &
       2*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0) +                              &
       (mu0 - Complex(0,1)*p0)**2)**4)
  fabmfFBkd2Bd1F=(k**3*(6*k**3*(1 + mf2)**1.5 +                                 &
      k**2*(12 - mb2 + 13*mf2)*(mu0 - Complex(0,1)*p0) +                        &
      8*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0)**2 +                            &
      (mu0 - Complex(0,1)*p0)**3))/                                             &
  ((1 + mf2)**1.5*(k**2*(-mb2 + mf2) +                                          &
       2*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0) +                              &
       (mu0 - Complex(0,1)*p0)**2)**4)
  fabmfBFMkd2Bd1F=-((k**3*(-6*k**3*(1 + mb2)**1.5 +                             &
        k**2*(12 + 13*mb2 - mf2)*(mu0 - Complex(0,1)*p0) -                      &
        8*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0)**2 +                          &
        (mu0 - Complex(0,1)*p0)**3))/                                           &
    ((1 + mb2)**1.5*(k**2*(mb2 - mf2) -                                         &
         2*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0) +                            &
         (mu0 - Complex(0,1)*p0)**2)**4))
  fabmfFBMkd2Bd1F=-((k**3*(-6*k**3*(1 + mf2)**1.5 +                             &
        k**2*(12 - mb2 + 13*mf2)*(mu0 - Complex(0,1)*p0) -                      &
        8*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0)**2 +                          &
        (mu0 - Complex(0,1)*p0)**3))/                                           &
    ((1 + mf2)**1.5*(k**2*(-mb2 + mf2) -                                        &
         2*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0) +                            &
         (mu0 - Complex(0,1)*p0)**2)**4))
  fabmfBFMkcd2Bd1F=-((k**3*(-6*k**3*(1 + mb2)**1.5 +                            &
        k**2*(12 + 13*mb2 - mf2)*(mu0 - Complex(0,1)*p0c) -                     &
        8*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0c)**2 +                         &
        (mu0 - Complex(0,1)*p0c)**3))/                                          &
    ((1 + mb2)**1.5*(k**2*(mb2 - mf2) -                                         &
         2*k*Sqrt(1 + mb2)*(mu0 - Complex(0,1)*p0c) +                           &
         (mu0 - Complex(0,1)*p0c)**2)**4))
  fabmfFBkcd2Bd1F=(k**3*(6*k**3*(1 + mf2)**1.5 +                                &
      k**2*(12 - mb2 + 13*mf2)*(mu0 - Complex(0,1)*p0c) +                       &
      8*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0c)**2 +                           &
      (mu0 - Complex(0,1)*p0c)**3))/                                            &
  ((1 + mf2)**1.5*(k**2*(-mb2 + mf2) +                                          &
       2*k*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0c) +                             &
       (mu0 - Complex(0,1)*p0c)**2)**4)

  fabmfBFkcd1Bd2F=(-6*k**5*(k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0c))/         &
  (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0c)**2)**4)
  fabmfFBMkcd1Bd2F=(6*k**5*(-(k*Sqrt(1 + mf2)) + mu0 - Complex(0,1)*p0c))/      &
  (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**4)
  fabmfBFkd1Bd2F=(-6*k**5*(k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0))/           &
  (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0)**2)**4)
  fabmfFBkd1Bd2F=(-6*k**5*(k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0))/           &
  (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                           &
       (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0)**2)**4)
  fabmfBFMkd1Bd2F=(6*k**5*(-(k*Sqrt(1 + mb2)) + mu0 - Complex(0,1)*p0))/        &
  (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0)**2)**4)
  fabmfFBMkd1Bd2F=(6*k**5*(-(k*Sqrt(1 + mf2)) + mu0 - Complex(0,1)*p0))/        &
  (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**4)
  fabmfBFMkcd1Bd2F=(6*k**5*(-(k*Sqrt(1 + mb2)) + mu0 - Complex(0,1)*p0c))/      &
  (Sqrt(1 + mb2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mb2) - mu0 + Complex(0,1)*p0c)**2)**4)
  fabmfFBkcd1Bd2F=(-6*k**5*(k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c))/         &
  (Sqrt(1 + mf2)*(-(k**2*(1 + mb2)) +                                           &
       (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c)**2)**4)

  fabmfFFk2=1/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mf2) + 2*mu0 - Complex(0,2)*p0)**2)
  fabmfFFMk2=1/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mf2) - 2*mu0 + Complex(0,2)*p0)**2)
  fabmfFFk2c=1/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mf2) + 2*mu0 - Complex(0,2)*p0c)**2)
  fabmfFFMk2c=1/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mf2) - 2*mu0 + Complex(0,2)*p0c)**2)

  fabmfFFk2d0Bd1F=k**2/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mf2) + 2*mu0 - Complex(0,2)*p0)**2)**2
  fabmfFFMk2d0Bd1F=k**2/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mf2) - 2*mu0 + Complex(0,2)*p0)**2)**2
  fabmfFFk2cd0Bd1F=k**2/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mf2) + 2*mu0 - Complex(0,2)*p0c)** &
      2)**2
  fabmfFFMk2cd0Bd1F=k**2/(-(k**2*(1 + mf2)) + (k*Sqrt(1 + mf2) - 2*mu0 + Complex(0,2)*p0c)** &
      2)**2

  fabmfFFk2d1Bd0F=-(k*(k*Sqrt(1 + mf2) + 2*mu0 - Complex(0,2)*p0))/             &
  (16.*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0)**2*                                &
    (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0)**2)
  fabmfFFMk2d1Bd0F=-(k*(k*Sqrt(1 + mf2) - 2*mu0 + Complex(0,2)*p0))/            &
  (16.*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0)**2*                                &
    (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)
  fabmfFFk2cd1Bd0F=-(k*(k*Sqrt(1 + mf2) + 2*mu0 - Complex(0,2)*p0c))/           &
  (16.*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0c)**2*                               &
    (k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c)**2)
  fabmfFFMk2cd1Bd0F=-(k*(k*Sqrt(1 + mf2) - 2*mu0 + Complex(0,2)*p0c))/          &
  (16.*Sqrt(1 + mf2)*(mu0 - Complex(0,1)*p0c)**2*                               &
    (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end

subroutine fm_com(mb2,mf2,T,k)
!Calculating the right hand side of differential equations

  implicit none


  real(8) mb2,mf2,T,k
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) mu0,p0,p0c

  real(8) finvEB,finvEBd1,finvEBd2
  real(8) finvEF,finvEFd1,finvEFd2
  complex(8) fm1Bk,fm1BMk,fm1Fk,fm1FMk,fm1Bkc,fm1FMkc,fm1BMkc,fm1Fkc
  real(8) fm1Bkd1,fm1BMkd1,fm1Bkcd1,fm1Fkd1,fm1FMkd1,fm1FMkcd1,fm1BMkcd1,fm1Fkcd1,&
          fm1Bkd2,fm1BMkd2,fm1Bkcd2,fm1BMkcd2,fm1Fkd2,fm1FMkd2,fm1Fkcd2,fm1FMkcd2
  complex(8) fm2Bk,fm2BMk,fm2Fk,fm2FMk,fm2Bkc,fm2FMkc,fm2BMkc,fm2Fkc,fm2Bkd1,fm2BMkd1,fm2Bkcd1,fm2Fkd1,&
             fm2FMkd1,fm2FMkcd1,fm2BMkcd1,fm2Fkcd1,fm2Bkd2,fm2BMkd2,fm2Bkcd2,fm2BMkcd2,fm2Fkd2,fm2FMkd2,&
             fm2Fkcd2,fm2FMkcd2
  complex(8) fm1Fk2,fm1FMk2,fm1Fk2c,fm1FMk2c
  real(8) fm1Fk2d1,fm1FMk2d1,fm1Fk2cd1,fm1FMk2cd1
  complex(8) fn1Bk,fn1Bkc,fn1Bkd1,fn1Bkcd1,fn1Bkd2,fn1Bkcd2
  complex(8) fpm1Fk,fpm1FMk,fpm1Fkc,fpm1Fkd1,fpm1FMkd1,fpm1Fkcd1
  complex(8) fpm2aBk,fpm2aBkc,fpm2aBkd1,fpm2aBkcd1,fpm2aBkd2,fpm2aBkcd2
  complex(8) fpm2bFk,fpm2bFMk,fpm2bFkc,fpm2bFMkc,fpm2bFkd1,fpm2bFMkd1,fpm2bFMkcd1

  common /finvEBd_com/ finvEB,finvEBd1,finvEBd2
  common /finvEFd_com/ finvEF,finvEFd1,finvEFd2
  common /fm_com_com/ fm1Bk,fm1BMk,fm1Fk,fm1FMk,fm1Bkc,fm1FMkc,fm1BMkc,fm1Fkc,&
                      fm1Bkd1,fm1BMkd1,fm1Bkcd1,fm1Fkd1,fm1FMkd1,fm1FMkcd1,fm1BMkcd1,fm1Fkcd1,&
                      fm1Bkd2,fm1BMkd2,fm1Bkcd2,fm1BMkcd2,fm1Fkd2,fm1FMkd2,fm1Fkcd2,fm1FMkcd2,&
                      fm2Bk,fm2BMk,fm2Fk,fm2FMk,fm2Bkc,fm2FMkc,fm2BMkc,fm2Fkc,fm2Bkd1,fm2BMkd1,fm2Bkcd1,fm2Fkd1,&
                      fm2FMkd1,fm2FMkcd1,fm2BMkcd1,fm2Fkcd1,fm2Bkd2,fm2BMkd2,fm2Bkcd2,fm2BMkcd2,fm2Fkd2,fm2FMkd2,&
                      fm2Fkcd2,fm2FMkcd2,&
                      fm1Fk2,fm1FMk2,fm1Fk2c,fm1FMk2c,fm1Fk2d1,fm1FMk2d1,fm1Fk2cd1,fm1FMk2cd1,&
                      fn1Bk,fn1Bkc,fn1Bkd1,fn1Bkcd1,fn1Bkd2,fn1Bkcd2,&
                      fpm1Fk,fpm1FMk,fpm1Fkc,fpm1Fkd1,fpm1FMkd1,fpm1Fkcd1,&
                      fpm2aBk,fpm2aBkc,fpm2aBkd1,fpm2aBkcd1,fpm2aBkd2,fpm2aBkcd2,&
                      fpm2bFk,fpm2bFMk,fpm2bFkc,fpm2bFMkc,fpm2bFkd1,fpm2bFMkd1,fpm2bFMkcd1

  common /mu0p0_com/ mu0,p0,p0c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  fm1Bk=k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0
  fm1BMk=-(k*Sqrt(1 + mb2)) + mu0 - Complex(0,1)*p0
  fm1Fk=k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0
  fm1FMk=-(k*Sqrt(1 + mf2)) + mu0 - Complex(0,1)*p0
  fm1Bkc=k*Sqrt(1 + mb2) + mu0 - Complex(0,1)*p0c
  fm1FMkc=-(k*Sqrt(1 + mf2)) + mu0 - Complex(0,1)*p0c
  fm1BMkc=-(k*Sqrt(1 + mb2)) + mu0 - Complex(0,1)*p0c
  fm1Fkc=k*Sqrt(1 + mf2) + mu0 - Complex(0,1)*p0c

  fm1Bkd1=k/(2.*Sqrt(1 + mb2))
  fm1BMkd1=-k/(2.*Sqrt(1 + mb2))
  fm1Bkcd1=k/(2.*Sqrt(1 + mb2))
  fm1Fkd1=k/(2.*Sqrt(1 + mf2))
  fm1FMkd1=-k/(2.*Sqrt(1 + mf2))
  fm1FMkcd1=-k/(2.*Sqrt(1 + mf2))
  fm1BMkcd1=-k/(2.*Sqrt(1 + mb2))
  fm1Fkcd1=k/(2.*Sqrt(1 + mf2))

  fm1Bkd2=-k/(4.*(1 + mb2)**1.5)
  fm1BMkd2=k/(4.*(1 + mb2)**1.5)
  fm1Bkcd2=-k/(4.*(1 + mb2)**1.5)
  fm1BMkcd2=k/(4.*(1 + mb2)**1.5)
  fm1Fkd2=-k/(4.*(1 + mf2)**1.5)
  fm1FMkd2=k/(4.*(1 + mf2)**1.5)
  fm1Fkcd2=-k/(4.*(1 + mf2)**1.5)
  fm1FMkcd2=k/(4.*(1 + mf2)**1.5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fm2Bk=fm1Bk**2*k*Sqrt(1 + mb2)
  fm2BMk=-(fm1BMk**2*k*Sqrt(1 + mb2))
  fm2Fk=fm1Fk**2*k*Sqrt(1 + mf2)
  fm2FMk=-(fm1FMk**2*k*Sqrt(1 + mf2))
  fm2Bkc=fm1Bkc**2*k*Sqrt(1 + mb2)
  fm2FMkc=-(fm1FMkc**2*k*Sqrt(1 + mf2))
  fm2BMkc=-(fm1BMkc**2*k*Sqrt(1 + mb2))
  fm2Fkc=fm1Fkc**2*k*Sqrt(1 + mf2)

  fm2Bkd1=(fm1Bk*k*(fm1Bk + 4*fm1Bkd1*(1 + mb2)))/(2.*Sqrt(1 + mb2))
  fm2BMkd1=-(fm1BMk*k*(fm1BMk + 4*fm1BMkd1*(1 + mb2)))/(2.*Sqrt(1 + mb2))
  fm2Bkcd1=(fm1Bkc*k*(fm1Bkc + 4*fm1Bkcd1*(1 + mb2)))/(2.*Sqrt(1 + mb2))
  fm2Fkd1=(fm1Fk*k*(fm1Fk + 4*fm1Fkd1*(1 + mf2)))/(2.*Sqrt(1 + mf2))
  fm2FMkd1=-(fm1FMk*k*(fm1FMk + 4*fm1FMkd1*(1 + mf2)))/(2.*Sqrt(1 + mf2))
  fm2FMkcd1=-(fm1FMkc*k*(fm1FMkc + 4*fm1FMkcd1*(1 + mf2)))/(2.*Sqrt(1 + mf2))
  fm2BMkcd1=-(fm1BMkc*k*(fm1BMkc + 4*fm1BMkcd1*(1 + mb2)))/(2.*Sqrt(1 + mb2))
  fm2Fkcd1=(fm1Fkc*k*(fm1Fkc + 4*fm1Fkcd1*(1 + mf2)))/(2.*Sqrt(1 + mf2))

  fm2Bkd2=(k*(-fm1Bk**2 + 8*fm1Bkd1**2*(1 + mb2)**2 +                           &
      8*fm1Bk*(1 + mb2)*(fm1Bkd1 + fm1Bkd2 + fm1Bkd2*mb2)))/                    &
  (4.*(1 + mb2)**1.5)
  fm2BMkd2=-(k*(-fm1BMk**2 + 8*fm1BMkd1**2*(1 + mb2)**2 +                       &
       8*fm1BMk*(1 + mb2)*(fm1BMkd1 + fm1BMkd2 + fm1BMkd2*mb2)))/               &
  (4.*(1 + mb2)**1.5)
  fm2Bkcd2=(k*(-fm1Bkc**2 + 8*fm1Bkcd1**2*(1 + mb2)**2 +                        &
      8*fm1Bkc*(1 + mb2)*(fm1Bkcd1 + fm1Bkcd2 + fm1Bkcd2*mb2)))/                &
  (4.*(1 + mb2)**1.5)
  fm2BMkcd2=-(k*(-fm1BMkc**2 + 8*fm1BMkcd1**2*(1 + mb2)**2 +                    &
       8*fm1BMkc*(1 + mb2)*(fm1BMkcd1 + fm1BMkcd2 + fm1BMkcd2*mb2)))/           &
  (4.*(1 + mb2)**1.5)
  fm2Fkd2=(k*(-fm1Fk**2 + 8*fm1Fkd1**2*(1 + mf2)**2 +                           &
      8*fm1Fk*(1 + mf2)*(fm1Fkd1 + fm1Fkd2 + fm1Fkd2*mf2)))/                    &
  (4.*(1 + mf2)**1.5)
  fm2FMkd2=-(k*(-fm1FMk**2 + 8*fm1FMkd1**2*(1 + mf2)**2 +                       &
       8*fm1FMk*(1 + mf2)*(fm1FMkd1 + fm1FMkd2 + fm1FMkd2*mf2)))/               &
  (4.*(1 + mf2)**1.5)
  fm2Fkcd2=(k*(-fm1Fkc**2 + 8*fm1Fkcd1**2*(1 + mf2)**2 +                        &
      8*fm1Fkc*(1 + mf2)*(fm1Fkcd1 + fm1Fkcd2 + fm1Fkcd2*mf2)))/                &
  (4.*(1 + mf2)**1.5)
  fm2FMkcd2=-(k*(-fm1FMkc**2 + 8*fm1FMkcd1**2*(1 + mf2)**2 +                    &
       8*fm1FMkc*(1 + mf2)*(fm1FMkcd1 + fm1FMkcd2 + fm1FMkcd2*mf2)))/           &
  (4.*(1 + mf2)**1.5)

  fm1Fk2=k*Sqrt(1 + mf2) + 2*mu0 - Complex(0,2)*p0
  fm1FMk2=-(k*Sqrt(1 + mf2)) + 2*mu0 - Complex(0,2)*p0
  fm1Fk2c=k*Sqrt(1 + mf2) + 2*mu0 - Complex(0,2)*p0c
  fm1FMk2c=-(k*Sqrt(1 + mf2)) + 2*mu0 - Complex(0,2)*p0c

  fm1Fk2d1=k/(2.*Sqrt(1 + mf2))
  fm1FMk2d1=-k/(2.*Sqrt(1 + mf2))
  fm1Fk2cd1=k/(2.*Sqrt(1 + mf2))
  fm1FMk2cd1=-k/(2.*Sqrt(1 + mf2))

  fn1Bk=finvEB*fm1Bk*fm1BMk
  fn1Bkc=finvEB*fm1Bkc*fm1BMkc

  fn1Bkd1=finvEBd1*fm1Bk*fm1BMk + finvEB*fm1Bkd1*fm1BMk + finvEB*fm1Bk*fm1BMkd1
  fn1Bkcd1=finvEBd1*fm1Bkc*fm1BMkc + finvEB*fm1Bkcd1*fm1BMkc + finvEB*fm1Bkc*fm1BMkcd1

  fn1Bkd2=finvEBd2*fm1Bk*fm1BMk + 2*finvEBd1*(fm1Bkd1*fm1BMk + fm1Bk*fm1BMkd1) + &
  finvEB*(fm1Bkd2*fm1BMk + 2*fm1Bkd1*fm1BMkd1 + fm1Bk*fm1BMkd2)
  fn1Bkcd2=finvEBd2*fm1Bkc*fm1BMkc + 2*finvEBd1*(fm1Bkcd1*fm1BMkc + fm1Bkc*fm1BMkcd1) + &
  finvEB*(fm1Bkcd2*fm1BMkc + 2*fm1Bkcd1*fm1BMkcd1 + fm1Bkc*fm1BMkcd2)

  fpm1Fk=finvEF*fm1Fk*fm1Fk2
  fpm1FMk=finvEF*fm1FMk*fm1FMk2
  fpm1Fkc=finvEF*fm1Fk2c*fm1Fkc

  fpm1Fkd1=finvEFd1*fm1Fk*fm1Fk2 + finvEF*fm1Fk*fm1Fk2d1 + finvEF*fm1Fk2*fm1Fkd1
  fpm1FMkd1=finvEFd1*fm1FMk*fm1FMk2 + finvEF*fm1FMk*fm1FMk2d1 + finvEF*fm1FMk2*fm1FMkd1
  fpm1Fkcd1=finvEFd1*fm1Fk2c*fm1Fkc + finvEF*fm1Fk2cd1*fm1Fkc + finvEF*fm1Fk2c*fm1Fkcd1

  fpm2aBk=fm1Bk*fm1BMk*k*Sqrt(1 + mb2)
  fpm2aBkc=fm1Bkc*fm1BMkc*k*Sqrt(1 + mb2)

  fpm2aBkd1=(k*(2*fm1Bkd1*fm1BMk*(1 + mb2) + fm1Bk*(fm1BMk + 2*fm1BMkd1*(1 + mb2))))/  &
  (2.*Sqrt(1 + mb2))
  fpm2aBkcd1=(k*(2*fm1Bkcd1*fm1BMkc*(1 + mb2) +                                 &
      fm1Bkc*(fm1BMkc + 2*fm1BMkcd1*(1 + mb2))))/(2.*Sqrt(1 + mb2))

  fpm2aBkd2=(k*(fm1Bk*(-fm1BMk + 4*(1 + mb2)*(fm1BMkd1 + fm1BMkd2 + fm1BMkd2*mb2)) +  &
      4*(1 + mb2)*(fm1Bkd2*fm1BMk*(1 + mb2) +                                   &
         fm1Bkd1*(fm1BMk + 2*fm1BMkd1*(1 + mb2)))))/(4.*(1 + mb2)**1.5)
  fpm2aBkcd2=(k*(fm1Bkc*(-fm1BMkc + 4*(1 + mb2)*                                &
          (fm1BMkcd1 + fm1BMkcd2 + fm1BMkcd2*mb2)) +                            &
      4*(1 + mb2)*(fm1Bkcd2*fm1BMkc*(1 + mb2) +                                 &
         fm1Bkcd1*(fm1BMkc + 2*fm1BMkcd1*(1 + mb2)))))/(4.*(1 + mb2)**1.5)

  fpm2bFk=fm1Fk**2*fm1Fk2
  fpm2bFMk=fm1FMk**2*fm1FMk2
  fpm2bFkc=fm1Fk2c*fm1Fkc**2
  fpm2bFMkc=fm1FMk2c*fm1FMkc**2

  fpm2bFkd1=fm1Fk*(fm1Fk*fm1Fk2d1 + 2*fm1Fk2*fm1Fkd1)
  fpm2bFMkd1=fm1FMk*(fm1FMk*fm1FMk2d1 + 2*fm1FMk2*fm1FMkd1)
  fpm2bFMkcd1=fm1FMkc*(fm1FMk2cd1*fm1FMkc + 2*fm1FMk2c*fm1FMkcd1)


end



