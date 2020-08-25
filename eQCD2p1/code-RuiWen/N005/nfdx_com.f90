subroutine nfdx_com(mf2,T,mu,l,lb,k)
!Calculating the right hand side of differential equations

  implicit none

  real(16) mf2,T,mu,l,lb,k
  real(16) x,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x
  real(16) nffFd0x,nffFd1x,nffFd2x,nffFd3x,nffFd4x,nffFd5x,nfaFd0x,nfaFd1x,nfaFd2x,nfaFd3x,nfaFd4x,nfaFd5x

  common /nffFdx_com/ nffFd0x,nffFd1x,nffFd2x,nffFd3x,nffFd4x,nffFd5x,nfaFd0x,nfaFd1x,nfaFd2x,nfaFd3x,nfaFd4x,nfaFd5x

  x=k*Sqrt(1 + mf2) - mu
  call Fnf012(x,T,l,lb,nf0,nf1,nf2)
  call nfdx(l,lb,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
  nffFd0x=nfd0x
  nffFd1x=nfd1x
  nffFd2x=nfd2x
  nffFd3x=nfd3x
  nffFd4x=nfd4x
  nffFd5x=nfd5x

  x=k*Sqrt(1 + mf2) + mu
  call Fnf012(x,T,lb,l,nf0,nf1,nf2)
  call nfdx(lb,l,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
  nfaFd0x=nfd0x
  nfaFd1x=nfd1x
  nfaFd2x=nfd2x
  nfaFd3x=nfd3x
  nfaFd4x=nfd4x
  nfaFd5x=nfd5x

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nfqmp_com(mf2,T,mu,l,lb,k,q,costhe)
!Calculating the right hand side of differential equations
!costhe is the cosine angle between the vectors q amd p

  implicit none

  real(16) mf2,T,mu,l,lb,k,q,costhe
  real(16) k2,q2_qmp
  real(16) x,nf0,nf1,nf2
  real(16) nffqmp,nfaqmp,Eqmp
  real(16) nffqmp_save,nfaqmp_save,Eqmp_save
  integer k_same_common

  common /nffqmp_com/ nffqmp,nfaqmp,Eqmp
  common /nffqmp_com_save/ nffqmp_save,nfaqmp_save,Eqmp_save
  common /k_same_com/k_same_common

  k2=k**2
  q2_qmp=q**2+k2-2.Q+00*q*k*costhe
  if(q2_qmp<k2)then
    if(k_same_common==0)then
      Eqmp=k*Sqrt(1.Q+0+mf2)
      x=Eqmp-mu
      call Fnf012(x,T,l,lb,nf0,nf1,nf2)
      nffqmp=nf0 + lb*nf1 + l*nf2
      x=Eqmp+mu
      call Fnf012(x,T,lb,l,nf0,nf1,nf2)
      nfaqmp=nf0 + l*nf1 + lb*nf2
      nffqmp_save=nffqmp
      nfaqmp_save=nfaqmp
      Eqmp_save=Eqmp
      k_same_common=1
    else
      nffqmp=nffqmp_save
      nfaqmp=nfaqmp_save
      Eqmp=Eqmp_save
    end if
  else
    Eqmp=Sqrt(q2_qmp+k2*mf2)

    x=Eqmp-mu
    call Fnf012(x,T,l,lb,nf0,nf1,nf2)
    nffqmp=nf0 + lb*nf1 + l*nf2

    x=Eqmp+mu
    call Fnf012(x,T,lb,l,nf0,nf1,nf2)
    nfaqmp=nf0 + l*nf1 + lb*nf2
  end if
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nbqmp_com(T,q,p,costhe)
!Calculating the right hand side of differential equations
!costhe is the cosine angle between the vectors q amd p

  implicit none

  real(16) T,q,p,costhe
  real(16) qmp2
  real(16) x,nb
  real(16) nbqGluon,nbd1xqGluon,nbqmpGluon,nbd1xqmpGluon,qmp

  common /nbqmp_comm/ nbqGluon,nbd1xqGluon,nbqmpGluon,nbd1xqmpGluon,qmp

  qmp2=q**2+p**2-2.*q*p*costhe
  qmp=Sqrt(qmp2)

  x=q
  call Fnb(x,T,nb)
  nbqGluon=nb
  nbd1xqGluon=-((nb*(1 + nb))/T)

  x=qmp
  call Fnb(x,T,nb)
  nbqmpGluon=nb
  nbd1xqmpGluon=-((nb*(1 + nb))/T)

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nbdx_com(mb2,T,k)
!Calculating the right hand side of differential equations

  implicit none

  real(16) mb2,T,k
  real(16) nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x
  real(16) nbBd0x,nbBd1x,nbBd2x,nbBd3x,nbBd4x,nbBd5x
  real(16) nbB,nbBd1,nbBd2,nbBd3,nbBd4,nbBd5
  real(16) finvEB,finvEBd1,finvEBd2

  common /nbBdx_com/ nbBd0x,nbBd1x,nbBd2x,nbBd3x,nbBd4x,nbBd5x
  common /nbBd_com/ nbB,nbBd1,nbBd2,nbBd3,nbBd4,nbBd5
  common /finvEBd_com/ finvEB,finvEBd1,finvEBd2


  call Fnb(k*Sqrt(1.Q+0 + mb2),T,nb)
  call nbdx(T,nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x)
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
