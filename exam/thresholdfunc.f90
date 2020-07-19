subroutine BF_mes_thr(mp2,ms2,mf2,T,k)
!Calculating the right hand side of differential equations

  implicit none


  real(8) mp2,ms2,mf2,T,k
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) mu0,p0,p0c
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  real(8) nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,nbSigma,nbd1xSigma,nbd2xSigma,nbd3xSigma,&
          nbd4xSigma,nbd5xSigma


  real(8) b1f2P,b1f2S,b1f3P,b1f3S,b2f1aP,b2f1aS,b2f2aP,b2f2aS,b3f1aP,b3f1aS,b2f3aP,b2f3aS,b3f2aP,b3f2aS

  complex(8) b1f2PC,b1f2SC,b1f3PC,b1f3SC,b2f1aPC,b2f1aSC,b2f2aPC,b2f2aSC,b3f1aPC,b3f1aSC,b2f3aPC,b2f3aSC,b3f2aPC,b3f2aSC

  common /mu0p0_com/ mu0,p0,p0c
  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  common /nbPS_com/ nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,&
                    nbSigma,nbd1xSigma,nbd2xSigma,nbd3xSigma,nbd4xSigma,nbd5xSigma
  common /BF_mes_thr_com/ b1f2P,b1f2S,b1f3P,b1f3S,b2f1aP,b2f1aS,b2f2aP,b2f2aS,b3f1aP,b3f1aS,b2f3aP,b2f3aS,b3f2aP,b3f2aS



  b1f2PC=(k**2*nfa)/(4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                     &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                      &
  (k**3*nfd1xa)/                                                                &
   (4.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                      &
  (k**3*nfd1xf)/                                                                &
   (4.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +                         &
  (k**2*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +      &
  (k**4*nbPion)/                                                                &
   (2.*Sqrt(1 + mp2)*(-(k**2*(1 + mf2)) +                                       &
        (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**2) +                  &
  (k**4*nbPion)/                                                                &
   (2.*Sqrt(1 + mp2)*(-(k**2*(1 + mf2)) +                                       &
        (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**2) -                     &
  (k**3*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                  &
  (k**3*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                     &
  k**2/(4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                  &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)) +                        &
  k**4/(2.*Sqrt(1 + mp2)*(-(k**2*(1 + mf2)) +                                   &
        (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**2) -                 &
  (k**3*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                            &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2)                      

  b1f2P=real(b1f2PC)


  b1f2SC=(k**2*nfa)/(4.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                     &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                      &
  (k**3*nfd1xa)/                                                                &
   (4.*(1 + mf2)*(-(k**2*(1 + ms2)) +                                           &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                      &
  (k**3*nfd1xf)/                                                                &
   (4.*(1 + mf2)*(-(k**2*(1 + ms2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +                         &
  (k**2*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +      &
  (k**4*nbSigma)/                                                               &
   (2.*Sqrt(1 + ms2)*(-(k**2*(1 + mf2)) +                                       &
        (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**2) +                  &
  (k**4*nbSigma)/                                                               &
   (2.*Sqrt(1 + ms2)*(-(k**2*(1 + mf2)) +                                       &
        (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**2) -                     &
  (k**3*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   (2.*(1 + mf2)*(-(k**2*(1 + ms2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                  &
  (k**3*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   (2.*(1 + mf2)*(-(k**2*(1 + ms2)) +                                           &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                     &
  k**2/(4.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                  &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)) +                        &
  k**4/(2.*Sqrt(1 + ms2)*(-(k**2*(1 + mf2)) +                                   &
        (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**2) -                 &
  (k**3*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                            &
   (2.*(1 + mf2)*(-(k**2*(1 + ms2)) +                                           &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2)                      

  b1f2S=real(b1f2SC)


  b1f3PC=(-(k**4*nfa)/(4.*(1 + mf2)**1.5*                                       &
       (-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**2) + (3*k**2*nfa)/                                               &
     (8.*(1 + mf2)**2.5*(-(k**2*(1 + mp2)) +                                    &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (3*k**3*nfd1xa)/                                                            &
     (8.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) +                    &
    (k**4*nfd2xa)/                                                              &
     (8.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (k**4*nff)/(4.*(1 + mf2)**1.5*                                              &
       (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        2) - (3*k**3*nfd1xf)/                                                   &
     (8.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +                       &
    (k**4*nfd2xf)/                                                              &
     (8.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +                       &
    (3*k**2*nff)/                                                               &
     (8.*(1 + mf2)**2.5*(-(k**2*(1 + mp2)) +                                    &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -                       &
    (k**6*nbPion)/                                                              &
     (Sqrt(1 + mp2)*(-(k**2*(1 + mf2)) +                                        &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**6*nbPion)/                                                              &
     (Sqrt(1 + mp2)*(-(k**2*(1 + mf2)) +                                        &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**3*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     (4.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**4*nfd1xa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                 &
     (2.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**4*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/                 &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**4*nfd1xf*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                    &
     (2.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (3*k**3*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                     &
     (4.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (k**4*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/                    &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    k**4/(4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                  &
    (3*k**2)/(8.*(1 + mf2)**2.5*                                                &
       (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2))     &
- k**6/(Sqrt(1 + mp2)*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**3) -               &
    (3*k**3*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     (4.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                  &
    (k**4*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)/                       &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3))/2.                

  b1f3P=real(b1f3PC)


  b1f3SC=(-(k**4*nfa)/(4.*(1 + mf2)**1.5*                                        &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**2) + (3*k**2*nfa)/                                               &
     (8.*(1 + mf2)**2.5*(-(k**2*(1 + ms2)) +                                    &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (3*k**3*nfd1xa)/                                                            &
     (8.*(1 + mf2)**2*(-(k**2*(1 + ms2)) +                                      &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) +                    &
    (k**4*nfd2xa)/                                                              &
     (8.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                    &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (k**4*nff)/(4.*(1 + mf2)**1.5*                                              &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        2) - (3*k**3*nfd1xf)/                                                   &
     (8.*(1 + mf2)**2*(-(k**2*(1 + ms2)) +                                      &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +                       &
    (k**4*nfd2xf)/                                                              &
     (8.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                    &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +                       &
    (3*k**2*nff)/                                                               &
     (8.*(1 + mf2)**2.5*(-(k**2*(1 + ms2)) +                                    &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -                       &
    (k**6*nbSigma)/                                                             &
     (Sqrt(1 + ms2)*(-(k**2*(1 + mf2)) +                                        &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**6*nbSigma)/                                                             &
     (Sqrt(1 + ms2)*(-(k**2*(1 + mf2)) +                                        &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**3*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     (4.*(1 + mf2)**2*(-(k**2*(1 + ms2)) +                                      &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**4*nfd1xa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                 &
     (2.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                    &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**4*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/                 &
     ((1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**4*nfd1xf*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                    &
     (2.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                    &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (3*k**3*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                     &
     (4.*(1 + mf2)**2*(-(k**2*(1 + ms2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (k**4*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/                    &
     ((1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    k**4/(4.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                  &
    (3*k**2)/(8.*(1 + mf2)**2.5*                                                &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2))     &
- k**6/(Sqrt(1 + ms2)*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**3) -               &
    (3*k**3*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     (4.*(1 + mf2)**2*(-(k**2*(1 + ms2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                  &
    (k**4*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)/                       &
     ((1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3))/2.

  b1f3S=real(b1f3SC)


  b2f1aPC=-(k**4*nfa)/(2.*Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) +                    &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) -                  &
  (k**4*nff)/(2.*Sqrt(1 + mf2)*                                                 &
     (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
+ (k**3*nbd1xPion)/                                                             &
   (4.*(1 + mp2)*(-(k**2*(1 + mf2)) +                                           &
       (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)) -                      &
  (k**2*nbPion)/                                                                &
   (4.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                      &
       (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)) +                      &
  (k**3*nbd1xPion)/                                                             &
   (4.*(1 + mp2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)) -                         &
  (k**2*nbPion)/                                                                &
   (4.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                      &
       (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)) +                         &
  (k**3*nbPion*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0))/                   &
   (2.*(1 + mp2)*(-(k**2*(1 + mf2)) +                                           &
        (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**2) -                  &
  (k**3*nbPion*(k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0))/                      &
   (2.*(1 + mp2)*(-(k**2*(1 + mf2)) +                                           &
        (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**2) +                     &
  k**4/(2.*Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) +                                   &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                    &
  k**2/(4.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                  &
       (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)) +                     &
  (k**3*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c))/                         &
   (2.*(1 + mp2)*(-(k**2*(1 + mf2)) +                                           &
        (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**2)

  b2f1aP=real(b2f1aPC)


  b2f1aSC=-(k**4*nfa)/(2.*Sqrt(1 + mf2)*(-(k**2*(1 + ms2)) +                    &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) -                  &
  (k**4*nff)/(2.*Sqrt(1 + mf2)*                                                 &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
+ (k**3*nbd1xSigma)/                                                            &
   (4.*(1 + ms2)*(-(k**2*(1 + mf2)) +                                           &
       (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)) -                      &
  (k**2*nbSigma)/                                                               &
   (4.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                      &
       (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)) +                      &
  (k**3*nbd1xSigma)/                                                            &
   (4.*(1 + ms2)*(-(k**2*(1 + mf2)) +                                           &
       (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)) -                         &
  (k**2*nbSigma)/                                                               &
   (4.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                      &
       (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)) +                         &
  (k**3*nbSigma*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0))/                  &
   (2.*(1 + ms2)*(-(k**2*(1 + mf2)) +                                           &
        (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**2) -                  &
  (k**3*nbSigma*(k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0))/                     &
   (2.*(1 + ms2)*(-(k**2*(1 + mf2)) +                                           &
        (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**2) +                     &
  k**4/(2.*Sqrt(1 + mf2)*(-(k**2*(1 + ms2)) +                                   &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                    &
  k**2/(4.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                  &
       (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)) +                     &
  (k**3*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c))/                         &
   (2.*(1 + ms2)*(-(k**2*(1 + mf2)) +                                           &
        (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**2)

  b2f1aS=real(b2f1aSC)


  b2f2aPC=-(k**4*nfa)/(4.*(1 + mf2)**1.5*                                       &
     (-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) + (k**5*nfd1xa)/                                                       &
   (4.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                  &
  (k**5*nfd1xf)/                                                                &
   (4.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                     &
  (k**4*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
- (k**5*nbd1xPion)/                                                             &
   (4.*(1 + mp2)*(-(k**2*(1 + mf2)) +                                           &
        (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**2) +                  &
  (k**4*nbPion)/                                                                &
   (4.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                      &
        (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**2) -                  &
  (k**5*nbd1xPion)/                                                             &
   (4.*(1 + mp2)*(-(k**2*(1 + mf2)) +                                           &
        (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**2) +                     &
  (k**4*nbPion)/                                                                &
   (4.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                      &
        (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**2) +                     &
  (k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   ((1 + mf2)*(-(k**2*(1 + mp2)) +                                              &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                  &
  (k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   ((1 + mf2)*(-(k**2*(1 + mp2)) +                                              &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                     &
  (k**5*nbPion*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0))/                   &
   ((1 + mp2)*(-(k**2*(1 + mf2)) +                                              &
        (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**3) +                  &
  (k**5*nbPion*(k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0))/                      &
   ((1 + mp2)*(-(k**2*(1 + mf2)) +                                              &
        (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**3) +                     &
  k**4/(4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                  &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) +                    &
  k**4/(4.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                  &
        (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**2) +                 &
  (k**5*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                            &
   ((1 + mf2)*(-(k**2*(1 + mp2)) +                                              &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) -                    &
  (k**5*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c))/                         &
   ((1 + mp2)*(-(k**2*(1 + mf2)) +                                              &
        (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**3)

  b2f2aP=real(b2f2aPC)
  

  b2f2aSC=-(k**4*nfa)/(4.*(1 + mf2)**1.5*                                       &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) + (k**5*nfd1xa)/                                                       &
   (4.*(1 + mf2)*(-(k**2*(1 + ms2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                  &
  (k**5*nfd1xf)/                                                                &
   (4.*(1 + mf2)*(-(k**2*(1 + ms2)) +                                           &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                     &
  (k**4*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
- (k**5*nbd1xSigma)/                                                            &
   (4.*(1 + ms2)*(-(k**2*(1 + mf2)) +                                           &
        (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**2) +                  &
  (k**4*nbSigma)/                                                               &
   (4.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                      &
        (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**2) -                  &
  (k**5*nbd1xSigma)/                                                            &
   (4.*(1 + ms2)*(-(k**2*(1 + mf2)) +                                           &
        (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**2) +                     &
  (k**4*nbSigma)/                                                               &
   (4.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                      &
        (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**2) +                     &
  (k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   ((1 + mf2)*(-(k**2*(1 + ms2)) +                                              &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                  &
  (k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   ((1 + mf2)*(-(k**2*(1 + ms2)) +                                              &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                     &
  (k**5*nbSigma*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0))/                  &
   ((1 + ms2)*(-(k**2*(1 + mf2)) +                                              &
        (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**3) +                  &
  (k**5*nbSigma*(k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0))/                     &
   ((1 + ms2)*(-(k**2*(1 + mf2)) +                                              &
        (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**3) +                     &
  k**4/(4.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                  &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) +                    &
  k**4/(4.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                  &
        (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**2) +                 &
  (k**5*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                            &
   ((1 + mf2)*(-(k**2*(1 + ms2)) +                                              &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) -                    &
  (k**5*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c))/                         &
   ((1 + ms2)*(-(k**2*(1 + mf2)) +                                              &
        (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**3)

  b2f2aS=real(b2f2aSC)


  b3f1aPC=((k**6*nfa)/(Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) +                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) +                &
    (k**6*nff)/(Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) +                              &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (k**4*nbPion)/                                                              &
     (4.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (3*k**3*nbd1xPion)/                                                         &
     (8.*(1 + mp2)**2*(-(k**2*(1 + mf2)) +                                      &
         (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (k**4*nbd2xPion)/                                                           &
     (8.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
         (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (3*k**2*nbPion)/                                                            &
     (8.*(1 + mp2)**2.5*(-(k**2*(1 + mf2)) +                                    &
         (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)) +                    &
    (k**4*nbPion)/                                                              &
     (4.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (3*k**3*nbd1xPion)/                                                         &
     (8.*(1 + mp2)**2*(-(k**2*(1 + mf2)) +                                      &
         (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)) -                       &
    (k**4*nbd2xPion)/                                                           &
     (8.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
         (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)) -                       &
    (3*k**2*nbPion)/                                                            &
     (8.*(1 + mp2)**2.5*(-(k**2*(1 + mf2)) +                                    &
         (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)) -                       &
    (k**4*nbd1xPion*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0))/              &
     (2.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (3*k**3*nbPion*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0))/               &
     (4.*(1 + mp2)**2*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**2) -                &
    (k**4*nbPion*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)/              &
     ((1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**3) +                &
    (k**4*nbd1xPion*(k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0))/                 &
     (2.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (3*k**3*nbPion*(k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0))/                  &
     (4.*(1 + mp2)**2*(-(k**2*(1 + mf2)) +                                      &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (k**4*nbPion*(k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)/                 &
     ((1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    k**6/(Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) +                                    &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    k**4/(4.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**2) -               &
    (3*k**2)/(8.*(1 + mp2)**2.5*                                                &
       (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**     &
          2)) + (3*k**3*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c))/         &
     (4.*(1 + mp2)**2*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**2) -               &
    (k**4*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)/                    &
     ((1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**3))/2.

  b3f1aP=real(b3f1aPC)


  b3f1aSC=((k**6*nfa)/(Sqrt(1 + mf2)*(-(k**2*(1 + ms2)) +                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) +                &
    (k**6*nff)/(Sqrt(1 + mf2)*(-(k**2*(1 + ms2)) +                              &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (k**4*nbSigma)/                                                             &
     (4.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (3*k**3*nbd1xSigma)/                                                        &
     (8.*(1 + ms2)**2*(-(k**2*(1 + mf2)) +                                      &
         (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (k**4*nbd2xSigma)/                                                          &
     (8.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
         (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (3*k**2*nbSigma)/                                                           &
     (8.*(1 + ms2)**2.5*(-(k**2*(1 + mf2)) +                                    &
         (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)) +                    &
    (k**4*nbSigma)/                                                             &
     (4.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (3*k**3*nbd1xSigma)/                                                        &
     (8.*(1 + ms2)**2*(-(k**2*(1 + mf2)) +                                      &
         (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)) -                       &
    (k**4*nbd2xSigma)/                                                          &
     (8.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
         (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)) -                       &
    (3*k**2*nbSigma)/                                                           &
     (8.*(1 + ms2)**2.5*(-(k**2*(1 + mf2)) +                                    &
         (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)) -                       &
    (k**4*nbd1xSigma*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0))/             &
     (2.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (3*k**3*nbSigma*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0))/              &
     (4.*(1 + ms2)**2*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**2) -                &
    (k**4*nbSigma*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)/             &
     ((1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**3) +                &
    (k**4*nbd1xSigma*(k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0))/                &
     (2.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (3*k**3*nbSigma*(k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0))/                 &
     (4.*(1 + ms2)**2*(-(k**2*(1 + mf2)) +                                      &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (k**4*nbSigma*(k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)/                &
     ((1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    k**6/(Sqrt(1 + mf2)*(-(k**2*(1 + ms2)) +                                    &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    k**4/(4.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**2) -               &
    (3*k**2)/(8.*(1 + ms2)**2.5*                                                &
       (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**     &
          2)) + (3*k**3*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c))/         &
     (4.*(1 + ms2)**2*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**2) -               &
    (k**4*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)/                    &
     ((1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**3))/2.

  b3f1aS=real(b3f1aSC)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  b2f3aPC=((k**6*nfa)/(2.*(1 + mf2)**1.5*                                       &
       (-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**3) - (3*k**4*nfa)/                                               &
     (8.*(1 + mf2)**2.5*(-(k**2*(1 + mp2)) +                                    &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (3*k**5*nfd1xa)/                                                            &
     (8.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) -                &
    (k**6*nfd2xa)/                                                              &
     (8.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**6*nff)/(2.*(1 + mf2)**1.5*                                              &
       (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        3) + (3*k**5*nfd1xf)/                                                   &
     (8.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (k**6*nfd2xf)/                                                              &
     (8.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (3*k**4*nff)/                                                               &
     (8.*(1 + mf2)**2.5*(-(k**2*(1 + mp2)) +                                    &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (k**7*nbd1xPion)/                                                           &
     (2.*(1 + mp2)*(-(k**2*(1 + mf2)) +                                         &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**6*nbPion)/                                                              &
     (2.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**3) +                &
    (k**7*nbd1xPion)/                                                           &
     (2.*(1 + mp2)*(-(k**2*(1 + mf2)) +                                         &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (k**6*nbPion)/                                                              &
     (2.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (3*k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     (2.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**6*nfd1xa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                 &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (3*k**6*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/               &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**4) +                &
    (k**6*nfd1xf*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                    &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                     &
     (2.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**6*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/                  &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**4) +                   &
    (3*k**7*nbPion*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0))/               &
     ((1 + mp2)*(-(k**2*(1 + mf2)) +                                            &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**4) -                &
    (3*k**7*nbPion*(k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0))/                  &
     ((1 + mp2)*(-(k**2*(1 + mf2)) +                                            &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**4) -                   &
    k**6/(2.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    (3*k**4)/(8.*(1 + mf2)**2.5*                                                &
       (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**        &
           2)**2) - k**6/                                                       &
     (2.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**3) +               &
    (3*k**5*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     (2.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    (3*k**6*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)/                     &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**4) +                  &
    (3*k**7*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c))/                     &
     ((1 + mp2)*(-(k**2*(1 + mf2)) +                                            &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**4))/2.

  b2f3aP=real(b2f3aPC)


  b2f3aSC=((k**6*nfa)/(2.*(1 + mf2)**1.5*                                       &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**3) - (3*k**4*nfa)/                                               &
     (8.*(1 + mf2)**2.5*(-(k**2*(1 + ms2)) +                                    &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (3*k**5*nfd1xa)/                                                            &
     (8.*(1 + mf2)**2*(-(k**2*(1 + ms2)) +                                      &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) -                &
    (k**6*nfd2xa)/                                                              &
     (8.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                    &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**6*nff)/(2.*(1 + mf2)**1.5*                                              &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        3) + (3*k**5*nfd1xf)/                                                   &
     (8.*(1 + mf2)**2*(-(k**2*(1 + ms2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (k**6*nfd2xf)/                                                              &
     (8.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                    &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (3*k**4*nff)/                                                               &
     (8.*(1 + mf2)**2.5*(-(k**2*(1 + ms2)) +                                    &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (k**7*nbd1xSigma)/                                                          &
     (2.*(1 + ms2)*(-(k**2*(1 + mf2)) +                                         &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**6*nbSigma)/                                                             &
     (2.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**3) +                &
    (k**7*nbd1xSigma)/                                                          &
     (2.*(1 + ms2)*(-(k**2*(1 + mf2)) +                                         &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (k**6*nbSigma)/                                                             &
     (2.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (3*k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     (2.*(1 + mf2)**2*(-(k**2*(1 + ms2)) +                                      &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**6*nfd1xa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                 &
     ((1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (3*k**6*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/               &
     ((1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**4) +                &
    (k**6*nfd1xf*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                    &
     ((1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                     &
     (2.*(1 + mf2)**2*(-(k**2*(1 + ms2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**6*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/                  &
     ((1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**4) +                   &
    (3*k**7*nbSigma*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0))/              &
     ((1 + ms2)*(-(k**2*(1 + mf2)) +                                            &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**4) -                &
    (3*k**7*nbSigma*(k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0))/                 &
     ((1 + ms2)*(-(k**2*(1 + mf2)) +                                            &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**4) -                   &
    k**6/(2.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    (3*k**4)/(8.*(1 + mf2)**2.5*                                                &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**        &
           2)**2) - k**6/                                                       &
     (2.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**3) +               &
    (3*k**5*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     (2.*(1 + mf2)**2*(-(k**2*(1 + ms2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    (3*k**6*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)/                     &
     ((1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**4) +                  &
    (3*k**7*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c))/                     &
     ((1 + ms2)*(-(k**2*(1 + mf2)) +                                            &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**4))/2.

  b2f3aS=real(b2f3aSC)


  b3f2aPC=((k**6*nfa)/(2.*(1 + mf2)**1.5*                                       &
       (-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**3) - (k**7*nfd1xa)/                                              &
     (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                         &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**7*nfd1xf)/                                                              &
     (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                         &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (k**6*nff)/(2.*(1 + mf2)**1.5*                                              &
       (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        3) - (k**6*nbPion)/                                                     &
     (2.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (3*k**5*nbd1xPion)/                                                         &
     (8.*(1 + mp2)**2*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**6*nbd2xPion)/                                                           &
     (8.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (3*k**4*nbPion)/                                                            &
     (8.*(1 + mp2)**2.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**2) -                &
    (k**6*nbPion)/                                                              &
     (2.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**5*nbd1xPion)/                                                         &
     (8.*(1 + mp2)**2*(-(k**2*(1 + mf2)) +                                      &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (k**6*nbd2xPion)/                                                           &
     (8.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (3*k**4*nbPion)/                                                            &
     (8.*(1 + mp2)**2.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (3*k**7*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     ((1 + mf2)*(-(k**2*(1 + mp2)) +                                            &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**4) +                &
    (3*k**7*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                     &
     ((1 + mf2)*(-(k**2*(1 + mp2)) +                                            &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**4) +                   &
    (k**6*nbd1xPion*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0))/              &
     ((1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (3*k**5*nbPion*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0))/               &
     (2.*(1 + mp2)**2*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**3) +                &
    (3*k**6*nbPion*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)/            &
     ((1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0)**2)**4) -                &
    (k**6*nbd1xPion*(k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0))/                 &
     ((1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (3*k**5*nbPion*(k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0))/                  &
     (2.*(1 + mp2)**2*(-(k**2*(1 + mf2)) +                                      &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (3*k**6*nbPion*(k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)/               &
     ((1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (k*Sqrt(1 + mp2) - mu0 + Complex(0,1)*p0)**2)**4) -                   &
    k**6/(2.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) -                  &
    k**6/(2.*(1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**3) +               &
    (3*k**4)/(8.*(1 + mp2)**2.5*                                                &
       (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) - mu0 +                         &
             Complex(0,1)*p0c)**2)**2) -                                        &
    (3*k**7*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     ((1 + mf2)*(-(k**2*(1 + mp2)) +                                            &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**4) -                  &
    (3*k**5*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c))/                     &
     (2.*(1 + mp2)**2*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**3) +               &
    (3*k**6*(-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)/                  &
     ((1 + mp2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (-(k*Sqrt(1 + mp2)) - mu0 + Complex(0,1)*p0c)**2)**4))/2.

  b3f2aP=real(b3f2aPC)


  b3f2aSC=((k**6*nfa)/(2.*(1 + mf2)**1.5*                                       &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**3) - (k**7*nfd1xa)/                                              &
     (2.*(1 + mf2)*(-(k**2*(1 + ms2)) +                                         &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**7*nfd1xf)/                                                              &
     (2.*(1 + mf2)*(-(k**2*(1 + ms2)) +                                         &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (k**6*nff)/(2.*(1 + mf2)**1.5*                                              &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        3) - (k**6*nbSigma)/                                                    &
     (2.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (3*k**5*nbd1xSigma)/                                                        &
     (8.*(1 + ms2)**2*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**6*nbd2xSigma)/                                                          &
     (8.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (3*k**4*nbSigma)/                                                           &
     (8.*(1 + ms2)**2.5*(-(k**2*(1 + mf2)) +                                    &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**2) -                &
    (k**6*nbSigma)/                                                             &
     (2.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**5*nbd1xSigma)/                                                        &
     (8.*(1 + ms2)**2*(-(k**2*(1 + mf2)) +                                      &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (k**6*nbd2xSigma)/                                                          &
     (8.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (3*k**4*nbSigma)/                                                           &
     (8.*(1 + ms2)**2.5*(-(k**2*(1 + mf2)) +                                    &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (3*k**7*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     ((1 + mf2)*(-(k**2*(1 + ms2)) +                                            &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**4) +                &
    (3*k**7*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                     &
     ((1 + mf2)*(-(k**2*(1 + ms2)) +                                            &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**4) +                   &
    (k**6*nbd1xSigma*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0))/             &
     ((1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (3*k**5*nbSigma*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0))/              &
     (2.*(1 + ms2)**2*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**3) +                &
    (3*k**6*nbSigma*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)/           &
     ((1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0)**2)**4) -                &
    (k**6*nbd1xSigma*(k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0))/                &
     ((1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (3*k**5*nbSigma*(k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0))/                 &
     (2.*(1 + ms2)**2*(-(k**2*(1 + mf2)) +                                      &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (3*k**6*nbSigma*(k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)/              &
     ((1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (k*Sqrt(1 + ms2) - mu0 + Complex(0,1)*p0)**2)**4) -                   &
    k**6/(2.*(1 + mf2)**1.5*(-(k**2*(1 + ms2)) +                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) -                  &
    k**6/(2.*(1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**3) +               &
    (3*k**4)/(8.*(1 + ms2)**2.5*                                                &
       (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) - mu0 +                         &
             Complex(0,1)*p0c)**2)**2) -                                        &
    (3*k**7*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     ((1 + mf2)*(-(k**2*(1 + ms2)) +                                            &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**4) -                  &
    (3*k**5*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c))/                     &
     (2.*(1 + ms2)**2*(-(k**2*(1 + mf2)) +                                      &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**3) +               &
    (3*k**6*(-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)/                  &
     ((1 + ms2)**1.5*(-(k**2*(1 + mf2)) +                                       &
          (-(k*Sqrt(1 + ms2)) - mu0 + Complex(0,1)*p0c)**2)**4))/2.

  b3f2aS=real(b3f2aSC)
  

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine F_thr(mf2,k)
!Calculating the right hand side of differential equations

  implicit none


  real(8) mf2,k
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  real(8) f2a,f3a,f4a

  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  common /F_thr_com/ f2a,f3a,f4a


  f2a=(1 - nfa + k*Sqrt(1 + mf2)*nfd1xa + k*Sqrt(1 + mf2)*nfd1xf - nff)/        &
  (4.*(1 + mf2)**1.5)

  f3a=(-(k**2*(1 + mf2)*(nfd2xa + nfd2xf)) +                                    &
    3*(1 - nfa + k*Sqrt(1 + mf2)*nfd1xa + k*Sqrt(1 + mf2)*nfd1xf - nff))/       &
  (16.*(1 + mf2)**2.5)

  f4a=(15 - 15*nfa + 15*k*Sqrt(1 + mf2)*nfd1xa + 15*k*Sqrt(1 + mf2)*nfd1xf -    &
    6*k**2*nfd2xa - 6*k**2*mf2*nfd2xa - 6*k**2*nfd2xf - 6*k**2*mf2*nfd2xf +     &
    k**3*Sqrt(1 + mf2)*nfd3xa + k**3*mf2*Sqrt(1 + mf2)*nfd3xa +                 &
    k**3*Sqrt(1 + mf2)*nfd3xf + k**3*mf2*Sqrt(1 + mf2)*nfd3xf - 15*nff)/        &
  (96.*(1 + mf2)**3.5)

  
end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FT_thr(mf2,k)
!finite temperature part of Fa

  implicit none


  real(8) mf2,k
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  real(8) f2aT,f3aT,f4aT

  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  common /FT_thr_com/ f2aT,f3aT,f4aT


  f2aT=(-nfa + k*Sqrt(1 + mf2)*(nfd1xa + nfd1xf) - nff)/(4.*(1 + mf2)**1.5)

  f3aT=(-3*nfa - k*(-3*Sqrt(1 + mf2)*nfd1xa - 3*Sqrt(1 + mf2)*nfd1xf +           &
       k*(1 + mf2)*(nfd2xa + nfd2xf)) - 3*nff)/(16.*(1 + mf2)**2.5)

  f4aT=(- 15*nfa + 15*k*Sqrt(1 + mf2)*nfd1xa + 15*k*Sqrt(1 + mf2)*nfd1xf -    &
  6*k**2*nfd2xa - 6*k**2*mf2*nfd2xa - 6*k**2*nfd2xf - 6*k**2*mf2*nfd2xf +     &
  k**3*Sqrt(1 + mf2)*nfd3xa + k**3*mf2*Sqrt(1 + mf2)*nfd3xa +                 &
  k**3*Sqrt(1 + mf2)*nfd3xf + k**3*mf2*Sqrt(1 + mf2)*nfd3xf - 15*nff)/        &
  (96.*(1 + mf2)**3.5)


end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FF_thr(mf2,T,mu,l,lb,k,q,costhe)
!Calculating the right hand side of differential equations

  implicit none

  real(8) mf2,T,mu,l,lb,k,q,costhe
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  real(8) Eq
  real(8) Fdnf1,Fdnf2,Fdnf1d1,Fdnf2d1
  real(8) nffqmp,nfaqmp,Eqmp
  real(8) F1F1,F2F1,F1F1phi,F2F1phi


  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  common /nffqmp_com/ nffqmp,nfaqmp,Eqmp
  common /FF_thr_com/ F1F1,F2F1,F1F1phi,F2F1phi



  call nfqmp_com(mf2,T,mu,l,lb,k,q,costhe)


  Eq=k*Sqrt(1 + mf2)



  if(abs(Eq - Eqmp)/Eq<1.d-5)then
    Fdnf1=nfd1xa
    Fdnf2=nfd1xf

    Fdnf1d1=(k**2*nfd2xa)/(4.*Eq)
    Fdnf2d1=(k**2*nfd2xf)/(4.*Eq)
  else
    Fdnf1=(-nfa + nfaqmp)/(-Eq + Eqmp)
    Fdnf2=(nff - nffqmp)/(Eq - Eqmp)

    Fdnf1d1=(k**2*(-nfa + nfaqmp))/(2.*Eq*(-Eq + Eqmp)**2) -                   &
          (k**2*nfd1xa)/(2.*Eq*(-Eq + Eqmp))
    Fdnf2d1=(k**2*nfd1xf)/(2.*Eq*(Eq - Eqmp)) -                                &
          (k**2*(nff - nffqmp))/(2.*Eq*(Eq - Eqmp)**2)
  end if

  F1F1phi=(k**3*(Fdnf1 + Fdnf2 + (1. - nfaqmp - nff)/(Eq + Eqmp) +                 &
      (-1 + nfa + nffqmp)/(-Eq - Eqmp)))/(4.*Eq*Eqmp)

  F1F1=(k**3*(Fdnf1 + Fdnf2 + (0.- nfaqmp - nff)/(Eq + Eqmp) +                 &
      (0. + nfa + nffqmp)/(-Eq - Eqmp)))/(4.*Eq*Eqmp)

  F2F1phi=(k**5*(Fdnf1 + Fdnf2 + (1. - nfaqmp - nff)/(Eq + Eqmp) +                 &
       (-1. + nfa + nffqmp)/(-Eq - Eqmp)))/(8.*Eq**3*Eqmp) -                    &
  (k**3*(Fdnf1d1 + Fdnf2d1 + (k**2*nfd1xa)/(2.*Eq*(-Eq - Eqmp)) -              &
       (k**2*nfd1xf)/(2.*Eq*(Eq + Eqmp)) -                                     &
       (k**2*(1 - nfaqmp - nff))/(2.*Eq*(Eq + Eqmp)**2) +                      &
       (k**2*(-1 + nfa + nffqmp))/(2.*Eq*(-Eq - Eqmp)**2)))/(4.*Eq*Eqmp)

  F2F1=(k**5*(Fdnf1 + Fdnf2 + (0. - nfaqmp - nff)/(Eq + Eqmp) +                 &
       (0. + nfa + nffqmp)/(-Eq - Eqmp)))/(8.*Eq**3*Eqmp) -                    &
  (k**3*(Fdnf1d1 + Fdnf2d1 + (k**2*nfd1xa)/(2.*Eq*(-Eq - Eqmp)) -              &
       (k**2*nfd1xf)/(2.*Eq*(Eq + Eqmp)) -                                     &
       (k**2*(0. - nfaqmp - nff))/(2.*Eq*(Eq + Eqmp)**2) +                      &
       (k**2*(0. + nfa + nffqmp))/(2.*Eq*(-Eq - Eqmp)**2)))/(4.*Eq*Eqmp)


end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FF_vacPT_thr(mf2,T,mu,l,lb,k,q,costhe)
!Calculating the right hand side of differential equations

  implicit none

  real(8) mf2,T,mu,l,lb,k,q,costhe
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  real(8) Eq
  real(8) Fdnf1,Fdnf2,Fdnf1d1,Fdnf2d1
  real(8) nffqmp,nfaqmp,Eqmp
  real(8) F1F1,F2F1


  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  common /nffqmp_com/ nffqmp,nfaqmp,Eqmp
  common /FF_thr_com/ F1F1,F2F1



  call nfqmp_com(mf2,T,mu,l,lb,k,q,costhe)


  Eq=k*Sqrt(1 + mf2)



  if(abs(Eq - Eqmp)/Eq<1.d-5)then
    Fdnf1=nfd1xa
    Fdnf2=nfd1xf

    Fdnf1d1=(k**2*nfd2xa)/(4.*Eq)
    Fdnf2d1=(k**2*nfd2xf)/(4.*Eq)
  else
    Fdnf1=(-nfa + nfaqmp)/(-Eq + Eqmp)
    Fdnf2=(nff - nffqmp)/(Eq - Eqmp)

    Fdnf1d1=(k**2*(-nfa + nfaqmp))/(2.*Eq*(-Eq + Eqmp)**2) -                   &
          (k**2*nfd1xa)/(2.*Eq*(-Eq + Eqmp))
    Fdnf2d1=(k**2*nfd1xf)/(2.*Eq*(Eq - Eqmp)) -                                &
          (k**2*(nff - nffqmp))/(2.*Eq*(Eq - Eqmp)**2)
  end if

  F1F1=(k**3*(Fdnf1 + Fdnf2 + (1 - nfaqmp - nff)/(Eq + Eqmp) +                 &
      (-1 + nfa + nffqmp)/(-Eq - Eqmp)))/(4.*Eq*Eqmp)


  F2F1=(k**5*(Fdnf1 + Fdnf2 + (1 - nfaqmp - nff)/(Eq + Eqmp) +                 &
       (-1 + nfa + nffqmp)/(-Eq - Eqmp)))/(8.*Eq**3*Eqmp) -                    &
  (k**3*(Fdnf1d1 + Fdnf2d1 + (k**2*nfd1xa)/(2.*Eq*(-Eq - Eqmp)) -              &
       (k**2*nfd1xf)/(2.*Eq*(Eq + Eqmp)) -                                     &
       (k**2*(1 - nfaqmp - nff))/(2.*Eq*(Eq + Eqmp)**2) +                      &
       (k**2*(-1 + nfa + nffqmp))/(2.*Eq*(-Eq - Eqmp)**2)))/(4.*Eq*Eqmp)

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BA_thr(T,k,q,p,costhe)
!vacuum terms are not included

  implicit none

  real(8) T,k,q,p,costhe
  real(8) q2t,qmp2t
  real(8) B2B0B1A,B2B1B1A,B3B1A,B2B1B2A,B3B1B1A,B4B1A,B3B2A,B2B0B2A,B3B0B1A
  real(8) F2F0F1Gh

  real(8) nbGluon,nbd1xGluon,nbd2xGluon,nbd3xGluon,nbd4xGluon,nbd5xGluon
  real(8) nfGhost,nfd1xGhost,nfd2xGhost,nfd3xGhost,nfd4xGhost
  real(8) nfqmpGhost
  real(8) nbqGluon,nbd1xqGluon,nbqmpGluon,nbd1xqmpGluon,qmp
  real(8) B3A,B4A,B4cA,B5A,B2B1A,B2A,F3Gh

  common /nbGluon_com/ nbGluon,nbd1xGluon,nbd2xGluon,nbd3xGluon,nbd4xGluon,nbd5xGluon
  common /nfGhost_com/ nfGhost,nfd1xGhost,nfd2xGhost,nfd3xGhost,nfd4xGhost
  common /nbqmp_comm/ nbqGluon,nbd1xqGluon,nbqmpGluon,nbd1xqmpGluon,qmp
  common /nfqmpGhost_comm/ nfqmpGhost
  common /BA_thr_com/ B3A,B4A,B4cA,B5A,B2B1A,B2A,F3Gh



  call nbqmp_com(T,q,p,costhe)
  call nfqmpGhost_com(T,q,p,costhe)

  q2t=(q/k)**2
  qmp2t=(qmp/k)**2


  B2B1A=(-nbGluon + nbqGluon/Sqrt(q2t))/(1 - q2t)**2 -                             &
        (-(k*nbd1xGluon)/2. + nbGluon/2.)/(1 - q2t)

  B2A=-(k*nbd1xGluon)/2. + nbGluon/2.


  if(qmp2t>1.d0)then

    B2B0B1A=(-nbGluon + nbqmpGluon/Sqrt(qmp2t))/(1 - qmp2t)**2 -                   &
            (-(k*nbd1xGluon)/2. + nbGluon/2.)/(1 - qmp2t)
    B3A=B2B0B1A

    B2B1B1A=nbGluon/((1 - q2t)*(1 - qmp2t)**2) + nbGluon/((1 - q2t)**2*(1 - qmp2t)) - &
            (k*nbd1xGluon)/(2.*(1 - q2t)*(1 - qmp2t)) +                            &
            nbGluon/(2.*(1 - q2t)*(1 - qmp2t)) -                                   &
            nbqGluon/((-1 + q2t)**2*Sqrt(q2t)*(q2t - qmp2t)) -                     &
            nbqmpGluon/((-1 + qmp2t)**2*Sqrt(qmp2t)*(-q2t + qmp2t))
    B4A=B2B1B1A

    B2B0B2A=(-2*(-nbGluon + nbqmpGluon/Sqrt(qmp2t)))/(1 - qmp2t)**3 +              &
            (-(k*nbd1xGluon)/2. + nbGluon/2.)/(1 - qmp2t)**2 -                     &
            (-nbqmpGluon/(2.*qmp2t**1.5) + (k*nbd1xqmpGluon)/(2.*qmp2t))/(1 - qmp2t)**2
    B4cA=B2B0B2A

    B2B1B2A=(-2*nbGluon)/((1 - q2t)*(1 - qmp2t)**3) -                              &
            nbGluon/((1 - q2t)**2*(1 - qmp2t)**2) +                                &
            (k*nbd1xGluon)/(2.*(1 - q2t)*(1 - qmp2t)**2) -                         &
            nbGluon/(2.*(1 - q2t)*(1 - qmp2t)**2) +                                &
            nbqGluon/((-1 + q2t)**2*Sqrt(q2t)*(q2t - qmp2t)**2) -                  &
            nbqmpGluon/((-1 + qmp2t)**2*Sqrt(qmp2t)*(-q2t + qmp2t)**2) -           &
            nbqmpGluon/(2.*(-1 + qmp2t)**2*qmp2t**1.5*(-q2t + qmp2t)) +            &
            (k*nbd1xqmpGluon)/(2.*(-1 + qmp2t)**2*qmp2t*(-q2t + qmp2t)) -          &
            (2*nbqmpGluon)/((-1 + qmp2t)**3*Sqrt(qmp2t)*(-q2t + qmp2t))
    B5A=B2B1B2A

    F2F0F1Gh=(nfGhost - nfqmpGhost/Sqrt(qmp2t))/(1 - qmp2t)**2 -                   &
             ((k*nfd1xGhost)/2. - nfGhost/2.)/(1 - qmp2t)
    F3Gh=F2F0F1Gh

  else

    B3A=((-3*k*nbd1xGluon)/4. + (k**2*nbd2xGluon)/4. + (3*nbGluon)/4.)/2.

    B3B1A=((2*(-nbGluon + nbqGluon/Sqrt(q2t)))/(1 - q2t)**3 -                     &
          (2*(-(k*nbd1xGluon)/2. + nbGluon/2.))/(1 - q2t)**2 +                    &
          ((3*k*nbd1xGluon)/4. - (k**2*nbd2xGluon)/4. - (3*nbGluon)/4.)/(1 - q2t))/ &
          2.
    B4A=B3B1A

    B3B0B1A=((2*(-nbGluon + nbqmpGluon/Sqrt(qmp2t)))/(1 - qmp2t)**3 -             &
            (2*(-(k*nbd1xGluon)/2. + nbGluon/2.))/(1 - qmp2t)**2 +                &
            ((3*k*nbd1xGluon)/4. - (k**2*nbd2xGluon)/4. - (3*nbGluon)/4.)/        &
            (1 - qmp2t))/2.
    B4cA=B3B0B1A

    B3B1B1A=((2*nbGluon)/((1 - q2t)*(1 - qmp2t)**3) +                             &
            (2*nbGluon)/((1 - q2t)**2*(1 - qmp2t)**2) -                           &
            (k*nbd1xGluon)/((1 - q2t)*(1 - qmp2t)**2) +                           &
            nbGluon/((1 - q2t)*(1 - qmp2t)**2) +                                  &
            (2*nbGluon)/((1 - q2t)**3*(1 - qmp2t)) -                              &
            (k*nbd1xGluon)/((1 - q2t)**2*(1 - qmp2t)) +                           &
            nbGluon/((1 - q2t)**2*(1 - qmp2t)) -                                  &
            (3*k*nbd1xGluon)/(4.*(1 - q2t)*(1 - qmp2t)) +                         &
            (k**2*nbd2xGluon)/(4.*(1 - q2t)*(1 - qmp2t)) +                        &
            (3*nbGluon)/(4.*(1 - q2t)*(1 - qmp2t)) +                              &
            (2*nbqGluon)/((-1 + q2t)**3*Sqrt(q2t)*(q2t - qmp2t)) +                &
            (2*nbqmpGluon)/((-1 + qmp2t)**3*Sqrt(qmp2t)*(-q2t + qmp2t)))/2.
    B5A=B3B1B1A

    F3Gh=((3*k*nfd1xGhost)/4. - (k**2*nfd2xGhost)/4. - (3*nfGhost)/4.)/2.
    
  end if

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine F_thr_corr(mf2,k,T)
!Calculating the right hand side of differential equations

  implicit none


  real(8) mf2,k
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  real(8) f2a,f3a,f4a
  real(8) T
  real(8) pi
  parameter(pi=3.1415926)

  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  common /F_thr_com/ f2a,f3a,f4a


  f2a=((1 + mf2)**(-1.5) - nfa/(1 + mf2)**1.5 + (k*nfd1xa)/(1 + mf2) +          &
    (k*nfd1xf)/(1 + mf2) - nff/(1 + mf2)**1.5 -                                 &
    (4*T)/(k*(1 + mf2 + (Pi**2*T**2)/k**2)**2))/4.

  f3a=(3/(1 + mf2)**2.5 - (3*nfa)/(1 + mf2)**2.5 + (3*k*nfd1xa)/(1 + mf2)**2 +  &
    (3*k*nfd1xf)/(1 + mf2)**2 - (k**2*nfd2xa)/(1 + mf2)**1.5 -                  &
    (k**2*nfd2xf)/(1 + mf2)**1.5 - (3*nff)/(1 + mf2)**2.5 -                     &
    (16*T)/(k*(1 + mf2 + (Pi**2*T**2)/k**2)**3))/16.

  f4a=(15/(1 + mf2)**3.5 - (15*nfa)/(1 + mf2)**3.5 + (15*k*nfd1xa)/(1 + mf2)**3 + &
    (15*k*nfd1xf)/(1 + mf2)**3 - (6*k**2*nfd2xa)/(1 + mf2)**2.5 -               &
    (6*k**2*nfd2xf)/(1 + mf2)**2.5 + (k**3*nfd3xa)/(1 + mf2)**2 +               &
    (k**3*nfd3xf)/(1 + mf2)**2 - (15*nff)/(1 + mf2)**3.5 -                      &
    (96*T)/(k*(1 + mf2 + (Pi**2*T**2)/k**2)**4))/96.


  
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BB_thr(mb2a,mb2b,k)
!Calculating the right hand side of differential equations

  implicit none

  real(8) mb2a,mb2b,k


  real(8) nbBa,nbBb,nbBad1,nbBbd1
  real(8) finvEBa,finvEBb,finvEBad1,finvEBbd1

  real(8) fabmbBaBbd1ad0b,fabmbBaBbd0ad1b,fabmbBaBbd1ad1b

  real(8) b2b2

  common /nbBad_com/ nbBa,nbBb,nbBad1,nbBbd1
  common /finvEBad_com/ finvEBa,finvEBb,finvEBad1,finvEBbd1
  common /BB_thr_com/ b2b2

  fabmbBaBbd1ad0b=-(1/(k**2*(mb2a - mb2b)**2))
  fabmbBaBbd0ad1b=1/(k**2*(mb2a - mb2b)**2)
  fabmbBaBbd1ad1b=-2/(k**2*(mb2a - mb2b)**3)


  b2b2=-(k**3*(fabmbBaBbd0ad1b*(finvEBad1 + 2*finvEBad1*nbBa + 2*finvEBa*nbBad1) + &
       fabmbBaBbd1ad1b*(finvEBa + 2*finvEBa*nbBa - finvEBb*(1 + 2*nbBb)) -      &
       fabmbBaBbd1ad0b*(finvEBbd1 + 2*finvEBbd1*nbBb + 2*finvEBb*nbBbd1)))/2.

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BBF_thr(mp2,ms2,mf2,T,k)
!Calculating the right hand side of differential equations

  implicit none


  real(8) mp2,ms2,mf2,T,k
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) mu0,p0,p0c
  real(8) nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,&
          nbSigma,nbd1xSigma,nbd2xSigma,nbd3xSigma,nbd4xSigma,nbd5xSigma
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  
  complex(8) B2B1F1aPSC,B2B1F1aSPC,B1B1F2aPSC,B1B1F3aPSC,B2B1F2aPSC,B2B1F2aSPC,B2F3aPC,B2F3aSC,B3F2aPC,B3F2aSC
  real(8) B2B1F1aPS,B2B1F1aSP,B1B1F2aPS,B1B1F3aPS,B2B1F2aPS,B2B1F2aSP

  common /mu0p0_com/ mu0,p0,p0c
  common /nbPS_com/ nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,&
                    nbSigma,nbd1xSigma,nbd2xSigma,nbd3xSigma,nbd4xSigma,nbd5xSigma
  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  common /BBF_thr_com/ B2B1F1aPS,B2B1F1aSP,B1B1F2aPS,B1B1F3aPS,B2B1F2aPS,B2B1F2aSP
 


  B2B1F1aPSC=-(k**3*nbd1xPion)/                                                 &
   (4.*(1 + mp2)*(mp2 - ms2)*(-(k**2*(1 + mf2)) +                               &
       (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2)) +                      &
  (k**2*nbPion)/                                                                &
   (2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                                            &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2))     &
+ (k**2*nbPion)/                                                                &
   (4.*(1 + mp2)**1.5*(mp2 - ms2)*                                              &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2))     &
- (k**3*nbd1xPion)/                                                             &
   (4.*(1 + mp2)*(mp2 - ms2)*(-(k**2*(1 + mf2)) +                               &
       (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)) +                         &
  (k**2*nbPion)/                                                                &
   (2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                                            &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)) +      &
  (k**2*nbPion)/                                                                &
   (4.*(1 + mp2)**1.5*(mp2 - ms2)*                                              &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)) -      &
  (k**2*nbSigma)/                                                               &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                           &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2))     &
- (k**2*nbSigma)/                                                               &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                           &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)) +      &
  (k**6*nfa)/(2.*Sqrt(1 + mf2)*                                                 &
     (-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**        &
         2)**2*(-(k**2*(1 + ms2)) +                                             &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) +                      &
  (k**6*nff)/(2.*Sqrt(1 + mf2)*                                                 &
     (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2*     &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -      &
  (k**3*nbPion*(-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0))/                   &
   (2.*(1 + mp2)*(mp2 - ms2)*(-(k**2*(1 + mf2)) +                               &
        (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2)**2) +                  &
  (k**3*nbPion*(k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0))/                      &
   (2.*(1 + mp2)*(mp2 - ms2)*(-(k**2*(1 + mf2)) +                               &
        (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)**2) +                     &
  k**2/(2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                                        &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c)**2))    &
+ k**2/(4.*(1 + mp2)**1.5*(mp2 - ms2)*                                          &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c)**2))    &
- k**2/(2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                       &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c)**2))    &
- k**6/(2.*Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) +                                   &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**2*                   &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2))    &
- (k**3*(-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c))/                         &
   (2.*(1 + mp2)*(mp2 - ms2)*(-(k**2*(1 + mf2)) +                               &
        (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c)**2)**2)                   

  B2B1F1aPS=real(B2B1F1aPSC)


  B2B1F1aSPC=-(k**2*nbPion)/(2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                   &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2))     &
- (k**2*nbPion)/                                                                &
   (2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                                            &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)) -      &
  (k**3*nbd1xSigma)/                                                            &
   (4.*(1 + ms2)*(-mp2 + ms2)*(-(k**2*(1 + mf2)) +                              &
       (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2)) +                      &
  (k**2*nbSigma)/                                                               &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                           &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2))     &
+ (k**2*nbSigma)/                                                               &
   (4.*(1 + ms2)**1.5*(-mp2 + ms2)*                                             &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2))     &
- (k**3*nbd1xSigma)/                                                            &
   (4.*(1 + ms2)*(-mp2 + ms2)*(-(k**2*(1 + mf2)) +                              &
       (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)) +                         &
  (k**2*nbSigma)/                                                               &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                           &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)) +      &
  (k**2*nbSigma)/                                                               &
   (4.*(1 + ms2)**1.5*(-mp2 + ms2)*                                             &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)) +      &
  (k**6*nfa)/(2.*Sqrt(1 + mf2)*                                                 &
     (-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*     &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) + (k**6*nff)/                                                          &
   (2.*Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) +                                       &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                           &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
- (k**3*nbSigma*(-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0))/                  &
   (2.*(1 + ms2)*(-mp2 + ms2)*(-(k**2*(1 + mf2)) +                              &
        (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2)**2) +                  &
  (k**3*nbSigma*(k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0))/                     &
   (2.*(1 + ms2)*(-mp2 + ms2)*(-(k**2*(1 + mf2)) +                              &
        (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)**2) -                     &
  k**2/(2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                                        &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c)**2))    &
+ k**2/(2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                       &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c)**2))    &
+ k**2/(4.*(1 + ms2)**1.5*(-mp2 + ms2)*                                         &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c)**2))    &
- k**6/(2.*Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) +                                   &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)*                       &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**       &
         2)**2) - (k**3*(-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c))/         &
   (2.*(1 + ms2)*(-mp2 + ms2)*(-(k**2*(1 + mf2)) +                              &
        (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c)**2)**2)                   

  B2B1F1aSP=real(B2B1F1aSPC)

  B1B1F2aPSC=-(k**4*nbPion)/(2.*Sqrt(1 + mp2)*(mp2 - ms2)*                      &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2)**    &
      2) - (k**4*nbPion)/                                                       &
   (2.*Sqrt(1 + mp2)*(mp2 - ms2)*                                               &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)**2)     &
- (k**4*nbSigma)/                                                               &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)*                                              &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2)**    &
      2) - (k**4*nbSigma)/                                                      &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)*                                              &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)**2)     &
- (k**4*nfa)/(4.*(1 + mf2)**1.5*                                                &
     (-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*     &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2))     &
+ (k**5*nfd1xa)/                                                                &
   (4.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                        &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2))     &
+ (k**5*nfd1xf)/                                                                &
   (4.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                           &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -      &
  (k**4*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*        &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +      &
  (k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                        &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) + (k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/             &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2*                    &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2))     &
- (k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                           &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
- (k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2*                       &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -      &
  k**4/(2.*Sqrt(1 + mp2)*(mp2 - ms2)*                                           &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c)**       &
         2)**2) - k**4/                                                         &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)*                                              &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c)**       &
         2)**2) + k**4/                                                         &
   (4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                      &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)*                       &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2))    &
- (k**5*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c))/                         &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)*                       &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**       &
         2)**2) - (k**5*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c))/         &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**2*                   &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2))    

  B1B1F2aPS=real(B1B1F2aPSC)
  

  B1B1F3aPSC=((k**6*nbPion)/(Sqrt(1 + mp2)*(mp2 - ms2)*                         &
       (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**      &
           2)**3) + (k**6*nbPion)/                                              &
     (Sqrt(1 + mp2)*(mp2 - ms2)*                                                &
       (-(k**2*(1 + mf2)) + (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)**     &
        3) + (k**6*nbSigma)/                                                    &
     (Sqrt(1 + ms2)*(-mp2 + ms2)*                                               &
       (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**      &
           2)**3) + (k**6*nbSigma)/                                             &
     (Sqrt(1 + ms2)*(-mp2 + ms2)*                                               &
       (-(k**2*(1 + mf2)) + (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)**     &
        3) + (k**6*nfa)/                                                        &
     (4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                      &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**2) + (k**6*nfa)/                                                 &
     (4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2*                  &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
          2)) - (3*k**4*nfa)/                                                   &
     (8.*(1 + mf2)**2.5*(-(k**2*(1 + mp2)) +                                    &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                      &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
          2)) + (3*k**5*nfd1xa)/                                                &
     (8.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                      &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
          2)) - (k**6*nfd2xa)/                                                  &
     (8.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                      &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
          2)) + (k**6*nff)/                                                     &
     (4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                         &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        2) + (k**6*nff)/                                                        &
     (4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2*                     &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2))      &
+ (3*k**5*nfd1xf)/                                                              &
     (8.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                         &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2))      &
- (k**6*nfd2xf)/                                                                &
     (8.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                         &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2))      &
- (3*k**4*nff)/(8.*(1 + mf2)**2.5*                                              &
       (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*      &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2))      &
+ (3*k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                    &
     (4.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                      &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**2) - (k**6*nfd1xa*                                               &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                           &
     (2.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                      &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**2) + (3*k**5*nfa*                                                &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                           &
     (4.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2*                  &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
          2)) - (k**6*nfd1xa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/     &
     (2.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2*                  &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
          2)) - (k**6*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/     &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                      &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**3) - (k**6*nfa*                                                  &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/                        &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2*                  &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**2) - (k**6*nfa*                                                  &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/                        &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3*                  &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
          2)) + (k**6*nfd1xf*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/        &
     (2.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                         &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        2) - (3*k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/            &
     (4.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                         &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        2) + (k**6*nfd1xf*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/           &
     (2.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                    &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2*                     &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2))      &
- (3*k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                       &
     (4.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2*                     &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2))      &
- (k**6*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/                      &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                         &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        3) - (k**6*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/           &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2*                     &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        2) - (k**6*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/           &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3*                     &
       (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2))      &
+ k**6/(Sqrt(1 + mp2)*(mp2 - ms2)*                                              &
       (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 -                         &
             Complex(0,1)*p0c)**2)**3) +                                        &
    k**6/(Sqrt(1 + ms2)*(-mp2 + ms2)*                                           &
       (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 -                         &
             Complex(0,1)*p0c)**2)**3) -                                        &
    k**6/(4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)*                     &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 +                         &
             Complex(0,1)*p0c)**2)**2) -                                        &
    k**6/(4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**2*                 &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**     &
          2)) + (3*k**4)/                                                       &
     (8.*(1 + mf2)**2.5*(-(k**2*(1 + mp2)) +                                    &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)*                     &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**     &
          2)) - (3*k**5*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c))/         &
     (4.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)*                     &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 +                         &
             Complex(0,1)*p0c)**2)**2) -                                        &
    (3*k**5*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c))/                     &
     (4.*(1 + mf2)**2*(-(k**2*(1 + mp2)) +                                      &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**2*                 &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**     &
          2)) + (k**6*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)/        &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)*                     &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 +                         &
             Complex(0,1)*p0c)**2)**3) +                                        &
    (k**6*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)/                    &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**2*                 &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 +                         &
             Complex(0,1)*p0c)**2)**2) +                                        &
    (k**6*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)/                    &
     ((1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                       &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**3*                 &
       (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**     &
          2)))/2.

  B1B1F3aPS=real(B1B1F3aPSC)
  

  B2B1F2aPSC=(k**5*nbd1xPion)/                                                  &
   (4.*(1 + mp2)*(mp2 - ms2)*(-(k**2*(1 + mf2)) +                               &
        (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2)**2) -                  &
  (k**4*nbPion)/                                                                &
   (2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                                            &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2)**    &
      2) - (k**4*nbPion)/                                                       &
   (4.*(1 + mp2)**1.5*(mp2 - ms2)*                                              &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2)**    &
      2) + (k**5*nbd1xPion)/                                                    &
   (4.*(1 + mp2)*(mp2 - ms2)*(-(k**2*(1 + mf2)) +                               &
        (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)**2) -                     &
  (k**4*nbPion)/                                                                &
   (2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                                            &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)**2)     &
- (k**4*nbPion)/                                                                &
   (4.*(1 + mp2)**1.5*(mp2 - ms2)*                                              &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)**2)     &
+ (k**4*nbSigma)/                                                               &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                           &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2)**    &
      2) + (k**4*nbSigma)/                                                      &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                           &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)**2)     &
+ (k**6*nfa)/(4.*(1 + mf2)**1.5*                                                &
     (-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**        &
         2)**2*(-(k**2*(1 + ms2)) +                                             &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                      &
  (k**7*nfd1xa)/                                                                &
   (4.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2*                    &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2))     &
- (k**7*nfd1xf)/                                                                &
   (4.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2*                       &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +      &
  (k**6*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2*     &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +      &
  (k**5*nbPion*(-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0))/                   &
   ((1 + mp2)*(mp2 - ms2)*(-(k**2*(1 + mf2)) +                                  &
        (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2)**3) -                  &
  (k**5*nbPion*(k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0))/                      &
   ((1 + mp2)*(mp2 - ms2)*(-(k**2*(1 + mf2)) +                                  &
        (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)**3) -                     &
  (k**7*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2*                    &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) - (k**7*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/             &
   ((1 + mf2)*(-(k**2*(1 + mp2)) +                                              &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3*                    &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2))     &
+ (k**7*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2*                       &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
+ (k**7*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   ((1 + mf2)*(-(k**2*(1 + mp2)) +                                              &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3*                       &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -      &
  k**4/(2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                                        &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c)**       &
         2)**2) - k**4/                                                         &
   (4.*(1 + mp2)**1.5*(mp2 - ms2)*                                              &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c)**       &
         2)**2) + k**4/                                                         &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                           &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c)**       &
         2)**2) - k**6/                                                         &
   (4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                      &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**2*                   &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2))    &
+ (k**5*(-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c))/                         &
   ((1 + mp2)*(mp2 - ms2)*(-(k**2*(1 + mf2)) +                                  &
        (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c)**2)**3) +                 &
  (k**7*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c))/                         &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**2*                   &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**       &
         2)**2) + (k**7*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c))/         &
   ((1 + mf2)*(-(k**2*(1 + mp2)) +                                              &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**3*                   &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2))

  B2B1F2aPS=real(B2B1F2aPSC)


  B2B1F2aSPC=(k**4*nbPion)/(2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                    &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0)**2)**    &
      2) + (k**4*nbPion)/                                                       &
   (2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                                            &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + mp2) + mu0 - Complex(0,1)*p0)**2)**2)     &
+ (k**5*nbd1xSigma)/                                                            &
   (4.*(1 + ms2)*(-mp2 + ms2)*(-(k**2*(1 + mf2)) +                              &
        (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2)**2) -                  &
  (k**4*nbSigma)/                                                               &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                           &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2)**    &
      2) - (k**4*nbSigma)/                                                      &
   (4.*(1 + ms2)**1.5*(-mp2 + ms2)*                                             &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2)**    &
      2) + (k**5*nbd1xSigma)/                                                   &
   (4.*(1 + ms2)*(-mp2 + ms2)*(-(k**2*(1 + mf2)) +                              &
        (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)**2) -                     &
  (k**4*nbSigma)/                                                               &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                           &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)**2)     &
- (k**4*nbSigma)/                                                               &
   (4.*(1 + ms2)**1.5*(-mp2 + ms2)*                                             &
     (-(k**2*(1 + mf2)) + (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)**2)     &
+ (k**6*nfa)/(4.*(1 + mf2)**1.5*                                                &
     (-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*     &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) - (k**7*nfd1xa)/                                                       &
   (4.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                        &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) - (k**7*nfd1xf)/                                                       &
   (4.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                           &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
+ (k**6*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*        &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
+ (k**5*nbSigma*(-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0))/                  &
   ((1 + ms2)*(-mp2 + ms2)*(-(k**2*(1 + mf2)) +                                 &
        (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0)**2)**3) -                  &
  (k**5*nbSigma*(k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0))/                     &
   ((1 + ms2)*(-mp2 + ms2)*(-(k**2*(1 + mf2)) +                                 &
        (k*Sqrt(1 + ms2) + mu0 - Complex(0,1)*p0)**2)**3) -                     &
  (k**7*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   ((1 + mf2)*(-(k**2*(1 + mp2)) +                                              &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)*                        &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      3) - (k**7*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/             &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2*                    &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) + (k**7*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                &
   ((1 + mf2)*(-(k**2*(1 + mp2)) +                                              &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)*                           &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3)     &
+ (k**7*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2*                       &
     (-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
+ k**4/(2.*Sqrt(1 + mp2)*(mp2 - ms2)**2*                                        &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + mp2)) + mu0 - Complex(0,1)*p0c)**       &
         2)**2) - k**4/                                                         &
   (2.*Sqrt(1 + ms2)*(-mp2 + ms2)**2*                                           &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c)**       &
         2)**2) - k**4/                                                         &
   (4.*(1 + ms2)**1.5*(-mp2 + ms2)*                                             &
     (-(k**2*(1 + mf2)) + (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c)**       &
         2)**2) - k**6/                                                         &
   (4.*(1 + mf2)**1.5*(-(k**2*(1 + mp2)) +                                      &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)*                       &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**       &
         2)**2) + (k**5*(-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c))/         &
   ((1 + ms2)*(-mp2 + ms2)*(-(k**2*(1 + mf2)) +                                 &
        (-(k*Sqrt(1 + ms2)) + mu0 - Complex(0,1)*p0c)**2)**3) +                 &
  (k**7*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c))/                         &
   ((1 + mf2)*(-(k**2*(1 + mp2)) +                                              &
       (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)*                       &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**       &
         2)**3) + (k**7*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c))/         &
   (2.*(1 + mf2)*(-(k**2*(1 + mp2)) +                                           &
        (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**2*                   &
     (-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0c)**2)**   &
      2)

  B2B1F2aSP=real(B2B1F2aSPC)


end


subroutine BF_A_thr(mf2,T,k)
!Calculating the right hand side of differential equations

  implicit none


  real(8) mf2,T,k
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) mu0,p0,p0c
  real(8) nbGluon,nbd1xGluon,nbd2xGluon,nbd3xGluon,nbd4xGluon,nbd5xGluon
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa

  real(8) b1f1A,b2f1aA,b1f2A,b2f2aA,b1f3A,b3f1aA,b3f2aA,b2f3aA
  complex(8) b1f1AC,b2f1aAC,b1f2AC,b2f2aAC,b1f3AC,b3f1aAC,b3f2aAC,b2f3aAC

  common /mu0p0_com/ mu0,p0,p0c
  common /nbGluon_com/ nbGluon,nbd1xGluon,nbd2xGluon,nbd3xGluon,nbd4xGluon,nbd5xGluon
  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  common /BF_A_thr_com/ b1f1A,b2f1aA,b1f2A,b2f2aA,b1f3A,b3f1aA,b3f2aA,b2f3aA


  b1f1AC=-(k**2*nbGluon)/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) - &
  (k**2*nbGluon)/(2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) +    &
  (k**2*nfa)/(2.*Sqrt(1 + mf2)*                                                 &
     (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) +               &
  (k**2*nff)/(2.*Sqrt(1 + mf2)*                                                 &
     (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -                  &
  k**2/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)) -            &
  k**2/(2.*Sqrt(1 + mf2)*(-k**2 +                                               &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2))                          

  b1f1A=real(b1f1AC)


  b2f1aAC=(k**3*nbd1xGluon)/                                                    &
   (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) -                 &
  (k**2*nbGluon)/                                                               &
   (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) +                 &
  (k**3*nbd1xGluon)/                                                            &
   (4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) -                  &
  (k**2*nbGluon)/(4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) -    &
  (k**4*nfa)/(2.*Sqrt(1 + mf2)*                                                 &
     (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) -            &
  (k**4*nff)/(2.*Sqrt(1 + mf2)*                                                 &
     (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +               &
  (k**3*nbGluon*(-k - mu0 + Complex(0,1)*p0))/                                  &
   (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) -              &
  (k**3*nbGluon*(k - mu0 + Complex(0,1)*p0))/                                   &
   (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) -               &
  k**2/(4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)) +            &
  k**4/(2.*Sqrt(1 + mf2)*(-k**2 +                                               &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) +                    &
  (k**3*(-k - mu0 + Complex(0,1)*p0c))/                                         &
   (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**2)               

  b2f1aA=real(b2f1aAC)
  

  b1f2AC=(k**4*nbGluon)/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)** &
      2) + (k**4*nbGluon)/                                                      &
   (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +               &
  (k**2*nfa)/(4.*(1 + mf2)**1.5*                                                &
     (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -               &
  (k**3*nfd1xa)/                                                                &
   (4.*(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2))     &
- (k**3*nfd1xf)/                                                                &
   (4.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +      &
  (k**2*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -                  &
  (k**3*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   (2.*(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) + (k**3*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                &
   (2.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
+ k**4/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**2) -         &
  k**2/(4.*(1 + mf2)**1.5*(-k**2 +                                              &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)) -                        &
  (k**3*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                            &
   (2.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2)    

  b1f2A=real(b1f2AC)


  b2f2aAC=-(k**5*nbd1xGluon)/                                                   &
   (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) +              &
  (k**4*nbGluon)/                                                               &
   (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) -              &
  (k**5*nbd1xGluon)/                                                            &
   (4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +               &
  (k**4*nbGluon)/                                                               &
   (4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) -               &
  (k**4*nfa)/(4.*(1 + mf2)**1.5*                                                &
     (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +            &
  (k**5*nfd1xa)/                                                                &
   (4.*(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) + (k**5*nfd1xf)/                                                       &
   (4.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
- (k**4*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -               &
  (k**5*nbGluon*(-k - mu0 + Complex(0,1)*p0))/                                  &
   (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3 +                   &
  (k**5*nbGluon*(k - mu0 + Complex(0,1)*p0))/                                   &
   (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3 +                    &
  (k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   ((1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3)     &
- (k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   ((1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +      &
  k**4/(4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**2) +         &
  k**4/(4.*(1 + mf2)**1.5*(-k**2 +                                              &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                    &
  (k**5*(-k - mu0 + Complex(0,1)*p0c))/                                         &
   (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3 +                  &
  (k**5*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                            &
   ((1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3)       

  b2f2aA=real(b2f2aAC)


  b1f3AC=(-((k**6*nbGluon)/                                                     &
       (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3) -              &
    (k**6*nbGluon)/                                                             &
     (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3 -                  &
    (k**4*nfa)/(4.*(1 + mf2)**1.5*                                              &
       (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +          &
    (3*k**2*nfa)/                                                               &
     (8.*(1 + mf2)**2.5*(-k**2 +                                                &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (3*k**3*nfd1xa)/                                                            &
     (8.*(1 + mf2)**2*(-k**2 +                                                  &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) +                    &
    (k**4*nfd2xa)/                                                              &
     (8.*(1 + mf2)**1.5*(-k**2 +                                                &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (k**4*nff)/(4.*(1 + mf2)**1.5*                                              &
       (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -             &
    (3*k**3*nfd1xf)/                                                            &
     (8.*(1 + mf2)**2*(-k**2 +                                                  &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +                       &
    (k**4*nfd2xf)/                                                              &
     (8.*(1 + mf2)**1.5*(-k**2 +                                                &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +                       &
    (3*k**2*nff)/                                                               &
     (8.*(1 + mf2)**2.5*(-k**2 +                                                &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -                       &
    (3*k**3*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     (4.*(1 + mf2)**2*(-k**2 +                                                  &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**4*nfd1xa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                 &
     (2.*(1 + mf2)**1.5*(-k**2 +                                                &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**4*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/                 &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**4*nfd1xf*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                    &
     (2.*(1 + mf2)**1.5*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (3*k**3*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                     &
     (4.*(1 + mf2)**2*(-k**2 +                                                  &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (k**4*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/                    &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    k**6/(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3 +            &
    k**4/(4.*(1 + mf2)**1.5*(-k**2 +                                            &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                  &
    (3*k**2)/(8.*(1 + mf2)**2.5*                                                &
       (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)) -               &
    (3*k**3*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     (4.*(1 + mf2)**2*(-k**2 +                                                  &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                  &
    (k**4*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)/                       &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3))/2.                

  b1f3A=real(b1f3AC)


  b3f1aAC=((k**4*nbGluon)/                                                      &
     (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) +            &
    (3*k**3*nbd1xGluon)/                                                        &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) -               &
    (k**4*nbd2xGluon)/                                                          &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) -               &
    (3*k**2*nbGluon)/                                                           &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) +               &
    (k**4*nbGluon)/                                                             &
     (4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +             &
    (3*k**3*nbd1xGluon)/                                                        &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) -                &
    (k**4*nbd2xGluon)/                                                          &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) -                &
    (3*k**2*nbGluon)/                                                           &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) +                &
    (k**6*nfa)/(Sqrt(1 + mf2)*(-k**2 +                                          &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) +                &
    (k**6*nff)/(Sqrt(1 + mf2)*(-k**2 +                                          &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (k**4*nbd1xGluon*(-k - mu0 + Complex(0,1)*p0))/                             &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) +            &
    (3*k**3*nbGluon*(-k - mu0 + Complex(0,1)*p0))/                              &
     (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) -            &
    (k**4*nbGluon*(-k - mu0 + Complex(0,1)*p0)**2)/                             &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3 +                 &
    (k**4*nbd1xGluon*(k - mu0 + Complex(0,1)*p0))/                              &
     (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) -             &
    (3*k**3*nbGluon*(k - mu0 + Complex(0,1)*p0))/                               &
     (4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) -             &
    (k**4*nbGluon*(k - mu0 + Complex(0,1)*p0)**2)/                              &
     (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3 +                  &
    k**4/(4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**2) -       &
    (3*k**2)/(8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)) -      &
    k**6/(Sqrt(1 + mf2)*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    (3*k**3*(-k - mu0 + Complex(0,1)*p0c))/                                     &
     (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**2) -           &
    (k**4*(-k - mu0 + Complex(0,1)*p0c)**2)/                                    &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3)/2.              

  b3f1aA=real(b3f1aAC)


  b3f2aAC=(-(k**6*nbGluon)/                                                     &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3) -            &
    (3*k**5*nbd1xGluon)/                                                        &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) +            &
    (k**6*nbd2xGluon)/                                                          &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) +            &
    (3*k**4*nbGluon)/                                                           &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) -            &
    (k**6*nbGluon)/                                                             &
     (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3) -             &
    (3*k**5*nbd1xGluon)/                                                        &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +             &
    (k**6*nbd2xGluon)/                                                          &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +             &
    (3*k**4*nbGluon)/                                                           &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +             &
    (k**6*nfa)/(2.*(1 + mf2)**1.5*                                              &
       (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -          &
    (k**7*nfd1xa)/                                                              &
     (2.*(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**3) - (k**7*nfd1xf)/                                              &
     (2.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        3) + (k**6*nff)/                                                        &
     (2.*(1 + mf2)**1.5*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (k**6*nbd1xGluon*(-k - mu0 + Complex(0,1)*p0))/                             &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3 -                 &
    (3*k**5*nbGluon*(-k - mu0 + Complex(0,1)*p0))/                              &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3) +            &
    (3*k**6*nbGluon*(-k - mu0 + Complex(0,1)*p0)**2)/                           &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**4 -                 &
    (k**6*nbd1xGluon*(k - mu0 + Complex(0,1)*p0))/                              &
     (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3 +                  &
    (3*k**5*nbGluon*(k - mu0 + Complex(0,1)*p0))/                               &
     (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3) +             &
    (3*k**6*nbGluon*(k - mu0 + Complex(0,1)*p0)**2)/                            &
     (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**4 -                  &
    (3*k**7*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     ((1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**     &
        4) + (3*k**7*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/            &
     ((1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**4)      &
- k**6/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3) +         &
    (3*k**4)/(8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**       &
        2) - k**6/                                                              &
     (2.*(1 + mf2)**1.5*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) -                  &
    (3*k**5*(-k - mu0 + Complex(0,1)*p0c))/                                     &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3) +           &
    (3*k**6*(-k - mu0 + Complex(0,1)*p0c)**2)/                                  &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**4 -                &
    (3*k**7*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     ((1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**4))/   &
  2.                                                                            

  b3f2aA=real(b3f2aAC)


  b2f3aAC=((k**7*nbd1xGluon)/                                                   &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3) -            &
    (k**6*nbGluon)/                                                             &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3) +            &
    (k**7*nbd1xGluon)/                                                          &
     (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3) -             &
    (k**6*nbGluon)/                                                             &
     (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3) +             &
    (k**6*nfa)/(2.*(1 + mf2)**1.5*                                              &
       (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -          &
    (3*k**4*nfa)/                                                               &
     (8.*(1 + mf2)**2.5*(-k**2 +                                                &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (3*k**5*nfd1xa)/                                                            &
     (8.*(1 + mf2)**2*(-k**2 +                                                  &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) -                &
    (k**6*nfd2xa)/                                                              &
     (8.*(1 + mf2)**1.5*(-k**2 +                                                &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**6*nff)/(2.*(1 + mf2)**1.5*                                              &
       (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +             &
    (3*k**5*nfd1xf)/                                                            &
     (8.*(1 + mf2)**2*(-k**2 +                                                  &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (k**6*nfd2xf)/                                                              &
     (8.*(1 + mf2)**1.5*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (3*k**4*nff)/                                                               &
     (8.*(1 + mf2)**2.5*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (3*k**7*nbGluon*(-k - mu0 + Complex(0,1)*p0))/                              &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**4 -                 &
    (3*k**7*nbGluon*(k - mu0 + Complex(0,1)*p0))/                               &
     (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**4 +                  &
    (3*k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     (2.*(1 + mf2)**2*(-k**2 +                                                  &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**6*nfd1xa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                 &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (3*k**6*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/               &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**4) +                &
    (k**6*nfd1xf*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                    &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                     &
     (2.*(1 + mf2)**2*(-k**2 +                                                  &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**6*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/                  &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**4) -                   &
    k**6/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3) -       &
    k**6/(2.*(1 + mf2)**1.5*(-k**2 +                                            &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    (3*k**4)/(8.*(1 + mf2)**2.5*                                                &
       (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) +            &
    (3*k**7*(-k - mu0 + Complex(0,1)*p0c))/                                     &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**4 +                &
    (3*k**5*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     (2.*(1 + mf2)**2*(-k**2 +                                                  &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    (3*k**6*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)/                     &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**4))/2.                

  b2f3aA=real(b2f3aAC)



end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine BF_AT_thr(mf2,T,k)
!Calculating the right hand side of differential equations

  implicit none


  real(8) mf2,T,k
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) mu0,p0,p0c
  real(8) nbGluon,nbd1xGluon,nbd2xGluon,nbd3xGluon,nbd4xGluon,nbd5xGluon
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa

  real(8) b1f1AT,b2f1aAT,b1f2AT,b2f2aAT,b1f3AT,b3f1aAT,b3f2aAT,b2f3aAT
  complex(8) b1f1AC,b2f1aAC,b1f2AC,b2f2aAC,b1f3AC,b3f1aAC,b3f2aAC,b2f3aAC

  common /mu0p0_com/ mu0,p0,p0c
  common /nbGluon_com/ nbGluon,nbd1xGluon,nbd2xGluon,nbd3xGluon,nbd4xGluon,nbd5xGluon
  common /nfdnx_com/ nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  common /BF_AT_thr_com/ b1f1AT,b2f1aAT,b1f2AT,b2f2aAT,b1f3AT,b3f1aAT,b3f2aAT,b2f3aAT


  b1f1AC=-(k**2*nbGluon)/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) - &
  (k**2*nbGluon)/(2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) +    &
  (k**2*nfa)/(2.*Sqrt(1 + mf2)*                                                 &
     (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) +               &
  (k**2*nff)/(2.*Sqrt(1 + mf2)*                                                 &
     (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -                  &
  (0.)*k**2/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)) -            &
  (0.)*k**2/(2.*Sqrt(1 + mf2)*(-k**2 +                                               &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2))                          

  b1f1AT=real(b1f1AC)


  b2f1aAC=(k**3*nbd1xGluon)/                                                    &
   (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) -                 &
  (k**2*nbGluon)/                                                               &
   (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) +                 &
  (k**3*nbd1xGluon)/                                                            &
   (4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) -                  &
  (k**2*nbGluon)/(4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) -    &
  (k**4*nfa)/(2.*Sqrt(1 + mf2)*                                                 &
     (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) -            &
  (k**4*nff)/(2.*Sqrt(1 + mf2)*                                                 &
     (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +               &
  (k**3*nbGluon*(-k - mu0 + Complex(0,1)*p0))/                                  &
   (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) -              &
  (k**3*nbGluon*(k - mu0 + Complex(0,1)*p0))/                                   &
   (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) -               &
  (0.)*k**2/(4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)) +            &
  (0.)*k**4/(2.*Sqrt(1 + mf2)*(-k**2 +                                               &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) +                    &
  (0.)*(k**3*(-k - mu0 + Complex(0,1)*p0c))/                                         &
   (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**2)               

  b2f1aAT=real(b2f1aAC)
  

  b1f2AC=(k**4*nbGluon)/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)** &
      2) + (k**4*nbGluon)/                                                      &
   (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +               &
  (k**2*nfa)/(4.*(1 + mf2)**1.5*                                                &
     (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -               &
  (k**3*nfd1xa)/                                                                &
   (4.*(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2))     &
- (k**3*nfd1xf)/                                                                &
   (4.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +      &
  (k**2*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -                  &
  (k**3*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   (2.*(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) + (k**3*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                &
   (2.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
+ (0.)*k**4/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**2) -         &
  (0.)*k**2/(4.*(1 + mf2)**1.5*(-k**2 +                                              &
       (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)) -                        &
  (0.)*(k**3*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                            &
   (2.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2)    

  b1f2AT=real(b1f2AC)


  b2f2aAC=-(k**5*nbd1xGluon)/                                                   &
   (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) +              &
  (k**4*nbGluon)/                                                               &
   (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) -              &
  (k**5*nbd1xGluon)/                                                            &
   (4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +               &
  (k**4*nbGluon)/                                                               &
   (4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) -               &
  (k**4*nfa)/(4.*(1 + mf2)**1.5*                                                &
     (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +            &
  (k**5*nfd1xa)/                                                                &
   (4.*(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**    &
      2) + (k**5*nfd1xf)/                                                       &
   (4.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2)     &
- (k**4*nff)/(4.*(1 + mf2)**1.5*                                                &
     (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -               &
  (k**5*nbGluon*(-k - mu0 + Complex(0,1)*p0))/                                  &
   (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3 +                   &
  (k**5*nbGluon*(k - mu0 + Complex(0,1)*p0))/                                   &
   (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3 +                    &
  (k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                      &
   ((1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3)     &
- (k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                         &
   ((1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +      &
  (0.)*k**4/(4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**2) +         &
  (0.)*k**4/(4.*(1 + mf2)**1.5*(-k**2 +                                              &
        (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                    &
  (0.)*(k**5*(-k - mu0 + Complex(0,1)*p0c))/                                         &
   (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3 +                  &
  (0.)*(k**5*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                            &
   ((1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3)       

  b2f2aAT=real(b2f2aAC)


  b1f3AC=(-((k**6*nbGluon)/                                                     &
       (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3) -              &
    (k**6*nbGluon)/                                                             &
     (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3 -                  &
    (k**4*nfa)/(4.*(1 + mf2)**1.5*                                              &
       (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +          &
    (3*k**2*nfa)/                                                               &
     (8.*(1 + mf2)**2.5*(-k**2 +                                                &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (3*k**3*nfd1xa)/                                                            &
     (8.*(1 + mf2)**2*(-k**2 +                                                  &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) +                    &
    (k**4*nfd2xa)/                                                              &
     (8.*(1 + mf2)**1.5*(-k**2 +                                                &
         (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)) -                    &
    (k**4*nff)/(4.*(1 + mf2)**1.5*                                              &
       (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -             &
    (3*k**3*nfd1xf)/                                                            &
     (8.*(1 + mf2)**2*(-k**2 +                                                  &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +                       &
    (k**4*nfd2xf)/                                                              &
     (8.*(1 + mf2)**1.5*(-k**2 +                                                &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) +                       &
    (3*k**2*nff)/                                                               &
     (8.*(1 + mf2)**2.5*(-k**2 +                                                &
         (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)) -                       &
    (3*k**3*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     (4.*(1 + mf2)**2*(-k**2 +                                                  &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**4*nfd1xa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                 &
     (2.*(1 + mf2)**1.5*(-k**2 +                                                &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**4*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/                 &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**4*nfd1xf*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                    &
     (2.*(1 + mf2)**1.5*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (3*k**3*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                     &
     (4.*(1 + mf2)**2*(-k**2 +                                                  &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (k**4*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/                    &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (0.)*k**6/(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3 +            &
    (0.)*k**4/(4.*(1 + mf2)**1.5*(-k**2 +                                            &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                  &
    (0.)*(3*k**2)/(8.*(1 + mf2)**2.5*                                                &
       (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)) -               &
    (0.)*(3*k**3*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     (4.*(1 + mf2)**2*(-k**2 +                                                  &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) -                  &
    (0.)*(k**4*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)/                       &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3))/2.                

  b1f3AT=real(b1f3AC)


  b3f1aAC=((k**4*nbGluon)/                                                      &
     (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) +            &
    (3*k**3*nbd1xGluon)/                                                        &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) -               &
    (k**4*nbd2xGluon)/                                                          &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) -               &
    (3*k**2*nbGluon)/                                                           &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)) +               &
    (k**4*nbGluon)/                                                             &
     (4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +             &
    (3*k**3*nbd1xGluon)/                                                        &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) -                &
    (k**4*nbd2xGluon)/                                                          &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) -                &
    (3*k**2*nbGluon)/                                                           &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)) +                &
    (k**6*nfa)/(Sqrt(1 + mf2)*(-k**2 +                                          &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) +                &
    (k**6*nff)/(Sqrt(1 + mf2)*(-k**2 +                                          &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (k**4*nbd1xGluon*(-k - mu0 + Complex(0,1)*p0))/                             &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) +            &
    (3*k**3*nbGluon*(-k - mu0 + Complex(0,1)*p0))/                              &
     (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) -            &
    (k**4*nbGluon*(-k - mu0 + Complex(0,1)*p0)**2)/                             &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3 +                 &
    (k**4*nbd1xGluon*(k - mu0 + Complex(0,1)*p0))/                              &
     (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) -             &
    (3*k**3*nbGluon*(k - mu0 + Complex(0,1)*p0))/                               &
     (4.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) -             &
    (k**4*nbGluon*(k - mu0 + Complex(0,1)*p0)**2)/                              &
     (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3 +                  &
    (0.)*k**4/(4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**2) -       &
    (0.)*(3*k**2)/(8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)) -      &
    (0.)*k**6/(Sqrt(1 + mf2)*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    (0.)*(3*k**3*(-k - mu0 + Complex(0,1)*p0c))/                                     &
     (4.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**2) -           &
    (0.)*(k**4*(-k - mu0 + Complex(0,1)*p0c)**2)/                                    &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3)/2.              

  b3f1aAT=real(b3f1aAC)


  b3f2aAC=(-(k**6*nbGluon)/                                                     &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3) -            &
    (3*k**5*nbd1xGluon)/                                                        &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) +            &
    (k**6*nbd2xGluon)/                                                          &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) +            &
    (3*k**4*nbGluon)/                                                           &
     (8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**2) -            &
    (k**6*nbGluon)/                                                             &
     (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3) -             &
    (3*k**5*nbd1xGluon)/                                                        &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +             &
    (k**6*nbd2xGluon)/                                                          &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +             &
    (3*k**4*nbGluon)/                                                           &
     (8.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**2) +             &
    (k**6*nfa)/(2.*(1 + mf2)**1.5*                                              &
       (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -          &
    (k**7*nfd1xa)/                                                              &
     (2.*(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**      &
           2)**3) - (k**7*nfd1xf)/                                              &
     (2.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**     &
        3) + (k**6*nff)/                                                        &
     (2.*(1 + mf2)**1.5*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +                   &
    (k**6*nbd1xGluon*(-k - mu0 + Complex(0,1)*p0))/                             &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3 -                 &
    (3*k**5*nbGluon*(-k - mu0 + Complex(0,1)*p0))/                              &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3) +            &
    (3*k**6*nbGluon*(-k - mu0 + Complex(0,1)*p0)**2)/                           &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**4 -                 &
    (k**6*nbd1xGluon*(k - mu0 + Complex(0,1)*p0))/                              &
     (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3 +                  &
    (3*k**5*nbGluon*(k - mu0 + Complex(0,1)*p0))/                               &
     (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3) +             &
    (3*k**6*nbGluon*(k - mu0 + Complex(0,1)*p0)**2)/                            &
     (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**4 -                  &
    (3*k**7*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     ((1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**     &
        4) + (3*k**7*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/            &
     ((1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**4)      &
- (0.)*k**6/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3) +         &
    (0.)*(3*k**4)/(8.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**       &
        2) - (0.)*k**6/                                                              &
     (2.*(1 + mf2)**1.5*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) -                  &
    (0.)*(3*k**5*(-k - mu0 + Complex(0,1)*p0c))/                                     &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3) +           &
    (0.)*(3*k**6*(-k - mu0 + Complex(0,1)*p0c)**2)/                                  &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**4 -                &
    (0.)*(3*k**7*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     ((1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**4))/   &
  2.                                                                            

  b3f2aAT=real(b3f2aAC)


  b2f3aAC=((k**7*nbd1xGluon)/                                                   &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3) -            &
    (k**6*nbGluon)/                                                             &
     (2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**3) +            &
    (k**7*nbd1xGluon)/                                                          &
     (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3) -             &
    (k**6*nbGluon)/                                                             &
     (2.*(-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**3) +             &
    (k**6*nfa)/(2.*(1 + mf2)**1.5*                                              &
       (-k**2 + (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -          &
    (3*k**4*nfa)/                                                               &
     (8.*(1 + mf2)**2.5*(-k**2 +                                                &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (3*k**5*nfd1xa)/                                                            &
     (8.*(1 + mf2)**2*(-k**2 +                                                  &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) -                &
    (k**6*nfd2xa)/                                                              &
     (8.*(1 + mf2)**1.5*(-k**2 +                                                &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**2) +                &
    (k**6*nff)/(2.*(1 + mf2)**1.5*                                              &
       (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) +             &
    (3*k**5*nfd1xf)/                                                            &
     (8.*(1 + mf2)**2*(-k**2 +                                                  &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (k**6*nfd2xf)/                                                              &
     (8.*(1 + mf2)**1.5*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) -                   &
    (3*k**4*nff)/                                                               &
     (8.*(1 + mf2)**2.5*(-k**2 +                                                &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**2) +                   &
    (3*k**7*nbGluon*(-k - mu0 + Complex(0,1)*p0))/                              &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0)**2)**4 -                 &
    (3*k**7*nbGluon*(k - mu0 + Complex(0,1)*p0))/                               &
     (-(k**2*(1 + mf2)) + (k - mu0 + Complex(0,1)*p0)**2)**4 +                  &
    (3*k**5*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                  &
     (2.*(1 + mf2)**2*(-k**2 +                                                  &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (k**6*nfd1xa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0))/                 &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**3) -                &
    (3*k**6*nfa*(-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)/               &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (-(k*Sqrt(1 + mf2)) - mu0 + Complex(0,1)*p0)**2)**4) +                &
    (k**6*nfd1xf*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                    &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**5*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0))/                     &
     (2.*(1 + mf2)**2*(-k**2 +                                                  &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**3) -                   &
    (3*k**6*nff*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)/                  &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0)**2)**4) -                   &
    (0.)*k**6/(2.*(-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**3) -       &
    (0.)*k**6/(2.*(1 + mf2)**1.5*(-k**2 +                                            &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    (0.)*(3*k**4)/(8.*(1 + mf2)**2.5*                                                &
       (-k**2 + (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**2) +            &
    (0.)*(3*k**7*(-k - mu0 + Complex(0,1)*p0c))/                                     &
     (-(k**2*(1 + mf2)) + (-k - mu0 + Complex(0,1)*p0c)**2)**4 +                &
    (0.)*(3*k**5*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c))/                        &
     (2.*(1 + mf2)**2*(-k**2 +                                                  &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**3) +                  &
    (0.)*(3*k**6*(k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)/                     &
     ((1 + mf2)**1.5*(-k**2 +                                                   &
          (k*Sqrt(1 + mf2) - mu0 + Complex(0,1)*p0c)**2)**4))/2.                

  b2f3aAT=real(b2f3aAC)



end
