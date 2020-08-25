subroutine lbtem(k,dmbo2drho,nbdnx,etaphi,lbt)
!input dimensionless masses


  implicit none
  integer N 
  !The number of mesons
  parameter(N=2)
  real(16) k
  real(16) dmbo2drho(N,12),nbdnx(N,6),lbt(N,12)
  real(16) mb2d00rho(N),mb2d10rho(N),mb2d20rho(N),mb2d30rho(N),mb2d40rho(N),mb2d50rho(N)
  real(16) mb2d01rho(N),mb2d11rho(N),mb2d21rho(N),mb2d31rho(N)
  real(16) mb2d02rho(N),mb2d12rho(N)
  real(16) nbd0x(N),nbd1x(N),nbd2x(N),nbd3x(N),nbd4x(N),nbd5x(N)
  real(16) etaphi
  real(16) lbt00(N),lbt10(N),lbt20(N),lbt30(N),lbt40(N),lbt50(N)
  real(16) lbt01(N),lbt11(N),lbt21(N),lbt31(N)
  real(16) lbt02(N),lbt12(N)

  mb2d00rho=dmbo2drho(:,1)
  mb2d10rho=dmbo2drho(:,2)
  mb2d20rho=dmbo2drho(:,3)
  mb2d30rho=dmbo2drho(:,4)
  mb2d40rho=dmbo2drho(:,5)
  mb2d50rho=dmbo2drho(:,6)
  mb2d01rho=dmbo2drho(:,7)
  mb2d11rho=dmbo2drho(:,8)
  mb2d21rho=dmbo2drho(:,9)
  mb2d31rho=dmbo2drho(:,10)
  mb2d02rho=dmbo2drho(:,11)
  mb2d12rho=dmbo2drho(:,12)

  nbd0x=nbdnx(:,1)
  nbd1x=nbdnx(:,2)
  nbd2x=nbdnx(:,3)
  nbd3x=nbdnx(:,4)
  nbd4x=nbdnx(:,5)
  nbd5x=nbdnx(:,6)

  lbt00=(2*(1 - etaphi/5.)*(0.5 + nbd0x))/(3.*Sqrt(1 + mb2d00rho))
  lbt10=-((1 - etaphi/5.)*mb2d10rho*(0.5 + nbd0x))/                             &
     (3.*(1 + mb2d00rho)**1.5) +                                                &
    ((1 - etaphi/5.)*k*mb2d10rho*nbd1x)/(3.*(1 + mb2d00rho))
  lbt20=((1 - etaphi/5.)*mb2d10rho**2*(0.5 + nbd0x))/                           &
     (2.*(1 + mb2d00rho)**2.5) -                                                &
    ((1 - etaphi/5.)*mb2d20rho*(0.5 + nbd0x))/(3.*(1 + mb2d00rho)**1.5) -       &
    ((1 - etaphi/5.)*k*mb2d10rho**2*nbd1x)/(2.*(1 + mb2d00rho)**2) +            &
    ((1 - etaphi/5.)*k*mb2d20rho*nbd1x)/(3.*(1 + mb2d00rho)) +                  &
    ((1 - etaphi/5.)*k**2*mb2d10rho**2*nbd2x)/(6.*(1 + mb2d00rho)**1.5)
  lbt30=(-5*(1 - etaphi/5.)*mb2d10rho**3*(0.5 + nbd0x))/                        &
     (4.*(1 + mb2d00rho)**3.5) +                                                &
    (3*(1 - etaphi/5.)*mb2d10rho*mb2d20rho*(0.5 + nbd0x))/                      &
     (2.*(1 + mb2d00rho)**2.5) -                                                &
    ((1 - etaphi/5.)*mb2d30rho*(0.5 + nbd0x))/(3.*(1 + mb2d00rho)**1.5) +       &
    (5*(1 - etaphi/5.)*k*mb2d10rho**3*nbd1x)/(4.*(1 + mb2d00rho)**3) -          &
    (3*(1 - etaphi/5.)*k*mb2d10rho*mb2d20rho*nbd1x)/                            &
     (2.*(1 + mb2d00rho)**2) +                                                  &
    ((1 - etaphi/5.)*k*mb2d30rho*nbd1x)/(3.*(1 + mb2d00rho)) -                  &
    ((1 - etaphi/5.)*k**2*mb2d10rho**3*nbd2x)/(2.*(1 + mb2d00rho)**2.5) +       &
    ((1 - etaphi/5.)*k**2*mb2d10rho*mb2d20rho*nbd2x)/                           &
     (2.*(1 + mb2d00rho)**1.5) +                                                &
    ((1 - etaphi/5.)*k**3*mb2d10rho**3*nbd3x)/(12.*(1 + mb2d00rho)**2)
  lbt40=(35*(1 - etaphi/5.)*mb2d10rho**4*(0.5 + nbd0x))/                        &
     (8.*(1 + mb2d00rho)**4.5) -                                                &
    (15*(1 - etaphi/5.)*mb2d10rho**2*mb2d20rho*(0.5 + nbd0x))/                  &
     (2.*(1 + mb2d00rho)**3.5) +                                                &
    (3*(1 - etaphi/5.)*mb2d20rho**2*(0.5 + nbd0x))/                             &
     (2.*(1 + mb2d00rho)**2.5) +                                                &
    (2*(1 - etaphi/5.)*mb2d10rho*mb2d30rho*(0.5 + nbd0x))/                      &
     (1 + mb2d00rho)**2.5 - ((1 - etaphi/5.)*mb2d40rho*(0.5 + nbd0x))/          &
     (3.*(1 + mb2d00rho)**1.5) -                                                &
    (35*(1 - etaphi/5.)*k*mb2d10rho**4*nbd1x)/(8.*(1 + mb2d00rho)**4) +         &
    (15*(1 - etaphi/5.)*k*mb2d10rho**2*mb2d20rho*nbd1x)/                        &
     (2.*(1 + mb2d00rho)**3) -                                                  &
    (3*(1 - etaphi/5.)*k*mb2d20rho**2*nbd1x)/(2.*(1 + mb2d00rho)**2) -          &
    (2*(1 - etaphi/5.)*k*mb2d10rho*mb2d30rho*nbd1x)/(1 + mb2d00rho)**2 +        &
    ((1 - etaphi/5.)*k*mb2d40rho*nbd1x)/(3.*(1 + mb2d00rho)) +                  &
    (15*(1 - etaphi/5.)*k**2*mb2d10rho**4*nbd2x)/                               &
     (8.*(1 + mb2d00rho)**3.5) -                                                &
    (3*(1 - etaphi/5.)*k**2*mb2d10rho**2*mb2d20rho*nbd2x)/                      &
     (1 + mb2d00rho)**2.5 + ((1 - etaphi/5.)*k**2*mb2d20rho**2*nbd2x)/          &
     (2.*(1 + mb2d00rho)**1.5) +                                                &
    (2*(1 - etaphi/5.)*k**2*mb2d10rho*mb2d30rho*nbd2x)/                         &
     (3.*(1 + mb2d00rho)**1.5) -                                                &
    (5*(1 - etaphi/5.)*k**3*mb2d10rho**4*nbd3x)/(12.*(1 + mb2d00rho)**3) +      &
    ((1 - etaphi/5.)*k**3*mb2d10rho**2*mb2d20rho*nbd3x)/                        &
     (2.*(1 + mb2d00rho)**2) +                                                  &
    ((1 - etaphi/5.)*k**4*mb2d10rho**4*nbd4x)/(24.*(1 + mb2d00rho)**2.5)
  lbt50=(-315*(1 - etaphi/5.)*mb2d10rho**5*(0.5 + nbd0x))/                      &
     (16.*(1 + mb2d00rho)**5.5) +                                               &
    (175*(1 - etaphi/5.)*mb2d10rho**3*mb2d20rho*(0.5 + nbd0x))/                 &
     (4.*(1 + mb2d00rho)**4.5) -                                                &
    (75*(1 - etaphi/5.)*mb2d10rho*mb2d20rho**2*(0.5 + nbd0x))/                  &
     (4.*(1 + mb2d00rho)**3.5) -                                                &
    (25*(1 - etaphi/5.)*mb2d10rho**2*mb2d30rho*(0.5 + nbd0x))/                  &
     (2.*(1 + mb2d00rho)**3.5) +                                                &
    (5*(1 - etaphi/5.)*mb2d20rho*mb2d30rho*(0.5 + nbd0x))/                      &
     (1 + mb2d00rho)**2.5 + (5*(1 - etaphi/5.)*mb2d10rho*mb2d40rho*             &
       (0.5 + nbd0x))/(2.*(1 + mb2d00rho)**2.5) -                               &
    ((1 - etaphi/5.)*mb2d50rho*(0.5 + nbd0x))/(3.*(1 + mb2d00rho)**1.5) +       &
    (315*(1 - etaphi/5.)*k*mb2d10rho**5*nbd1x)/(16.*(1 + mb2d00rho)**5) -       &
    (175*(1 - etaphi/5.)*k*mb2d10rho**3*mb2d20rho*nbd1x)/                       &
     (4.*(1 + mb2d00rho)**4) +                                                  &
    (75*(1 - etaphi/5.)*k*mb2d10rho*mb2d20rho**2*nbd1x)/                        &
     (4.*(1 + mb2d00rho)**3) +                                                  &
    (25*(1 - etaphi/5.)*k*mb2d10rho**2*mb2d30rho*nbd1x)/                        &
     (2.*(1 + mb2d00rho)**3) -                                                  &
    (5*(1 - etaphi/5.)*k*mb2d20rho*mb2d30rho*nbd1x)/(1 + mb2d00rho)**2 -        &
    (5*(1 - etaphi/5.)*k*mb2d10rho*mb2d40rho*nbd1x)/                            &
     (2.*(1 + mb2d00rho)**2) +                                                  &
    ((1 - etaphi/5.)*k*mb2d50rho*nbd1x)/(3.*(1 + mb2d00rho)) -                  &
    (35*(1 - etaphi/5.)*k**2*mb2d10rho**5*nbd2x)/                               &
     (4.*(1 + mb2d00rho)**4.5) +                                                &
    (75*(1 - etaphi/5.)*k**2*mb2d10rho**3*mb2d20rho*nbd2x)/                     &
     (4.*(1 + mb2d00rho)**3.5) -                                                &
    (15*(1 - etaphi/5.)*k**2*mb2d10rho*mb2d20rho**2*nbd2x)/                     &
     (2.*(1 + mb2d00rho)**2.5) -                                                &
    (5*(1 - etaphi/5.)*k**2*mb2d10rho**2*mb2d30rho*nbd2x)/                      &
     (1 + mb2d00rho)**2.5 + (5*(1 - etaphi/5.)*k**2*mb2d20rho*mb2d30rho*        &
       nbd2x)/(3.*(1 + mb2d00rho)**1.5) +                                       &
    (5*(1 - etaphi/5.)*k**2*mb2d10rho*mb2d40rho*nbd2x)/                         &
     (6.*(1 + mb2d00rho)**1.5) +                                                &
    (35*(1 - etaphi/5.)*k**3*mb2d10rho**5*nbd3x)/(16.*(1 + mb2d00rho)**4) -     &
    (25*(1 - etaphi/5.)*k**3*mb2d10rho**3*mb2d20rho*nbd3x)/                     &
     (6.*(1 + mb2d00rho)**3) +                                                  &
    (5*(1 - etaphi/5.)*k**3*mb2d10rho*mb2d20rho**2*nbd3x)/                      &
     (4.*(1 + mb2d00rho)**2) +                                                  &
    (5*(1 - etaphi/5.)*k**3*mb2d10rho**2*mb2d30rho*nbd3x)/                      &
     (6.*(1 + mb2d00rho)**2) -                                                  &
    (5*(1 - etaphi/5.)*k**4*mb2d10rho**5*nbd4x)/                                &
     (16.*(1 + mb2d00rho)**3.5) +                                               &
    (5*(1 - etaphi/5.)*k**4*mb2d10rho**3*mb2d20rho*nbd4x)/                      &
     (12.*(1 + mb2d00rho)**2.5) +                                               &
    ((1 - etaphi/5.)*k**5*mb2d10rho**5*nbd5x)/(48.*(1 + mb2d00rho)**3)
  lbt01=-((1 - etaphi/5.)*mb2d01rho*(0.5 + nbd0x))/                             &
     (3.*(1 + mb2d00rho)**1.5) +                                                &
    ((1 - etaphi/5.)*k*mb2d01rho*nbd1x)/(3.*(1 + mb2d00rho))
  lbt11=((1 - etaphi/5.)*mb2d01rho*mb2d10rho*(0.5 + nbd0x))/                    &
     (2.*(1 + mb2d00rho)**2.5) -                                                &
    ((1 - etaphi/5.)*mb2d11rho*(0.5 + nbd0x))/(3.*(1 + mb2d00rho)**1.5) -       &
    ((1 - etaphi/5.)*k*mb2d01rho*mb2d10rho*nbd1x)/(2.*(1 + mb2d00rho)**2) +     &
    ((1 - etaphi/5.)*k*mb2d11rho*nbd1x)/(3.*(1 + mb2d00rho)) +                  &
    ((1 - etaphi/5.)*k**2*mb2d01rho*mb2d10rho*nbd2x)/                           &
     (6.*(1 + mb2d00rho)**1.5)
  lbt21=(-5*(1 - etaphi/5.)*mb2d01rho*mb2d10rho**2*(0.5 + nbd0x))/              &
     (4.*(1 + mb2d00rho)**3.5) +                                                &
    ((1 - etaphi/5.)*mb2d10rho*mb2d11rho*(0.5 + nbd0x))/                        &
     (1 + mb2d00rho)**2.5 + ((1 - etaphi/5.)*mb2d01rho*mb2d20rho*               &
       (0.5 + nbd0x))/(2.*(1 + mb2d00rho)**2.5) -                               &
    ((1 - etaphi/5.)*mb2d21rho*(0.5 + nbd0x))/(3.*(1 + mb2d00rho)**1.5) +       &
    (5*(1 - etaphi/5.)*k*mb2d01rho*mb2d10rho**2*nbd1x)/                         &
     (4.*(1 + mb2d00rho)**3) -                                                  &
    ((1 - etaphi/5.)*k*mb2d10rho*mb2d11rho*nbd1x)/(1 + mb2d00rho)**2 -          &
    ((1 - etaphi/5.)*k*mb2d01rho*mb2d20rho*nbd1x)/(2.*(1 + mb2d00rho)**2) +     &
    ((1 - etaphi/5.)*k*mb2d21rho*nbd1x)/(3.*(1 + mb2d00rho)) -                  &
    ((1 - etaphi/5.)*k**2*mb2d01rho*mb2d10rho**2*nbd2x)/                        &
     (2.*(1 + mb2d00rho)**2.5) +                                                &
    ((1 - etaphi/5.)*k**2*mb2d10rho*mb2d11rho*nbd2x)/                           &
     (3.*(1 + mb2d00rho)**1.5) +                                                &
    ((1 - etaphi/5.)*k**2*mb2d01rho*mb2d20rho*nbd2x)/                           &
     (6.*(1 + mb2d00rho)**1.5) +                                                &
    ((1 - etaphi/5.)*k**3*mb2d01rho*mb2d10rho**2*nbd3x)/                        &
     (12.*(1 + mb2d00rho)**2)
  lbt31=(35*(1 - etaphi/5.)*mb2d01rho*mb2d10rho**3*(0.5 + nbd0x))/              &
     (8.*(1 + mb2d00rho)**4.5) -                                                &
    (15*(1 - etaphi/5.)*mb2d10rho**2*mb2d11rho*(0.5 + nbd0x))/                  &
     (4.*(1 + mb2d00rho)**3.5) -                                                &
    (15*(1 - etaphi/5.)*mb2d01rho*mb2d10rho*mb2d20rho*(0.5 + nbd0x))/           &
     (4.*(1 + mb2d00rho)**3.5) +                                                &
    (3*(1 - etaphi/5.)*mb2d11rho*mb2d20rho*(0.5 + nbd0x))/                      &
     (2.*(1 + mb2d00rho)**2.5) +                                                &
    (3*(1 - etaphi/5.)*mb2d10rho*mb2d21rho*(0.5 + nbd0x))/                      &
     (2.*(1 + mb2d00rho)**2.5) +                                                &
    ((1 - etaphi/5.)*mb2d01rho*mb2d30rho*(0.5 + nbd0x))/                        &
     (2.*(1 + mb2d00rho)**2.5) -                                                &
    ((1 - etaphi/5.)*mb2d31rho*(0.5 + nbd0x))/(3.*(1 + mb2d00rho)**1.5) -       &
    (35*(1 - etaphi/5.)*k*mb2d01rho*mb2d10rho**3*nbd1x)/                        &
     (8.*(1 + mb2d00rho)**4) +                                                  &
    (15*(1 - etaphi/5.)*k*mb2d10rho**2*mb2d11rho*nbd1x)/                        &
     (4.*(1 + mb2d00rho)**3) +                                                  &
    (15*(1 - etaphi/5.)*k*mb2d01rho*mb2d10rho*mb2d20rho*nbd1x)/                 &
     (4.*(1 + mb2d00rho)**3) -                                                  &
    (3*(1 - etaphi/5.)*k*mb2d11rho*mb2d20rho*nbd1x)/                            &
     (2.*(1 + mb2d00rho)**2) -                                                  &
    (3*(1 - etaphi/5.)*k*mb2d10rho*mb2d21rho*nbd1x)/                            &
     (2.*(1 + mb2d00rho)**2) -                                                  &
    ((1 - etaphi/5.)*k*mb2d01rho*mb2d30rho*nbd1x)/(2.*(1 + mb2d00rho)**2) +     &
    ((1 - etaphi/5.)*k*mb2d31rho*nbd1x)/(3.*(1 + mb2d00rho)) +                  &
    (15*(1 - etaphi/5.)*k**2*mb2d01rho*mb2d10rho**3*nbd2x)/                     &
     (8.*(1 + mb2d00rho)**3.5) -                                                &
    (3*(1 - etaphi/5.)*k**2*mb2d10rho**2*mb2d11rho*nbd2x)/                      &
     (2.*(1 + mb2d00rho)**2.5) -                                                &
    (3*(1 - etaphi/5.)*k**2*mb2d01rho*mb2d10rho*mb2d20rho*nbd2x)/               &
     (2.*(1 + mb2d00rho)**2.5) +                                                &
    ((1 - etaphi/5.)*k**2*mb2d11rho*mb2d20rho*nbd2x)/                           &
     (2.*(1 + mb2d00rho)**1.5) +                                                &
    ((1 - etaphi/5.)*k**2*mb2d10rho*mb2d21rho*nbd2x)/                           &
     (2.*(1 + mb2d00rho)**1.5) +                                                &
    ((1 - etaphi/5.)*k**2*mb2d01rho*mb2d30rho*nbd2x)/                           &
     (6.*(1 + mb2d00rho)**1.5) -                                                &
    (5*(1 - etaphi/5.)*k**3*mb2d01rho*mb2d10rho**3*nbd3x)/                      &
     (12.*(1 + mb2d00rho)**3) +                                                 &
    ((1 - etaphi/5.)*k**3*mb2d10rho**2*mb2d11rho*nbd3x)/                        &
     (4.*(1 + mb2d00rho)**2) +                                                  &
    ((1 - etaphi/5.)*k**3*mb2d01rho*mb2d10rho*mb2d20rho*nbd3x)/                 &
     (4.*(1 + mb2d00rho)**2) +                                                  &
    ((1 - etaphi/5.)*k**4*mb2d01rho*mb2d10rho**3*nbd4x)/                        &
     (24.*(1 + mb2d00rho)**2.5)
  lbt02=((1 - etaphi/5.)*mb2d01rho**2*(0.5 + nbd0x))/                           &
     (2.*(1 + mb2d00rho)**2.5) -                                                &
    ((1 - etaphi/5.)*mb2d02rho*(0.5 + nbd0x))/(3.*(1 + mb2d00rho)**1.5) -       &
    ((1 - etaphi/5.)*k*mb2d01rho**2*nbd1x)/(2.*(1 + mb2d00rho)**2) +            &
    ((1 - etaphi/5.)*k*mb2d02rho*nbd1x)/(3.*(1 + mb2d00rho)) +                  &
    ((1 - etaphi/5.)*k**2*mb2d01rho**2*nbd2x)/(6.*(1 + mb2d00rho)**1.5)
  lbt12=(-5*(1 - etaphi/5.)*mb2d01rho**2*mb2d10rho*(0.5 + nbd0x))/              &
     (4.*(1 + mb2d00rho)**3.5) +                                                &
    ((1 - etaphi/5.)*mb2d02rho*mb2d10rho*(0.5 + nbd0x))/                        &
     (2.*(1 + mb2d00rho)**2.5) +                                                &
    ((1 - etaphi/5.)*mb2d01rho*mb2d11rho*(0.5 + nbd0x))/                        &
     (1 + mb2d00rho)**2.5 - ((1 - etaphi/5.)*mb2d12rho*(0.5 + nbd0x))/          &
     (3.*(1 + mb2d00rho)**1.5) +                                                &
    (5*(1 - etaphi/5.)*k*mb2d01rho**2*mb2d10rho*nbd1x)/                         &
     (4.*(1 + mb2d00rho)**3) -                                                  &
    ((1 - etaphi/5.)*k*mb2d02rho*mb2d10rho*nbd1x)/(2.*(1 + mb2d00rho)**2) -     &
    ((1 - etaphi/5.)*k*mb2d01rho*mb2d11rho*nbd1x)/(1 + mb2d00rho)**2 +          &
    ((1 - etaphi/5.)*k*mb2d12rho*nbd1x)/(3.*(1 + mb2d00rho)) -                  &
    ((1 - etaphi/5.)*k**2*mb2d01rho**2*mb2d10rho*nbd2x)/                        &
     (2.*(1 + mb2d00rho)**2.5) +                                                &
    ((1 - etaphi/5.)*k**2*mb2d02rho*mb2d10rho*nbd2x)/                           &
     (6.*(1 + mb2d00rho)**1.5) +                                                &
    ((1 - etaphi/5.)*k**2*mb2d01rho*mb2d11rho*nbd2x)/                           &
     (3.*(1 + mb2d00rho)**1.5) +                                                &
    ((1 - etaphi/5.)*k**3*mb2d01rho**2*mb2d10rho*nbd3x)/                        &
     (12.*(1 + mb2d00rho)**2)

  lbt(:,1)=lbt00
  lbt(:,2)=lbt10
  lbt(:,3)=lbt20
  lbt(:,4)=lbt30
  lbt(:,5)=lbt40
  lbt(:,6)=lbt50
  lbt(:,7)=lbt01
  lbt(:,8)=lbt11
  lbt(:,9)=lbt21
  lbt(:,10)=lbt31
  lbt(:,11)=lbt02
  lbt(:,12)=lbt12
end

subroutine lftem(k,dmfe2drho,nffdnx,nfadnx,etapsi,lft)

  implicit none
  integer N 
  !The number of quarks
  parameter(N=2)
  real(16) k
  real(16) dmfe2drho(N,4),nffdnx(N,6),nfadnx(N,6),lft(N,12)
  real(16) etapsi
  real(16) mf2d00rho(N),mf2d10rho(N),mf2d01rho(N),mf2d02rho(N)
  real(16) nffd0x(N),nffd1x(N),nffd2x(N),nffd3x(N),nffd4x(N),nffd5x(N)
  real(16) nfad0x(N),nfad1x(N),nfad2x(N),nfad3x(N),nfad4x(N),nfad5x(N)
  real(16) lft00(N),lft10(N),lft20(N),lft30(N),lft40(N),lft50(N)
  real(16) lft01(N),lft11(N),lft21(N),lft31(N)
  real(16) lft02(N),lft12(N)

  mf2d00rho=dmfe2drho(:,1)
  mf2d10rho=dmfe2drho(:,2)
  mf2d01rho=dmfe2drho(:,3)
  mf2d02rho=dmfe2drho(:,4)

  nffd0x=nffdnx(:,1)
  nffd1x=nffdnx(:,2)
  nffd2x=nffdnx(:,3)
  nffd3x=nffdnx(:,4)
  nffd4x=nffdnx(:,5)
  nffd5x=nffdnx(:,6)

  nfad0x=nfadnx(:,1)
  nfad1x=nfadnx(:,2)
  nfad2x=nfadnx(:,3)
  nfad3x=nfadnx(:,4)
  nfad4x=nfadnx(:,5)
  nfad5x=nfadnx(:,6)

  lft00=((1 - etapsi/4.)*(1 - nfad0x - nffd0x))/                                &
    (3.*Sqrt(1 + mf2d00rho))
  lft10= -((1 - etapsi/4.)*mf2d10rho*(1 - nfad0x - nffd0x))/                    &
     (6.*(1 + mf2d00rho)**1.5) +                                                &
    ((1 - etapsi/4.)*(-(k*mf2d10rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -          &
         (k*mf2d10rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (3.*Sqrt(1 + mf2d00rho))
  lft20=((1 - etapsi/4.)*mf2d10rho**2*(1 - nfad0x - nffd0x))/                   &
     (4.*(1 + mf2d00rho)**2.5) -                                                &
    ((1 - etapsi/4.)*mf2d10rho*                                                 &
       (-(k*mf2d10rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d10rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (3.*(1 + mf2d00rho)**1.5) +                                                &
    ((1 - etapsi/4.)*((k*mf2d10rho**2*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -       &
         (k**2*mf2d10rho**2*nfad2x)/(4.*(1 + mf2d00rho)) +                      &
         (k*mf2d10rho**2*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -                    &
         (k**2*mf2d10rho**2*nffd2x)/(4.*(1 + mf2d00rho))))/                     &
     (3.*Sqrt(1 + mf2d00rho))
  lft30=(-5*(1 - etapsi/4.)*mf2d10rho**3*(1 - nfad0x - nffd0x))/                &
     (8.*(1 + mf2d00rho)**3.5) +                                                &
    (3*(1 - etapsi/4.)*mf2d10rho**2*                                            &
       (-(k*mf2d10rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d10rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (4.*(1 + mf2d00rho)**2.5) -                                                &
    ((1 - etapsi/4.)*mf2d10rho*                                                 &
       ((k*mf2d10rho**2*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -                     &
         (k**2*mf2d10rho**2*nfad2x)/(4.*(1 + mf2d00rho)) +                      &
         (k*mf2d10rho**2*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -                    &
         (k**2*mf2d10rho**2*nffd2x)/(4.*(1 + mf2d00rho))))/                     &
     (2.*(1 + mf2d00rho)**1.5) +                                                &
    ((1 - etapsi/4.)*((-3*k*mf2d10rho**3*nfad1x)/                               &
          (8.*(1 + mf2d00rho)**2.5) +                                           &
         (3*k**2*mf2d10rho**3*nfad2x)/(8.*(1 + mf2d00rho)**2) -                 &
         (k**3*mf2d10rho**3*nfad3x)/(8.*(1 + mf2d00rho)**1.5) -                 &
         (3*k*mf2d10rho**3*nffd1x)/(8.*(1 + mf2d00rho)**2.5) +                  &
         (3*k**2*mf2d10rho**3*nffd2x)/(8.*(1 + mf2d00rho)**2) -                 &
         (k**3*mf2d10rho**3*nffd3x)/(8.*(1 + mf2d00rho)**1.5)))/                &
     (3.*Sqrt(1 + mf2d00rho))
  lft40= (35*(1 - etapsi/4.)*mf2d10rho**4*(1 - nfad0x - nffd0x))/               &
     (16.*(1 + mf2d00rho)**4.5) -                                               &
    (5*(1 - etapsi/4.)*mf2d10rho**3*                                            &
       (-(k*mf2d10rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d10rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (2.*(1 + mf2d00rho)**3.5) +                                                &
    (3*(1 - etapsi/4.)*mf2d10rho**2*                                            &
       ((k*mf2d10rho**2*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -                     &
         (k**2*mf2d10rho**2*nfad2x)/(4.*(1 + mf2d00rho)) +                      &
         (k*mf2d10rho**2*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -                    &
         (k**2*mf2d10rho**2*nffd2x)/(4.*(1 + mf2d00rho))))/                     &
     (2.*(1 + mf2d00rho)**2.5) -                                                &
    (2*(1 - etapsi/4.)*mf2d10rho*                                               &
       ((-3*k*mf2d10rho**3*nfad1x)/(8.*(1 + mf2d00rho)**2.5) +                  &
         (3*k**2*mf2d10rho**3*nfad2x)/(8.*(1 + mf2d00rho)**2) -                 &
         (k**3*mf2d10rho**3*nfad3x)/(8.*(1 + mf2d00rho)**1.5) -                 &
         (3*k*mf2d10rho**3*nffd1x)/(8.*(1 + mf2d00rho)**2.5) +                  &
         (3*k**2*mf2d10rho**3*nffd2x)/(8.*(1 + mf2d00rho)**2) -                 &
         (k**3*mf2d10rho**3*nffd3x)/(8.*(1 + mf2d00rho)**1.5)))/                &
     (3.*(1 + mf2d00rho)**1.5) +                                                &
    ((1 - etapsi/4.)*((15*k*mf2d10rho**4*nfad1x)/                               &
          (16.*(1 + mf2d00rho)**3.5) -                                          &
         (15*k**2*mf2d10rho**4*nfad2x)/(16.*(1 + mf2d00rho)**3) +               &
         (3*k**3*mf2d10rho**4*nfad3x)/(8.*(1 + mf2d00rho)**2.5) -               &
         (k**4*mf2d10rho**4*nfad4x)/(16.*(1 + mf2d00rho)**2) +                  &
         (15*k*mf2d10rho**4*nffd1x)/(16.*(1 + mf2d00rho)**3.5) -                &
         (15*k**2*mf2d10rho**4*nffd2x)/(16.*(1 + mf2d00rho)**3) +               &
         (3*k**3*mf2d10rho**4*nffd3x)/(8.*(1 + mf2d00rho)**2.5) -               &
         (k**4*mf2d10rho**4*nffd4x)/(16.*(1 + mf2d00rho)**2)))/                 &
     (3.*Sqrt(1 + mf2d00rho))
  lft50=(-315*(1 - etapsi/4.)*mf2d10rho**5*(1 - nfad0x - nffd0x))/              &
     (32.*(1 + mf2d00rho)**5.5) +                                               &
    (175*(1 - etapsi/4.)*mf2d10rho**4*                                          &
       (-(k*mf2d10rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d10rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (16.*(1 + mf2d00rho)**4.5) -                                               &
    (25*(1 - etapsi/4.)*mf2d10rho**3*                                           &
       ((k*mf2d10rho**2*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -                     &
         (k**2*mf2d10rho**2*nfad2x)/(4.*(1 + mf2d00rho)) +                      &
         (k*mf2d10rho**2*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -                    &
         (k**2*mf2d10rho**2*nffd2x)/(4.*(1 + mf2d00rho))))/                     &
     (4.*(1 + mf2d00rho)**3.5) +                                                &
    (5*(1 - etapsi/4.)*mf2d10rho**2*                                            &
       ((-3*k*mf2d10rho**3*nfad1x)/(8.*(1 + mf2d00rho)**2.5) +                  &
         (3*k**2*mf2d10rho**3*nfad2x)/(8.*(1 + mf2d00rho)**2) -                 &
         (k**3*mf2d10rho**3*nfad3x)/(8.*(1 + mf2d00rho)**1.5) -                 &
         (3*k*mf2d10rho**3*nffd1x)/(8.*(1 + mf2d00rho)**2.5) +                  &
         (3*k**2*mf2d10rho**3*nffd2x)/(8.*(1 + mf2d00rho)**2) -                 &
         (k**3*mf2d10rho**3*nffd3x)/(8.*(1 + mf2d00rho)**1.5)))/                &
     (2.*(1 + mf2d00rho)**2.5) -                                                &
    (5*(1 - etapsi/4.)*mf2d10rho*                                               &
       ((15*k*mf2d10rho**4*nfad1x)/(16.*(1 + mf2d00rho)**3.5) -                 &
         (15*k**2*mf2d10rho**4*nfad2x)/(16.*(1 + mf2d00rho)**3) +               &
         (3*k**3*mf2d10rho**4*nfad3x)/(8.*(1 + mf2d00rho)**2.5) -               &
         (k**4*mf2d10rho**4*nfad4x)/(16.*(1 + mf2d00rho)**2) +                  &
         (15*k*mf2d10rho**4*nffd1x)/(16.*(1 + mf2d00rho)**3.5) -                &
         (15*k**2*mf2d10rho**4*nffd2x)/(16.*(1 + mf2d00rho)**3) +               &
         (3*k**3*mf2d10rho**4*nffd3x)/(8.*(1 + mf2d00rho)**2.5) -               &
         (k**4*mf2d10rho**4*nffd4x)/(16.*(1 + mf2d00rho)**2)))/                 &
     (6.*(1 + mf2d00rho)**1.5) +                                                &
    ((1 - etapsi/4.)*((-105*k*mf2d10rho**5*nfad1x)/                             &
          (32.*(1 + mf2d00rho)**4.5) +                                          &
         (105*k**2*mf2d10rho**5*nfad2x)/(32.*(1 + mf2d00rho)**4) -              &
         (45*k**3*mf2d10rho**5*nfad3x)/(32.*(1 + mf2d00rho)**3.5) +             &
         (5*k**4*mf2d10rho**5*nfad4x)/(16.*(1 + mf2d00rho)**3) -                &
         (k**5*mf2d10rho**5*nfad5x)/(32.*(1 + mf2d00rho)**2.5) -                &
         (105*k*mf2d10rho**5*nffd1x)/(32.*(1 + mf2d00rho)**4.5) +               &
         (105*k**2*mf2d10rho**5*nffd2x)/(32.*(1 + mf2d00rho)**4) -              &
         (45*k**3*mf2d10rho**5*nffd3x)/(32.*(1 + mf2d00rho)**3.5) +             &
         (5*k**4*mf2d10rho**5*nffd4x)/(16.*(1 + mf2d00rho)**3) -                &
         (k**5*mf2d10rho**5*nffd5x)/(32.*(1 + mf2d00rho)**2.5)))/               &
     (3.*Sqrt(1 + mf2d00rho))
  lft01= -((1 - etapsi/4.)*mf2d01rho*(1 - nfad0x - nffd0x))/                    &
     (6.*(1 + mf2d00rho)**1.5) +                                                &
         ((1 - etapsi/4.)*(-(k*mf2d01rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -     &
                  (k*mf2d01rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/              &
                       (3.*Sqrt(1 + mf2d00rho))
  lft11= ((1 - etapsi/4.)*mf2d01rho*mf2d10rho*(1 - nfad0x - nffd0x))/           &
     (4.*(1 + mf2d00rho)**2.5) -                                                &
    ((1 - etapsi/4.)*mf2d10rho*                                                 &
       (-(k*mf2d01rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d01rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (6.*(1 + mf2d00rho)**1.5) -                                                &
    ((1 - etapsi/4.)*mf2d01rho*                                                 &
       (-(k*mf2d10rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d10rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (6.*(1 + mf2d00rho)**1.5) +                                                &
    ((1 - etapsi/4.)*((k*mf2d01rho*mf2d10rho*nfad1x)/                           &
          (4.*(1 + mf2d00rho)**1.5) -                                           &
         (k**2*mf2d01rho*mf2d10rho*nfad2x)/(4.*(1 + mf2d00rho)) +               &
         (k*mf2d01rho*mf2d10rho*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -             &
         (k**2*mf2d01rho*mf2d10rho*nffd2x)/(4.*(1 + mf2d00rho))))/              &
     (3.*Sqrt(1 + mf2d00rho))
  lft21= (-5*(1 - etapsi/4.)*mf2d01rho*mf2d10rho**2*(1 - nfad0x - nffd0x))/     &
     (8.*(1 + mf2d00rho)**3.5) +                                                &
    ((1 - etapsi/4.)*mf2d10rho**2*                                              &
       (-(k*mf2d01rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d01rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (4.*(1 + mf2d00rho)**2.5) +                                                &
    ((1 - etapsi/4.)*mf2d01rho*mf2d10rho*                                       &
       (-(k*mf2d10rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d10rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (2.*(1 + mf2d00rho)**2.5) -                                                &
    ((1 - etapsi/4.)*mf2d10rho*                                                 &
       ((k*mf2d01rho*mf2d10rho*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -              &
         (k**2*mf2d01rho*mf2d10rho*nfad2x)/(4.*(1 + mf2d00rho)) +               &
         (k*mf2d01rho*mf2d10rho*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -             &
         (k**2*mf2d01rho*mf2d10rho*nffd2x)/(4.*(1 + mf2d00rho))))/              &
     (3.*(1 + mf2d00rho)**1.5) -                                                &
    ((1 - etapsi/4.)*mf2d01rho*                                                 &
       ((k*mf2d10rho**2*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -                     &
         (k**2*mf2d10rho**2*nfad2x)/(4.*(1 + mf2d00rho)) +                      &
         (k*mf2d10rho**2*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -                    &
         (k**2*mf2d10rho**2*nffd2x)/(4.*(1 + mf2d00rho))))/                     &
     (6.*(1 + mf2d00rho)**1.5) +                                                &
    ((1 - etapsi/4.)*((-3*k*mf2d01rho*mf2d10rho**2*nfad1x)/                     &
          (8.*(1 + mf2d00rho)**2.5) +                                           &
         (3*k**2*mf2d01rho*mf2d10rho**2*nfad2x)/(8.*(1 + mf2d00rho)**2) -       &
         (k**3*mf2d01rho*mf2d10rho**2*nfad3x)/(8.*(1 + mf2d00rho)**1.5) -       &
         (3*k*mf2d01rho*mf2d10rho**2*nffd1x)/(8.*(1 + mf2d00rho)**2.5) +        &
         (3*k**2*mf2d01rho*mf2d10rho**2*nffd2x)/(8.*(1 + mf2d00rho)**2) -       &
         (k**3*mf2d01rho*mf2d10rho**2*nffd3x)/(8.*(1 + mf2d00rho)**1.5)))/      &
     (3.*Sqrt(1 + mf2d00rho))
  lft31= (35*(1 - etapsi/4.)*mf2d01rho*mf2d10rho**3*(1 - nfad0x - nffd0x))/     &
     (16.*(1 + mf2d00rho)**4.5) -                                               &
    (5*(1 - etapsi/4.)*mf2d10rho**3*                                            &
       (-(k*mf2d01rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d01rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (8.*(1 + mf2d00rho)**3.5) -                                                &
    (15*(1 - etapsi/4.)*mf2d01rho*mf2d10rho**2*                                 &
       (-(k*mf2d10rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d10rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (8.*(1 + mf2d00rho)**3.5) +                                                &
    (3*(1 - etapsi/4.)*mf2d10rho**2*                                            &
       ((k*mf2d01rho*mf2d10rho*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -              &
         (k**2*mf2d01rho*mf2d10rho*nfad2x)/(4.*(1 + mf2d00rho)) +               &
         (k*mf2d01rho*mf2d10rho*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -             &
         (k**2*mf2d01rho*mf2d10rho*nffd2x)/(4.*(1 + mf2d00rho))))/              &
     (4.*(1 + mf2d00rho)**2.5) +                                                &
    (3*(1 - etapsi/4.)*mf2d01rho*mf2d10rho*                                     &
       ((k*mf2d10rho**2*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -                     &
         (k**2*mf2d10rho**2*nfad2x)/(4.*(1 + mf2d00rho)) +                      &
         (k*mf2d10rho**2*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -                    &
         (k**2*mf2d10rho**2*nffd2x)/(4.*(1 + mf2d00rho))))/                     &
     (4.*(1 + mf2d00rho)**2.5) -                                                &
    ((1 - etapsi/4.)*mf2d10rho*                                                 &
       ((-3*k*mf2d01rho*mf2d10rho**2*nfad1x)/(8.*(1 + mf2d00rho)**2.5) +        &
         (3*k**2*mf2d01rho*mf2d10rho**2*nfad2x)/(8.*(1 + mf2d00rho)**2) -       &
         (k**3*mf2d01rho*mf2d10rho**2*nfad3x)/(8.*(1 + mf2d00rho)**1.5) -       &
         (3*k*mf2d01rho*mf2d10rho**2*nffd1x)/(8.*(1 + mf2d00rho)**2.5) +        &
         (3*k**2*mf2d01rho*mf2d10rho**2*nffd2x)/(8.*(1 + mf2d00rho)**2) -       &
         (k**3*mf2d01rho*mf2d10rho**2*nffd3x)/(8.*(1 + mf2d00rho)**1.5)))/      &
     (2.*(1 + mf2d00rho)**1.5) -                                                &
    ((1 - etapsi/4.)*mf2d01rho*                                                 &
       ((-3*k*mf2d10rho**3*nfad1x)/(8.*(1 + mf2d00rho)**2.5) +                  &
         (3*k**2*mf2d10rho**3*nfad2x)/(8.*(1 + mf2d00rho)**2) -                 &
         (k**3*mf2d10rho**3*nfad3x)/(8.*(1 + mf2d00rho)**1.5) -                 &
         (3*k*mf2d10rho**3*nffd1x)/(8.*(1 + mf2d00rho)**2.5) +                  &
         (3*k**2*mf2d10rho**3*nffd2x)/(8.*(1 + mf2d00rho)**2) -                 &
         (k**3*mf2d10rho**3*nffd3x)/(8.*(1 + mf2d00rho)**1.5)))/                &
     (6.*(1 + mf2d00rho)**1.5) +                                                &
    ((1 - etapsi/4.)*((15*k*mf2d01rho*mf2d10rho**3*nfad1x)/                     &
          (16.*(1 + mf2d00rho)**3.5) -                                          &
         (15*k**2*mf2d01rho*mf2d10rho**3*nfad2x)/                               &
          (16.*(1 + mf2d00rho)**3) +                                            &
         (3*k**3*mf2d01rho*mf2d10rho**3*nfad3x)/                                &
          (8.*(1 + mf2d00rho)**2.5) -                                           &
         (k**4*mf2d01rho*mf2d10rho**3*nfad4x)/(16.*(1 + mf2d00rho)**2) +        &
         (15*k*mf2d01rho*mf2d10rho**3*nffd1x)/(16.*(1 + mf2d00rho)**3.5) -      &
         (15*k**2*mf2d01rho*mf2d10rho**3*nffd2x)/                               &
          (16.*(1 + mf2d00rho)**3) +                                            &
         (3*k**3*mf2d01rho*mf2d10rho**3*nffd3x)/                                &
          (8.*(1 + mf2d00rho)**2.5) -                                           &
         (k**4*mf2d01rho*mf2d10rho**3*nffd4x)/(16.*(1 + mf2d00rho)**2)))/       &
     (3.*Sqrt(1 + mf2d00rho))
  lft02= ((1 - etapsi/4.)*mf2d01rho**2*(1 - nfad0x - nffd0x))/                  &
     (4.*(1 + mf2d00rho)**2.5) -                                                &
    ((1 - etapsi/4.)*mf2d02rho*(1 - nfad0x - nffd0x))/                          &
     (6.*(1 + mf2d00rho)**1.5) -                                                &
    ((1 - etapsi/4.)*mf2d01rho*                                                 &
       (-(k*mf2d01rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d01rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (3.*(1 + mf2d00rho)**1.5) +                                                &
    ((1 - etapsi/4.)*((k*mf2d01rho**2*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -       &
         (k*mf2d02rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k**2*mf2d01rho**2*nfad2x)/(4.*(1 + mf2d00rho)) +                      &
         (k*mf2d01rho**2*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -                    &
         (k*mf2d02rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k**2*mf2d01rho**2*nffd2x)/(4.*(1 + mf2d00rho))))/                     &
     (3.*Sqrt(1 + mf2d00rho))
  lft12= (-5*(1 - etapsi/4.)*mf2d01rho**2*mf2d10rho*(1 - nfad0x - nffd0x))/     &
     (8.*(1 + mf2d00rho)**3.5) +                                                &
    ((1 - etapsi/4.)*mf2d02rho*mf2d10rho*(1 - nfad0x - nffd0x))/                &
     (4.*(1 + mf2d00rho)**2.5) +                                                &
    ((1 - etapsi/4.)*mf2d01rho*mf2d10rho*                                       &
       (-(k*mf2d01rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d01rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (2.*(1 + mf2d00rho)**2.5) +                                                &
    ((1 - etapsi/4.)*mf2d01rho**2*                                              &
       (-(k*mf2d10rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d10rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (4.*(1 + mf2d00rho)**2.5) -                                                &
    ((1 - etapsi/4.)*mf2d02rho*                                                 &
       (-(k*mf2d10rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k*mf2d10rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho))))/                       &
     (6.*(1 + mf2d00rho)**1.5) -                                                &
    ((1 - etapsi/4.)*mf2d10rho*                                                 &
       ((k*mf2d01rho**2*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -                     &
         (k*mf2d02rho*nfad1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k**2*mf2d01rho**2*nfad2x)/(4.*(1 + mf2d00rho)) +                      &
         (k*mf2d01rho**2*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -                    &
         (k*mf2d02rho*nffd1x)/(2.*Sqrt(1 + mf2d00rho)) -                        &
         (k**2*mf2d01rho**2*nffd2x)/(4.*(1 + mf2d00rho))))/                     &
     (6.*(1 + mf2d00rho)**1.5) -                                                &
    ((1 - etapsi/4.)*mf2d01rho*                                                 &
       ((k*mf2d01rho*mf2d10rho*nfad1x)/(4.*(1 + mf2d00rho)**1.5) -              &
         (k**2*mf2d01rho*mf2d10rho*nfad2x)/(4.*(1 + mf2d00rho)) +               &
         (k*mf2d01rho*mf2d10rho*nffd1x)/(4.*(1 + mf2d00rho)**1.5) -             &
         (k**2*mf2d01rho*mf2d10rho*nffd2x)/(4.*(1 + mf2d00rho))))/              &
     (3.*(1 + mf2d00rho)**1.5) +                                                &
    ((1 - etapsi/4.)*((-3*k*mf2d01rho**2*mf2d10rho*nfad1x)/                     &
          (8.*(1 + mf2d00rho)**2.5) +                                           &
         (k*mf2d02rho*mf2d10rho*nfad1x)/(4.*(1 + mf2d00rho)**1.5) +             &
         (3*k**2*mf2d01rho**2*mf2d10rho*nfad2x)/(8.*(1 + mf2d00rho)**2) -       &
         (k**2*mf2d02rho*mf2d10rho*nfad2x)/(4.*(1 + mf2d00rho)) -               &
         (k**3*mf2d01rho**2*mf2d10rho*nfad3x)/(8.*(1 + mf2d00rho)**1.5) -       &
         (3*k*mf2d01rho**2*mf2d10rho*nffd1x)/(8.*(1 + mf2d00rho)**2.5) +        &
         (k*mf2d02rho*mf2d10rho*nffd1x)/(4.*(1 + mf2d00rho)**1.5) +             &
         (3*k**2*mf2d01rho**2*mf2d10rho*nffd2x)/(8.*(1 + mf2d00rho)**2) -       &
         (k**2*mf2d02rho*mf2d10rho*nffd2x)/(4.*(1 + mf2d00rho)) -               &
         (k**3*mf2d01rho**2*mf2d10rho*nffd3x)/(8.*(1 + mf2d00rho)**1.5)))/      &
     (3.*Sqrt(1 + mf2d00rho)) 

  lft(:,1)=lft00
  lft(:,2)=lft10
  lft(:,3)=lft20
  lft(:,4)=lft30
  lft(:,5)=lft40
  lft(:,6)=lft50
  lft(:,7)=lft01
  lft(:,8)=lft11
  lft(:,9)=lft21
  lft(:,10)=lft31
  lft(:,11)=lft02
  lft(:,12)=lft12
end 
