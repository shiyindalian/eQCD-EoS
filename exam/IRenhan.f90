subroutine IRenha(k,ashift,dashift)
!this subroutine gives the enhanced factor for the IR
  implicit none

  real(8) pi,hc
  parameter(pi=3.141592653589793d0)
  parameter(hc=197.33)

  real(8) k,ashift,dashift
  real(8) IRstrength,cut
  real(8) kGeV
  integer n

  kGeV=k*hc*1.d-3
!  IRstrength=0.1515
!  IRstrength=0.160
!   IRstrength=0.0375
!  IRstrength=0.07
!  IRstrength=0.068
  IRstrength=0.034


!  cut=1.3            !GeV
!  n=3

!  cut=5.            !GeV
  cut=2.            !GeV
  n=2

  if((kGeV/cut)**n>150.)then
    ashift=1.
    dashift=0.
  else
    ashift=1. + (IRstrength*(kGeV/cut)**n)/(-1. + exp((kGeV/cut)**n))
    dashift=-((IRstrength*(kGeV/cut)**n*(1. + exp((kGeV/cut)**n)*(-1. + (kGeV/cut)**n))*n)/ &
    (-1. + exp((kGeV/cut)**n))**2)
  end if

end


subroutine massAScreen(k,T,dmass2S,ptdmass2S)
!this subroutine gives the enhanced factor for the IR
  implicit none

  real(8) pi,hc
  parameter(pi=3.141592653589793d0)
  parameter(hc=197.33)

  real(8) k,T,dmass2S,ptdmass2S

  real(8) cs,kT
  integer ns


  kT=2.*pi*T*(1./2.)
  cs=2.

  ns=2

  dmass2S=cs**2*T**2*exp(-(k/kT)**ns)

  ptdmass2S=-((cs**2*k*(k/kT)**(-1 + ns)*ns*T**2*exp(-(k/kT)**ns))/kT)


end




