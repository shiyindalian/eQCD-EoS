subroutine IRenha(k,ashift,dashift)
!this subroutine gives the enhanced factor for the IR
  implicit none

  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) k,ashift,dashift
  real(16) IRstrength,cut
  real(16) kGeV
  integer n

  kGeV=k*hc*1.Q-3
  IRstrength=0.034Q+0
  cut=2.Q+0            
  !GeV
  n=2

  if((kGeV/cut)**n>150.Q+0)then
    ashift=1.Q+0
    dashift=0.Q+0
  else
    ashift=1.Q+0 + (IRstrength*(kGeV/cut)**n)/(-1. + exp((kGeV/cut)**n))
    dashift=-((IRstrength*(kGeV/cut)**n*(1. + exp((kGeV/cut)**n)*(-1. + (kGeV/cut)**n))*n)/ &
    (-1. + exp((kGeV/cut)**n))**2)
  end if

end


subroutine massAScreen(k,T,dmass2S,ptdmass2S)
!this subroutine gives the enhanced factor for the IR
  implicit none

  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) k,T,dmass2S,ptdmass2S
  real(16) cs,kT
  integer ns

  kT=2.*pi*T*(1.Q+0/2.Q+0)
  cs=2.Q+0

  ns=2
  dmass2S=cs**2*T**2*exp(-(k/kT)**ns)
  ptdmass2S=-((cs**2*k*(k/kT)**(-1 + ns)*ns*T**2*exp(-(k/kT)**ns))/kT)
end
