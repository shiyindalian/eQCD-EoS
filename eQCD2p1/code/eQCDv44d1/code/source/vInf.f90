subroutine vInf(k_UV,Vinfi)
!calculate the contribution to thermaldynamical potential arising from
!momenta above UV cutoff

  implicit none
  real(8) k_UV,Vinfi
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) T,mu
  real(8) l_com,lb_com
  real(8) l,lb,l_con,lb_con
  integer npoint, in
  parameter(npoint=256)
  real(8) w(npoint),y(npoint) !Guass integral
  real(8) k,k_infi
  real(8) Nc,Nf
  parameter(Nc=3.,Nf=2.)
  real(8) v3
  parameter(v3=1./(2.*pi**2))
  real(8) Fnf0,Fnf1,Fnf2
  external Fnf0,Fnf1,Fnf2
  real(8) nff,nfa
  real(8) nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x
  real(8) xff,xfa



  common /Tmu/ T,mu
  common /polyakov_com/ l_com,lb_com


  l=l_com
  lb=lb_com

  l_con=l
  lb_con=lb


  k_infi=k_UV*2.

  call gauleg(k_UV, k_infi, y, w, npoint)

  Vinfi=0.
  do in=1, npoint
    k=y(in)

    xff=-mu + k
    xfa=mu + k

    l=l_con
    lb=lb_con

    nf0=Fnf0(xff,T,l,lb)
    nf1=Fnf1(xff,T,l,lb)
    nf2=Fnf2(xff,T,l,lb)

    nff=nf0 + lb*nf1 + l*nf2

    l=lb_con
    lb=l_con

    nf0=Fnf0(xfa,T,l,lb)
    nf1=Fnf1(xfa,T,l,lb)
    nf2=Fnf2(xfa,T,l,lb)

    nfa=nf0 + lb*nf1 + l*nf2

    l=l_con
    lb=lb_con


    Vinfi=Vinfi+w(in)*k**3*(0.-nfa-nff)
  end do
  Vinfi=Vinfi*v3*4*Nc*Nf/3./2.

end
