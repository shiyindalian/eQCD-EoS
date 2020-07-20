program QM

  implicit none

  real(16) pi
  parameter(pi=3.141592653589793238462643383279Q+00)
  integer i,j,jl,Tnum,munum,muup
  parameter(Tnum=250,munum=99)
  real(16) dt
  real(16) T,mu
  !temperature and chemical potential
  real(16) mu_up,mu_down
  real(16) l,lb 
  !polyakov loop
  real(16) l_i1,l_i2,lb_i1,lb_i2
  !polyakov loop
  real(16) Vtotal0
  real(16) sigma_UV_i(2),sigma_UV(2)
  real(16) k_UV_i(2,munum),k_UV(2,munum),k_UV_0(2)
  real(16) V_res(Tnum,munum),T_res(Tnum,munum),P_res(Tnum,munum),mu_res(Tnum,munum)
  real(16) pre_com(munum)
  integer  mmax
  parameter(mmax=10)
  !order of chebyshev polynomial
  real(16) dcd0mu(mmax),dcd1mu(mmax),dcd2mu(mmax),dcd3mu(mmax),dcd4mu(mmax),dcd5mu(mmax),dcd6mu(mmax)
  real(16) chi(mmax,Tnum)
  real(16) factorial
  real(16) chebev
  external chebev

  character(14) fileV1
  character(15) fileV2

  common /Tmu/ T,mu
  common /prefit/ pre_com

  dT=1.Q+00/197.33Q+0

  mu_up=250./197.33Q+0
  mu_down=-mu_up

  do jl=1,9
    j=jl
    write(fileV1,'(a9,i1,a4)')'../data/V',jl,'.dat'
    open(unit=51,file=Trim(fileV1))
    do i=1,Tnum
      read(51,*)V_res(i,j)
    end do
    close(51)
  end do

  do jl=10,munum
    j=jl
    write(fileV2,'(a9,i2,a4)')'../data/V',jl,'.dat'
    open(unit=51,file=Trim(fileV2))
    do i=1,Tnum
      read(51,*)V_res(i,j)
    end do
    close(51)
  end do

  Vtotal0=V_res(1,(munum+1)/2)

  P_res=-V_res+Vtotal0

  do i=1,Tnum
    T=real(i,kind=16)*dt
    pre_com=P_res(i,:)
    call chebft(dcd0mu,munum,mmax)
    call chder(mu_down,mu_up,dcd0mu,dcd1mu,mmax)
    call chder(mu_down,mu_up,dcd1mu,dcd2mu,mmax)
    call chder(mu_down,mu_up,dcd2mu,dcd3mu,mmax)
    call chder(mu_down,mu_up,dcd3mu,dcd4mu,mmax)
    call chder(mu_down,mu_up,dcd4mu,dcd5mu,mmax)
    call chder(mu_down,mu_up,dcd5mu,dcd6mu,mmax)

    chi(1,i)=chebev(mu_down,mu_up,dcd0mu,mmax,0.Q+00)/T**4
    chi(2,i)=chebev(mu_down,mu_up,dcd1mu,mmax,0.Q+00)/T**3
    chi(3,i)=chebev(mu_down,mu_up,dcd2mu,mmax,0.Q+00)/T**2
    chi(4,i)=chebev(mu_down,mu_up,dcd3mu,mmax,0.Q+00)/T
    chi(5,i)=chebev(mu_down,mu_up,dcd4mu,mmax,0.Q+00)
    chi(6,i)=chebev(mu_down,mu_up,dcd5mu,mmax,0.Q+00)*T
    chi(7,i)=chebev(mu_down,mu_up,dcd6mu,mmax,0.Q+00)*T**2
  end do
  open(unit=51,file='../data/chi0.dat')
  do i=1,Tnum
    write(51, *)chi(1,i)
  end do
  close(51)

  open(unit=51,file='../data/chi1.dat')
  do i=1,Tnum
    write(51, *)chi(2,i)
  end do
  close(51)

  open(unit=51,file='../data/chi2.dat')
  do i=1,Tnum
    write(51, *)chi(3,i)
  end do
  close(51)

  open(unit=51,file='../data/chi3.dat')
  do i=1,Tnum
    write(51, *)chi(4,i)
  end do
  close(51)

  open(unit=51,file='../data/chi4.dat')
  do i=1,Tnum
    write(51, *)chi(5,i)
  end do
  close(51)

  open(unit=51,file='../data/chi5.dat')
  do i=1,Tnum
    write(51, *)chi(6,i)
  end do
  close(51)

  open(unit=51,file='../data/chi6.dat')
  do i=1,Tnum
    write(51, *)chi(7,i)
  end do
  close(51)

  open(unit=51,file='../data/R42.dat')
  do i=1,Tnum
    write(51, *)chi(5,i)/chi(3,i)
  end do
  close(51)

end
