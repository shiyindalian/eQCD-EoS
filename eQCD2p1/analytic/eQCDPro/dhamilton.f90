  hd00p22=lam10 + (ck*Sqrt(rho + Sqrt(6)*Sqrt(rho2)))/Sqrt(3) - 
    Sqrt(0.6666666666666666)*lam01*Sqrt(rho2),
  hd10p22=lam20 + (ck/Sqrt(rho + Sqrt(6)*Sqrt(rho2)) - 
       2*Sqrt(2)*lam11*Sqrt(rho2))/(2.*Sqrt(3)),
  hd20p22=lam30 - ck/(4.*Sqrt(3)*(rho + Sqrt(6)*Sqrt(rho2))**1.5) - 
    Sqrt(0.6666666666666666)*lam21*Sqrt(rho2),
  hd30p22=lam40 + (Sqrt(3)*ck)/(8.*(rho + Sqrt(6)*Sqrt(rho2))**2.5) - 
    Sqrt(0.6666666666666666)*lam31*Sqrt(rho2),
  hd40p22=lam50 - (5*Sqrt(3)*ck)/(16.*(rho + Sqrt(6)*Sqrt(rho2))**3.5),
  hd50p22=(35*Sqrt(3)*ck)/(32.*(rho + Sqrt(6)*Sqrt(rho2))**4.5),
  hd01p22=lam11 + (-2*Sqrt(3)*lam01 + 
       (3*ck)/Sqrt(rho + Sqrt(6)*Sqrt(rho2)) - 4*Sqrt(3)*lam02*rho2)/
     (6.*Sqrt(2)*Sqrt(rho2)),hd11p22=
   lam21 + (-4*Sqrt(3)*lam11 - (3*ck)/(rho + Sqrt(6)*Sqrt(rho2))**1.5 - 
       8*Sqrt(3)*lam12*rho2)/(12.*Sqrt(2)*Sqrt(rho2)),
  hd21p22=lam31 + (-8*Sqrt(3)*lam21 + 
       (9*ck)/(rho + Sqrt(6)*Sqrt(rho2))**2.5)/(24.*Sqrt(2)*Sqrt(rho2)),
  hd31p22=(-16*Sqrt(3)*lam31 - (45*ck)/(rho + Sqrt(6)*Sqrt(rho2))**3.5)/
    (48.*Sqrt(2)*Sqrt(rho2)),hd02p22=
   (24*lam12 + (-3*ck*(Sqrt(2)*rho + 3*Sqrt(3)*Sqrt(rho2)) + 
         2*(Sqrt(6)*rho + 6*Sqrt(rho2))*Sqrt(rho + Sqrt(6)*Sqrt(rho2))*
          (lam01 - 4*lam02*rho2))/
       ((rho + Sqrt(6)*Sqrt(rho2))**1.5*rho2**1.5))/24.,
  hd12p22=(4*Sqrt(6)*lam11 + 
      (3*ck*(Sqrt(2)*rho + 5*Sqrt(3)*Sqrt(rho2)))/
       (rho + Sqrt(6)*Sqrt(rho2))**2.5 - 16*Sqrt(6)*lam12*rho2)/
    (48.*rho2**1.5),hd00s19=
   (9*lam10 + 12*lam20*rho + lam01*(4*rho - 5*Sqrt(6)*Sqrt(rho2)) - 
      3*Sqrt(3)*ck*Sqrt(rho + Sqrt(6)*Sqrt(rho2)) - 
      6*Sqrt(6)*lam20*Sqrt(rho2) - 8*Sqrt(6)*lam11*rho*Sqrt(rho2) + 
      24*lam11*rho2 + lam02*(8*rho*rho2 - 4*Sqrt(6)*rho2**1.5))/9.,
  hd10s19=(4*lam01 + 21*lam20 + 4*lam11*rho + 12*lam30*rho - 
      (3*Sqrt(3)*ck)/(2.*Sqrt(rho + Sqrt(6)*Sqrt(rho2))) - 
      13*Sqrt(6)*lam11*Sqrt(rho2) - 6*Sqrt(6)*lam30*Sqrt(rho2) - 
      8*Sqrt(6)*lam21*rho*Sqrt(rho2) + 8*lam02*rho2 + 24*lam21*rho2 + 
      8*lam12*rho*rho2 - 4*Sqrt(6)*lam12*rho2**1.5)/9.,
  hd20s19=(32*lam11 + 132*lam30 + 16*lam21*rho + 48*lam40*rho + 
      (3*Sqrt(3)*ck)/(rho + Sqrt(6)*Sqrt(rho2))**1.5 - 
      84*Sqrt(6)*lam21*Sqrt(rho2) - 24*Sqrt(6)*lam40*Sqrt(rho2) - 
      32*Sqrt(6)*lam31*rho*Sqrt(rho2) + 64*lam12*rho2 + 96*lam31*rho2)/36.,
  hd30s19=(96*lam21 + 360*lam40 + 32*lam31*rho + 96*lam50*rho - 
      (9*Sqrt(3)*ck)/(rho + Sqrt(6)*Sqrt(rho2))**2.5 - 
      232*Sqrt(6)*lam31*Sqrt(rho2) - 48*Sqrt(6)*lam50*Sqrt(rho2))/72.,
  hd40s19=(16*lam31)/9. + (19*lam50)/3. + 
    (5*Sqrt(3)*ck)/(16.*(rho + Sqrt(6)*Sqrt(rho2))**3.5),
  hd50s19=(-35*Sqrt(3)*ck)/(32.*(rho + Sqrt(6)*Sqrt(rho2))**4.5),
  hd01s19=(-9*Sqrt(2)*ck - 2*Sqrt(rho + Sqrt(6)*Sqrt(rho2))*
       (5*Sqrt(6)*lam01 + 6*Sqrt(6)*lam20 + 8*Sqrt(6)*lam11*rho - 
         66*lam11*Sqrt(rho2) - 24*lam02*rho*Sqrt(rho2) - 
         24*lam21*rho*Sqrt(rho2) + 22*Sqrt(6)*lam02*rho2 + 
         12*Sqrt(6)*lam21*rho2 + 16*Sqrt(6)*lam12*rho*rho2 - 
         48*lam12*rho2**1.5))/
    (36.*Sqrt(rho + Sqrt(6)*Sqrt(rho2))*Sqrt(rho2)),
  hd11s19=(96*lam02 + 96*lam12*rho + 96*lam31*rho + 
      lam21*(360 - (32*Sqrt(6)*rho)/Sqrt(rho2)) - 
      (52*Sqrt(6)*lam11)/Sqrt(rho2) - (24*Sqrt(6)*lam30)/Sqrt(rho2) + 
      (9*Sqrt(2)*ck)/((rho + Sqrt(6)*Sqrt(rho2))**1.5*Sqrt(rho2)) - 
      152*Sqrt(6)*lam12*Sqrt(rho2) - 48*Sqrt(6)*lam31*Sqrt(rho2))/72.,
  hd21s19=(8*lam12)/3. + (lam31*(57 - (4*Sqrt(6)*rho)/Sqrt(rho2)))/9. + 
    (-56*Sqrt(3)*lam21 - 16*Sqrt(3)*lam40 - 
       (9*ck)/(rho + Sqrt(6)*Sqrt(rho2))**2.5)/(24.*Sqrt(2)*Sqrt(rho2)),
  hd31s19=(-464*Sqrt(3)*lam31 - 96*Sqrt(3)*lam50 + 
      (135*ck)/(rho + Sqrt(6)*Sqrt(rho2))**3.5)/(144.*Sqrt(2)*Sqrt(rho2)),
  hd02s19=(lam12*(57 - (8*Sqrt(6)*rho)/Sqrt(rho2)) + 
      (10*Sqrt(6)*lam01 + 12*Sqrt(6)*lam20 + 16*Sqrt(6)*lam11*rho + 
         (9*Sqrt(2)*ck)/Sqrt(rho + Sqrt(6)*Sqrt(rho2)) + 
         (9*Sqrt(3)*ck*Sqrt(rho2))/(rho + Sqrt(6)*Sqrt(rho2))**1.5 - 
         64*Sqrt(6)*lam02*rho2 - 48*Sqrt(6)*lam21*rho2)/(8.*rho2**1.5))/9.,
  hd12s19=(52*Sqrt(6)*lam11 + 24*Sqrt(6)*lam30 + 32*Sqrt(6)*lam21*rho - 
      (9*Sqrt(2)*ck)/(rho + Sqrt(6)*Sqrt(rho2))**1.5 - 
      (27*Sqrt(3)*ck*Sqrt(rho2))/(rho + Sqrt(6)*Sqrt(rho2))**2.5 - 
      256*Sqrt(6)*lam12*rho2 - 96*Sqrt(6)*lam31*rho2)/(144.*rho2**1.5)
