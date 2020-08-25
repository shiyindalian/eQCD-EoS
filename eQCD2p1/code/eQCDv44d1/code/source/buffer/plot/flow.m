%figure(1)
%plot(T,fpi,'k-','LineWidth',2.);hold on;
%plot(T1,fpi1,'r-.','LineWidth',2);hold on;
%plot(T2,fpi2,'m-. ','LineWidth',2);hold on;
%plot(T3,fpi3,'b-. ','LineWidth',2);hold on;
%plot(T4,fpi4,'g-. ','LineWidth',2);hold on;
%xlabel('\fontsize{18} x');
%ylabel('\fontsize{18} y');
figure(1)
loglog(k_pk,mpion_pk,'k-','LineWidth',2);hold on;
loglog(k_p0,mpion_p0,'b-.','LineWidth',2);hold on;
loglog(k_pk,mpion,'g-.','LineWidth',2);hold on;
%semilogx(k4,fpi4,'m-.','LineWidth',2);hold on;
%semilogx(k_new1,dtg_new1,'r.','LineWidth',2);hold on;
%semilogx(k_new1,dtg3A_new1,'g.','LineWidth',2);hold on;
%semilogx(k_new1,etaAQLsT0_new1,'b-.','LineWidth',2);hold on;
%plot(T,fpi450,'b-. ','LineWidth',2);hold on;
%plot(T,fpi470,'m-. ','LineWidth',2);hold on;
xlim([1,25000])
xlabel('\fontsize{18} x');
ylabel('\fontsize{18} y');

figure(2)
semilogx(k_T150_p0_mp238,mf_T150_p0_mp238,'k-','LineWidth',2);hold on;
semilogx(k_T150_pk,mf_T150_pk,'b-.','LineWidth',2);hold on;
%semilogx(k3,fpi3,'g-.','LineWidth',2);hold on;
%semilogx(k4,fpi4,'m-.','LineWidth',2);hold on;
%semilogx(k_new1,dtg_new1,'r.','LineWidth',2);hold on;
%semilogx(k_new1,dtg3A_new1,'g.','LineWidth',2);hold on;
%semilogx(k_new1,etaAQLsT0_new1,'b-.','LineWidth',2);hold on;
%plot(T,fpi450,'b-. ','LineWidth',2);hold on;
%plot(T,fpi470,'m-. ','LineWidth',2);hold on;
xlim([1,25000])
xlabel('\fontsize{18} x');
ylabel('\fontsize{18} y');