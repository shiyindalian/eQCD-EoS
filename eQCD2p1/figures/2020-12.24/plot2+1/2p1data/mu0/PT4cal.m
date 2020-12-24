clear
hc=197.33;
load('VTOTAL.DAT');
T=load('TMEV.DAT');
p=-VTOTAL+VTOTAL(1);
pT4=p./(T/hc).^4;
plot(T,pT4,'g')
trace=(pT4(2:end)-pT4(1:end-1))/(2/hc);
hold on
%plot(T(2:end),trace)
axis([0,200,-0.5,3.5])