clc
clear all
close all

%length=0.15;
%freq=10^9;
%k=2*pi*freq/(3*(10^8));
%kl=k*length;
%k_l=k/length;


C=0.577215664901532860606512090082402431042; % Euler constant
eta=377;
a=10^-3; % Dipoles diameter


kl=linspace(eps,2*pi*3,50);
imp_real=(eta/(2*pi)).*(C+log(kl)-cosint(kl) + 0.5*sin(kl).*(sinint(2*kl)- sinint(kl)*2) + ...
     0.5*cos(kl).*(C +log(0.5*kl) +cosint(2*kl)- cosint(kl)*2))./(sin(0.5*kl).^2); 
 %figure(1);plot(kl,imp_real);axis([0 6*pi 0 1000])
imp_imag1=(eta/(4*pi)).*(2*sinint(kl)+ cos(kl).*(2*sinint(kl) - sinint(2*kl)) - sin(kl).*(2*cosint(kl)-cosint(2*kl)))./(sin(0.5*kl).^2);
%figure(2);plot(kl,imp_imag1);axis([0 6*pi -1000 1000])


impreal_fit = fit(kl(~excludedata(kl,imp_real,'range',[0 2000]))',imp_real(~excludedata(kl,imp_real,'range',[0 2000]))','smoothingspline')
%figure(3);plot(linspace(eps,6*pi,5000),impreal_fit(linspace(eps,6*pi,5000)));axis([0 6*pi 0 1000])
impimag1_fit = fit(kl(~excludedata(kl,imp_imag1,'range',[0 2000]))',imp_imag1(~excludedata(kl,imp_imag1,'range',[0 2000]))','smoothingspline')
%figure(4);plot(linspace(eps,6*pi,5000),impimag1_fit(linspace(eps,6*pi,5000)));axis([0 6*pi -1000 1000])
clear k_l kl k eta a C imp_real imp_imag1
save('impfits')

%imp_imag2=(eta/(4*pi)).*sin(kl).*cosint(2*k_l*a*a);

%imp_imag=1i*(imp_imag1+imp_imag2)
