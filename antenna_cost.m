function [cost] = antenna_cost(variables)
%clc
variables;
freq_up=1500; % Frequency range
freq_low=1000;
numb_of_freq_points=201; 
Z0=50; % Characteristic impedance of transmission line
l=variables(1:length(variables)/3);
d=variables(((length(variables)/3)+1):length(variables)*2/3);
z_normalized=variables((length(variables)*(2/3) +1):end);
z=z_normalized*1;
freq_up=freq_up*10^6;
freq_low=freq_low*10^6;
load impfits

C=0.577215664901532860606512090082402431042; % Euler constant
eta=377; % Intrinsic impedance
a=10^-4; % Dipoles diameter

N=length(d); % Number of elements
kl=zeros(1,N); 
k_l=zeros(1,N);
impedances=zeros(1,N);
reflection=zeros(1,numb_of_freq_points);
cint_par=zeros(1,N);
index=1;



for freq=linspace(freq_low,freq_up,numb_of_freq_points) % Reflection and impedance are calculated for each frequency point
    k=2*pi*freq/(3*(10^8));
    kl(:)=k*l(:);
    k_l(:)=k./l(:);
    cint_par(:)=2*a*a*k_l(:);
    C=0.577215664901532860606512090082402431042; % Euler constant
    v=1:50;
    temp1=zeros(1,length(v));
    cint=zeros(1,N);
    for i=1:N
        temp1(1,:)=((-1).^v(1,:)).*(cint_par(i).^(2.*v(1,:)))./((2.*v(1,:)).*factorial(2.*v(1,:)));
        if cint_par(i)<6
            cint(i)=(C+log(cint_par(i))+sum(temp1));
        end
    end
    impedances(:)=impreal_fit(kl(:)) + (impimag1_fit(kl(:))+  ((eta/(4*pi)).*sin(kl(:))  .*cint(:)   )./(sin(0.5*kl(:)).^2))*1i;
    Z=zeros(N);
    if N>1
        for i=1:N
            for j=(i+1):N
                l1=l(i);
                l2=l(j);
                di=sum(d(i:(j-1)));
                z12m  = @(z) 1i*30*sin(k*(l2/2 - abs(z)))*(exp(-1i*k*sqrt(di*di + (l1/2 - z).^2))/sqrt(di*di + (l1/2 - z).^2) + exp(-1i*k*sqrt(di*di + (-l1/2 - z).^2))/sqrt(di*di + (-l1/2 - z).^2) - 2*cos(k*l1/2)*exp(-1i*k*sqrt(di*di + z.*z))/sqrt(di*di + z.*z));
                out=(quad(z12m,-l2/2,l2/2,10^-1));
                Z(i,j)=out;
            end
        end
        Z=Z+Z'+diag(impedances);
        impedances=sum(Z);
    end
    
    Zeq=Inf;
    for i=1:N
        Zeq=(1/(z(i)*(impedances(i)+ (1i*z(i)*tan(k*d(i))))/(z(i) + 1i*impedances(i)*tan(k*d(i)))))+(1/Zeq);
        Zeq=1/Zeq;
    end
    % Antenna arrangement
    % feed---z(N)---Nth-----(N-1)th---- ... -z(2)----2nd---z(1)---1st
    Zeq;
    reflection(index)=((abs((Zeq-Z0)/(Zeq+Z0)))^2);
    
    index=index+1;
end
freq=linspace(freq_low,freq_up,numb_of_freq_points);
reflection;
%figure(1)
%plot(freq,20*log10(sqrt(reflection)));%axis([10^9 1.5*10^9 -22 -6]);
%grid on;
cost=sum((reflection))/numb_of_freq_points;

end
