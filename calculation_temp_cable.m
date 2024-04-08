%% Calculations for noisy cable paper
clear all; 
%% constants
k_b=1.38e-23; %Boltzmann const
B=12.5e3; %bandwidth 

%% cable prop constant

R = 17.6e-3; % (*Resistance, per unit length may have to adjust but OK*)
L = 0.2e-9; %(*Inductance*)
G = 10e-9; %(*Conductance*)
C = 80.3e-12; %(*Capacitance*)
w = 2*pi*150e6;% (*Angular frequency*)

kc=sqrt((R+1i*w*L)*(G+1i*w*C)); %prop constant 
k=real(kc); 
y=imag(kc); 
Zc=sqrt((R+1i*w*L)./(G+1i*w*C)); 

a=0;
b=.1; %length of cable in meters 

%% Temperature along cable for non uniform T 

% linear temp profile 
m=90; 
T0=290; 

% Gaussian temp profile 
mu=300; 
sig=10; 

% sinusodial temp profile 

%% numerical integration for noise temp 

%R matrix bit linear temp profile: 
funr11 = @(x) R.*(m.*x + T0).*cos(k*x).*cos(conj(k)*x);
funr12 = @(x) R.*(m.*x + T0).*(-1i.*cos(k*x).*sin(conj(k)*x))./conj(Zc);
funr21 = @(x) R.*(m.*x + T0).*(1i.*sin(k*x).*cos(conj(k)*x))./Zc;
funr22 = @(x) R.*(m.*x + T0).*(sin(k*x).*sin(conj(k)*x))./(Zc.*conj(Zc));

%numerical integration for R matrix bit 
F={funr11 funr12; funr21 funr22};
r=zeros(2,2);
for i=1:2
    for j=1:2
        r(i,j)=integral(F{i,j},a,b);
    end
end



%G matrix bit with linear temp profile 
fung11 = @(x) G.*(m.*x + T0).*(Zc.*conj(Zc)).*sin(k*x).*sin(conj(k)*x);
fung12 = @(x) G.*(m.*x + T0).*(1i.*Zc.*sin(k*x).*cos(conj(k)*x))./conj(Zc);
fung21 = @(x) G.*(m.*x + T0).*(-1i.*conj(Zc).*cos(k*x).*sin(conj(k)*x));
fung22 = @(x) G.*(m.*x + T0).*(cos(k*x).*cos(conj(k)*x));


%calculating the integral for G matrix bit 
H={fung11 fung12; fung21 fung22};
g=zeros(2,2);
for i=1:2
    for j=1:2
        g(i,j)=integral(H{i,j},a,b);
    end
end

%% correlation matrix 

cor=4.*k_b.*B.*(g+r);





