function [ ratio ] = rovert( ratio0,t )
%ROVERT this function calculates the corresponding position of HIO with the
%given time t. ratio0 is the amplitude of the oscillation, t must be less than 4T (the total period of the HIO),
%and the initial position of HIO is considered as 0.
%%%%%%%%%%%%%%%%%%%%%%%
%constants
%%%%%%%%%%%%%%%%%%%%%%%
G=6.67e-11;
%the density of Earth: rho=A-Br^2,according to the Earth model
A=13.0885e3;
B=8.8381e3;
%alpha and K
alpha=1/sqrt(8*pi*G*(A/6-B/20*ratio0^2));
k=B/20*ratio0^2/(A/6-B/20*ratio0^2);
%%%%%%%%%%%%%%%%%%%%%%%
%The calculation of the equation
%a simple dichotomy is used
%%%%%%%%%%%%%%%%%%%%%%%
tmin=0;
tmax=pi/2;
for n=1:20
    tint=(tmin+tmax)/2;
    eq=t-alpha*ellipticF(tint,k);
    if eq<0
        tmax=tint;
    else
        tmin=tint;
    end
end
ratio=sin(tint)*ratio0;
end

