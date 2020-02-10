%this transcipt first calculates the r(t) sequence, then convert it into
%deltag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=6.67e-11;
A=13.0885e3;
B=8.8381e3;
ratio0max=1221.5/6371; %doesn't use in the calculation
ratio0=0.1;
alpha=1/sqrt(8*pi*G*(A/6-B/20*ratio0^2));
k=B/20*ratio0^2/(A/6-B/20*ratio0^2);
%one period of the HIO is equal to 4T
T=alpha*ellipticF(pi/2,k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%settings of the required data set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number=4000; %number of total points
frequency=60; %sampling time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the r(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=linspace(0,(number-1)*frequency,number);
r=zeros(1,number);
for n=1:number
    tt=t(n);
    %the function "rovert" calculates the position of 0-T
    %thus a convertion is needed
    m=floor(tt/4/T);
    tleft=tt-m*4*T;
    mleft=floor(tleft/T);
    tleft=tleft-mleft*T;
    if mleft==0
        r(n)=rovert(ratio0,tleft);
    end
    if mleft==1
        r(n)=rovert(ratio0,T-tleft);
    end
    if mleft==2
        r(n)=-rovert(ratio0,tleft);
    end
    if mleft==3
        r(n)=-rovert(ratio0,T-tleft);
    end
end
%convert r into real coordinates
r=r*6371e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the deltag(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rotation of the Earth: omega0
omega0=2*pi/3600/24;
%mass of Earth: me
me=5.972e24;
%mass of HIO: mh
mh=1e15;
deltag=zeros(1,number);
for n=1:number
    %alpha and beta are described in fig.2
    alpha=omega0*t(n);
    rtemp=sqrt(r(n)^2+6371e3^2-2*r(n)*6371e3*cos(alpha));
    cosbeta=(rtemp^2+6371e3^2-r(n)^2)/2/rtemp/6371e3;
    deltag1r=G*mh/rtemp/rtemp*cosbeta;
    deltag1t=G*mh/rtemp/rtemp*sqrt(1-cosbeta^2);
    mt=4*pi*(A/3*r(n)^3-B/5*r(n)^5/6371e3^2);
    if r(n)==0
        deltag2=0;
    else
        deltag2r=G*mh*mt/me/r(n)/r(n)*cos(alpha);
        deltag2t=G*mh*mt/me/r(n)/r(n)*sin(alpha);
    end
    deltagr=deltag1r+deltag2r;
    deltagt=deltag1t-deltag2t;
    deltag(n)=sqrt((9.78+deltagr)^2+deltagt^2)-9.78;
end
deltag=deltag-ones(1,number)*deltag(1);