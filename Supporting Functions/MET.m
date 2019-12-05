function [T,TZZ,N2,O2,OO,Ar,He,H,EM,dens,dl,p]=MET(z,xlat,xlng,iyr,mn,ida,ihr,min,i1,f10,f10b,gi)

%% Marshall Engineering Thermosphere Model
% Hand converted by Dr. Alexandra Long from the Fortran 77 model written by Mike Hickey
%
%This program is a driving program for the following subroutines
%
%                            TIME                                            
%                           J70                                             
%                                                                            
% The atmospheric model is a modified Jacchia 1970 model and is given in   
% the subroutine J70.  All of the other subroutines were designed to       
% allow flexible use of this model so that various input parameters could  
% be varied within a driving program with very little software development.
%
%   Input Data: 
%       Z = altitude (km) = indata(1)
%       xlat = geocentric latitude (deg) = indata(2)
%       xlng = longitude (deg) = indata(3)
%       iyr = year (yy) = indata(4)
%       mn = month (mm) = indata(5)
%       ida = day (dd) = indata(6)
%       ihr = hour (hh) = indata(7)
%       min = minutes (mm) = indata(8)
%       i1 = geomagnetic index flag = indata(9)
%       f10 = solar radio noise flux = indata(10)
%       f10b = 162-day average f10 = indata(11)
%       gi = geomagnetic activity index = indata(12)
%
%   Output Data: (All output in MKS units)
%       T = Exospheric temperature = outdata(1), K
%       Tzz = Temperature = outdata(2), K
%       A(1) = N2 number density = outdata(3), /m^3
%       A(2) = O2 number density = outdata(4), /m^3
%       A(3) = O number density = outdata(5), /m^3
%       A(4) = A number density = outdata(6), /m^3
%       A(5) = He number density = outdata(7), /m^3
%       A(6) = H number density = outdata(8), /m^3
%       em = Average molecular mass = outdata(9)
%       dens = Total mass density = outdata(10), kg/m^3
%       dl = Log10 mass density = outdata(11)
%       P = Total pressure = outdata(12), Pa

%% Parameters
Rgas=8.31432e3; %J/kmol-K
bfh=440.0;
z = z/1000; % convert to km from m

% Calculations performed for only one latitude, one longitude, and one
% altitude

%% Calculations
[sda,sha,dd,dy,xlat]=TME(mn,ida,iyr,ihr,min,xlat,xlng);

T = TINF(f10,f10b,gi,xlat,sda,sha,dy,i1);

[Tz,A(1),A(2),A(3),A(4),A(5),A(6),EM,dens,DL]=JAC(z,T);

denlg=0;
den=DL;

if z<=170
   denlg=SLV(z,xlat,dd);
end

% 'Fair' helium number density between base fairing height (BFH) and 500 km
if z>=500
    [DL,A(5)]=SLVH(den,A(5),xlat,sda);
elseif z>bfh
    dHel1=A(5);
    dHel2=A(5);
    dLG1=DL;
    dLG2=DL;
    [dLG2,dHel2]=SLVH(dLG2,dHel2,xlat,sda);
    IH=z;
    [A(5),DL]=FAIR5(dHel1,dHel2,dLG1,dLG2,IH);
end
DL=DL+denlg;
dens=10^DL;
xlat=xlat*57.29577951;

%% Set output values
TZZ=Tz;
N2=1e6*(10^A(1));
O2=1e6*(10^A(2));
OO=1e6*(10^A(3));
Ar=1e6*(10^A(4));
He=1e6*(10^A(5));
H=1e6*(10^A(6));
dens=dens*1000;
dl=DL;
p=dens*Rgas*Tz/EM;
end

function [sda,sha,dd,dy,xlat]=TME(mn,ida,iyr,ihr,min,xlat,xlng)
% Inputs:
%   mn = month
%   ida = day
%   iyr = year
%   ihr = hour
%   min = minute
%   xlat = latitude (input-geocentric latitude)
%   xlng = longitude (input-geocentric longitude, -180, +180)
%
% Outputs:
%   sda = solar declination angle (rad)
%   sha = solar hour angle (rad)
%   dd = day number from 1 Jan
%   dy = dd/tropical year

%% Parameters
year=365.2422;
A1=99.6909833;
A2=36000.76892;
A3=0.00038708;
A4=0.250684477;
B1=0.0172028;
B2=0.0335;
B3=1.407;

% pia=3.14159265;
% tpi=6.28318531;
% pi2=1.57079633;
% pi32=4.71238898;
% rad_deg=0.017453293;

pia=pi;
tpi=2*pi;
pi2=pi/2;
pi32=3*pi/2;
rad_deg=pi/180;

iday=[31,28,31,30,31,30,31,31,30,31,30,31];

%% Calculations
xlat=xlat/57.29577951;
yr=iyr;

if mod(iyr,4) == 0
    if mod(iyr,100) ~= 0
        iday(2)=29; % Century not a leap year
    end
else
    iday(2)=28;
end

id=0;
if mn>1
    for iter=1:mn-1
        id=id+iday(iter);
    end
end

dd=id+ida;
dy=dd/year;

% Compute mean Julian date
xmjd=2415020+365*(yr-1900)+dd+(iyr-1901)/4;

% Compute Greenwich mean time in minutes GMT
xhr=ihr;
xmin=min;
gmt=60*xhr+xmin;
fmjd=xmjd-2435839+gmt/1440;

% Compute Greenwich mean position - GP (in rad)
xj=(xmjd-2415020.5)/36525;
gp=mod((A1+A2*xj+A3*xj*xj+A4*gmt),360);

% Compute right ascension point - RAP (in rad)
% first convert geocentric longitude to deg longitude - west neg, + east
if xlng>180
    xlng=xlng-360;
end
RAP=mod((gp+xlng),360);

% Compute celestial longitude - XLS (in rad) - 0 to 2PI
y1=B1*fmjd;
y2=0.017202*(fmjd-3);
xls=mod((y1+B2*sin(y2)-B3),tpi);

% Compute solar declination angle - sda (in rad)
B4=rad_deg*(23.4523-0.013*xj);
sda=asin(sin(xls)*sin(B4));

% Compute right ascension of Sun - ras (in rad) - 0 to 2 pi
% These next few lines do not appear in NASA CR-179359
% They are added here to ensure that argument of asin stays bounded between
% -1 and +1, which could otherwise be effected by roundoff error.
arg = tan(sda)/tan(B4);
if arg>1
    arg=1;
end
if arg<-1
    arg=-1;
end
ras=asin(arg);

% Put ras in same quadrant as xls
ras=abs(ras);
temp=abs(xls);

if temp<=pia && temp>(pi2)
    ras=pia-ras;
elseif temp<=pi32 && temp>pia
    ras=pia+ras;
elseif temp>pi32
    ras=tpi-ras;
end

if xls<0
    ras=-ras;
end

% Compute solar hour angle - sha (in deg)
sha = RAP*rad_deg - ras;
if sha>pia
    sha=sha-tpi;
end
if sha<-pia
    sha=sha+tpi;
end

end

function TE = TINF(f10,f10b,gi,xlat,sda,sha,dy,i1)
% Parameters
xm=2.5;
xnn=3.0;
% pia=3.14159265;
% tpi=6.28318531;
pia=pi;
tpi=2*pi;


%Ci are solar activity variation variable
C1=383.0;
C2=3.32;
C3=1.80;

%Di are geomagnetic variation variables
D1=28.0;
D2=0.03;
D3=1.0;
D4=100.0;
D5=-0.08;

%Ei are semiannual variation variables
E1=2.41;
E2=0.349;
E3=0.206;
E4=6.2831853;
E5=3.9531708;
E6=12.5663706;
E7=4.3214352;
E8=0.1145;
E9=0.5;
E10=6.2831853;
E11=5.9742620;
E12=2.16;

beta=-0.6457718;
gamma=0.7504916;
P=0.1047198;
RE=0.31; %Diurnal factor, KP, F10b, avg

% Solar Activity cariation
TC=C1+C2*f10b+C3*(f10-f10b);

% diurnal variation
eta=0.5*abs(xlat-sda);
theta=0.5*abs(xlat+sda);
tau=sha+beta+P*sin(sha+gamma);

if tau>pia
    tau=tau-tpi;
end
if tau < -pia
    tau=tau+tpi;
end

A1=(sin(theta)).^xm;
A2=(cos(eta)).^xm;
A3=(cos(tau/2)).^xnn;
B1=1+RE*A1;
B2=(A2-A1)/B1;
TV=B1*(1+RE*B2*A3);
TL=TC*TV;

% Geomagnetic variation
if i1==1
    TG=D1*gi+D2*exp(gi);
else
    TG=D3*gi+D4*(1-exp(D5*gi));
end

% Semiannual variation
G3=0.5*(1+sin(E10*dy+E11));
G3=G3.^E12;
tau1=dy+E8*(G3-E9);
G1=E2+E3*(sin(E4*tau1+E5));
G2=sin(E6*tau1+E7);
TS=E1+f10b*G1*G2;

% Exospheric Temperature
TE=TL+TG+TS;
end

function [TZ,AN,AO2,AO,AA,AHE,AH,EM,dens,DL]=JAC(Z,T)
% Parameters
AV=6.02257e23;  % Avogadro's number [1/mol]
QN=0.78110;
QO2=0.20955;
QA=0.009343;
QHE=1.289e-5;
Rgas=8.31432;
T0=183;
% pia=3.14159265;
pia=pi;

alpha=[0,0,0,0,-0.380,0];
EI=[28.0134,31.9988,15.9994,39.948,4.0026,1.00797];

TX=444.3807+0.02385*T-392.8292*exp(-.0021357*T);
A2=2.*(T-TX)/pia;
TX_T0=TX-T0;
T1=1.9*TX_T0/35;
T3=-1.7*TX_T0/(35.^3);
T4=-0.8*TX_T0/(35^4);
TZ=TEMP(Z,TX,T1,T3,T4,A2);

%% Section 1
A=90;
D=min(Z,105);

% Integrate gM/T from 90 to minimum of Z or 105 km
R=GAUSS(A,D,1,TX,T1,T3,T4,A2);

% The number 2.1926e-8=density*temp/mean molecular weight at 90 km

EM=Mol_WT(D);
TD=TEMP(D,TX,T1,T3,T4,A2);

dens=2.1926e-8*EM*exp(-R/Rgas)/TD;
factor=AV*dens;
par=factor/EM;
factor=factor/28.96;

% For altitudes below and at 105 km calculate the individual specie number
% densities from the mean molecular weight and total density

if Z <= 105
    DL=log10(dens);
    AN=log10(QN*factor);
    AA=log10(QA*factor);
    AHE=log10(QHE*factor);
    AO=log10(2*par*(1-EM/28.96));
    AO2=log10(par*(EM*(1+QO2)/28.96-1));
    AH=0;
    return
end

%% Section 2
% This section is only performed for altitudes above 105 km
% Note that having reached this section means that D in section 1 is 105 km
% Calculate individual specie number densities from the total density and
% mean molecular weight at 105 km

DI(1)=QN*factor;
DI(2)=par*(EM*(1+QO2)/28.96-1);
DI(3)=2*par*(1-EM/28.96);
DI(4)=QA*factor;
DI(5)=QHE*factor;

%Integrate g/T from 105 km to Z km
R=GAUSS(D,Z,2,TX,T1,T3,T4,A2);

for iter=1:5
    DIT(iter)=DI(iter)*(TD/TZ)^(1+alpha(iter))*exp(-EI(iter)*R/Rgas);
    if DIT(iter)<=0
        DIT(iter)=1e-6;       
    end
end

%% Section 3
% This section calculates atomic hydrogen densities above 500 km altitude
% Below this altitude, H densities are set to 10^-6

if Z>500
   A1=500;
   S=TEMP(A1,TX,T1,T3,T4,A2);
   DI(6)=10^(73.13-39.4*log10(S)+5.5*log10(S)*log10(S));
   R=GAUSS(A1,Z,7,TX,T1,T3,T4,A2);
   DIT(6)=DI(6)*(S/TZ)*exp(-EI(6)*R/Rgas);
else
    DIT(6)=1;
end

% For altitudes greater than 105 km, calculate total density and mean
% molecular weight from individual specie number densities.
dens=0;
for iter2=1:6
    dens=dens+EI(iter2)*DIT(iter2)/AV;
end

EM=dens*AV/(DIT(1)+DIT(2)+DIT(3)+DIT(4)+DIT(5)+DIT(6));
DL=log10(dens);
AN=log10(DIT(1));
AO2=log10(DIT(2));
AO=log10(DIT(3));
AA=log10(DIT(4));
AHE=log10(DIT(5));
AH=log10(DIT(6));
end

function T=TEMP(alt,TX,T1,T3,T4,A2)
% Calculates the temperature at altitude atl using equation (10) for
% altitudes between 90 and 125 km and equation (13) for altitudes greater
% than 125 km, from SAO Report

BB=4.5e-6;
U=alt-125;

if U>0
    T=TX+A2*atan(T1*U*(1+BB*(U^2.5))/A2);
else
    T=TX+T1*U+T3*U^3+T4*U^4;
end
end

function MW=Mol_WT(A)
% Calculates the molecular weight for altitudes between 90 and 105 km
% according to equation (1) of SAO report 313. Otherwise, MOL_WT is set to
% unity.
B=[28.15204,-0.085586,1.284e-4,-1.0056e-5,-1.021e-5,1.5044e-6,9.9826e-8];

if A>105
    MW=1;
else
    U=A-100;
    MW=B(1);
    for iter=2:7
        MW=MW+B(iter)*U.^(iter-1);
    end
end
end

function R=GAUSS(Z1,Z2,Nmin,TX,T1,T3,T4,A2)
% Subdivide total integration-altitude range into intervals suitable for
% applying Gaussian Quadrature, set the number of points for integration
% for each sub-interval, and then perform Gaussian Quadrature

altmin=[90,105,125,160,200,300,500,1500,2500];
NG=[4,5,6,6,6,6,6,6];

% Coefficients for Gaussian Quadrature
C=[.5555556,.8888889,.5555556,.00,.00,.00,.00,.00;...
    .3478548,.6521452,.6521452,.3478548,.00,.00,.00,.00;...
    .2369269,.4786287,.5688889,.4786287,.2369269,.00,.00,.00;...
    .1713245,.3607616,.4679139,.4679139,.3607616,.1713245,.00,.00;...
    .1294850,.2797054,.3818301,.4179592,.3818301,.2797054,.1294850,.00;...
    .1012285,.2223810,.3137067,.3626838,.3626838,.3137067,.2223810,.1012285];

% Abscussas for Gaussian Quadrature
x=[-.7745967,.00,.7745967,.00,.00,.00,.00,.00;...
    -.8611363,-.3399810,.3399810,.8611363,.00,.00,.00,.00;...
    -.9061798,-.5384693,.00,.5384693,.9061798,.00,.00,.00;...
    -.9324695,-.6612094,-.2386192,.2386192,.6612094,.9324695,.00,.00;...
    -.9491079,-.7415312,-.4058452,.00,.4058452,.7415312,.9491079,.00;...
    -.9602899,-.7966665,-.5255324,-.1834346,.1834346,.5255324,.7966665,.9602899];

C=C'; % Need to account for Fortran being column major in array population
x=x';

R=0;
for k=Nmin:8;
    ngauss=NG(k);
    A=altmin(k);
    D=min(Z2,altmin(k+1));
    RR=0;
    del=0.5*(D-A);
    J=ngauss-2;
    for it=1:ngauss
        Z=del*(x(it,J)+1)+A;
        grav=9.80665/((1+Z/6.356766e3)^2);
        MW=Mol_WT(Z);
        T=TEMP(Z,TX,T1,T3,T4,A2);
        RR=RR+C(it,J)*MW*grav/T;
    end
    RR=del*RR;
    R=R+RR;
    if D == Z2
        return
    end
end
end

function den=SLV(alt,xlat,day)
% Computes the seasonal-latitudinal variation of density in the lower
% thermosphere in accordance with L. Jacchia, SAO 332, 1971.
%This affects the densitites between 90 and 170 km. This function need not
%be called for densitites above 170 km, because no effect is observed.
% 
% The variation should be computed after the calculation of density due to
% temperature variations and the density (den) must be in the form of a
% base 10 log. No adjustements are made to the temperature or constituent
% number densitites in the region affected by this variation.

%initialize density (den) = 0
den=0.0;

% check if altitude exceeds 170 km
if alt>170
    return
end

% Compute density change in lower thermosphere
Z=alt-90;
X=-0.0013*Z*Z;
Y=0.0172*day+1.72;
P=sin(Y);
SP=(sin(xlat))^2;
S=0.014*Z*exp(X);
D=S*P*SP;

% Check to compute absolute value of 'xlat'
if xlat<0
    D=-D;
end
den=D;
end

function [den,denHE]=SLVH(den,denHE,xlat,SDA)
% Computes the seasonal-latitudinal variation of the helium number density
% according to L. Jacchia, SAO 332, 1971. This correction is not important
% below 500 km.

D0=10^denHE;
A=abs(0.65*(SDA/0.40909079));
B=0.5*xlat;

% Check to compute absolute value of 'B'
if SDA<0
    B=-B;
end

% Compute X, Y, DHE, and denHE
X=0.7854-B;
Y=(sin(X))^3;
DHE=A*(Y-0.35356);
denHE=denHE+DHE;

% Compute helium number density change
D1=10^denHE;
del=D1-D0;
rho=10^den;
drho=6.646e-24*del;
rho=rho+drho;
den=log10(rho);
end

function [FdHel,FdLG]=FAIR5(dHel1,dHel2,DLG1,DLG2,IH)
% Fairs between the region above 500 km, which invokes the
% seasonal-latitudinal variation of the helium number density (SLVH), and
% the region below, which does not invoke any seasonal-latitudinal
% variation at all.

cz=[1,0.9045085,0.6545085,0.3454915,0.0954915,0];
IBFH=440;

% Height Index
I=floor((IH-IBFH)/10+1);

% Non-SLVH fairing coefficient
czi=cz(I);

% SLVH fairing coefficient
szi=1-czi;

% Faired density
FdLG=(DLG1*czi)+(DLG2*szi);

% Faired helium number density
FdHel=(dHel1*czi)+(dHel2*szi);

end