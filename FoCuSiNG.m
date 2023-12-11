% FoCuSiNG
% FoCuSiNG implements the Range-Doppler Algorithm on the rawdata obtained
% from SARrawSim.m improved for visualization purposes. 10/24/2023. A. J.
% Zozaya

clear all
% close all
clc

%% Loading raw data and metadata and general parameters definition
c=3e8;                                  % [m/s] speed of light
load rawdata.mat;
lambda0=c/f0;                           % [m] carrer wavelength

h3=figure(3);
set(gcf, 'WindowState', 'maximized');

%% raw image

subplot(151)
pcolor(r,u,abs(raw_image))
xlabel('$r$','Interpreter','LaTeX')
ylabel('$u$','Interpreter','LaTeX')
title(['raw data $(r,u)$ domain image'],'Interpreter','LaTeX')
shading interp
axis equal tight

%% Range compression

eR=RaNGeC(P,e,NofS);

Dr=c*Ts/2;                          % range bin size
r=ti*c/2:Dr:ti*c/2+(NofS-1)*Dr;     % [m] ultimate range support vector
subplot(152)
pcolor(r,u,abs(eR))
xlabel('$r$','Interpreter','LaTeX')
ylabel('$u$','Interpreter','LaTeX')
title(['range-compressed $(r,u)$ domain  image'],'Interpreter','LaTeX')
shading interp
axis equal tight

%% Range Cell Migration Correction

[ERA,ERCA,fD]=RCMC(eR,NofR,Dr,tR,r,f0,v);

subplot(153)
pcolor(r,fD.*1,abs(ERCA))
xlabel('$r$','Interpreter','LaTeX')
ylabel('$f^\prime$','Interpreter','LaTeX')
title(['CRMC $(r,f^\prime)$ domain  image'],'Interpreter','LaTeX')
shading interp
% axis equal tight

%% Azimuth compression

eRA=azimuthC(ERCA,v,fD,r,lambda0);

eRA=eRA(1:NofR,1:NofC);
r=ti*c/2:Dr:ti*c/2+(NofC-1)*Dr;     % [m] image range support vector
u=u(1:NofR);

subplot(154)
pcolor(r,u,abs(eRA))
% pcolor(abs(eRA))
xlabel('$r$','Interpreter','LaTeX')
ylabel('$u$','Interpreter','LaTeX')
title('$(r,u)$ domain image','Interpreter','LaTeX')
shading interp
axis equal tight

subplot(155)
plot(scene(1,:),scene(2,:),targets(1,:),targets(2,:),'.')
% pcolor(abs(eRA))
xlabel('$r$','Interpreter','LaTeX')
ylabel('$u$','Interpreter','LaTeX')
title('target $(r,u)$ domain','Interpreter','LaTeX')
% shading interp
axis equal tight

exportgraphics(h3,'focusing.jpg','Resolution',300)


function eR=RaNGeC(P,e,NofS)
% eR=RaNGeC(P,e,TofC) compresses range-lines (rows) from SAR raw data
% e matrix. 10/20/2023. A. J. Zozaya

nfft=size(e,2);
ER=fft(e,nfft,2);                % DFT of rows
H=conj(P);                       % range matched filter frequency response
eR=ifft(H.*ER,nfft,2);           % range compression in the frequency domain and returning back to range domain

eR=eR(:,1:NofS);
end

function [ERA,ERCA,fD]=RCMC(eR,NofR,dr,tR,r,f0,v)
% [ERA,ERAC,fD]=RCMC(eR,NofR,dr,tR,r,f0,v) corrects the range cell
% migration of raw data focused in range. 
% 10/20/2023. A. J. Zozaya 

lambda0=(3e8)/f0;               % [m] carrier wavelength 

%% DFT of columns of e
nfftA=2^(ceil(log2(NofR)));     % nfftA is the size of the column DFT
% nfftA=2*nfftA;
ERA=fft(eR,nfftA,1);            % DFT of eR's columns with fD: 0 <= fD <= fsA.
ERA=fftshift(ERA,1);            % DFT of eR's columns with fD: -fsA/2 <= fD <= fsA/2.
fsA=1/tR;                       % [Hz] sampling frequency in the slow-time domain.
dfD=fsA/nfftA;                  % [Hz] frequency step in the Doppler frequency fD domain.
fD=-fsA/2:dfD:fsA/2-dfD;        % [Hz] Doppler frequency support vector.

%% Interpolation function definition

for n=1:16
hr(n,:)=sinc([-3:4]-n/16);
end

%% Range cell migration correction itself

fD=fD(:);
ncc=size(ERA,2)-20;
ERCA=ERA;
 for j=4:ncc
       Delta_R=r(j).*(lambda0.^2).*(fD(:).^2)./(8*v^2);
       Dj=Delta_R./dr;
       for i=1:nfftA
           m=floor(Dj(i));
           D=Dj(i)-m;
           n=round(D*16);
           if n~=0
            ERCA(i,j)=ERA(i,j+m-3:j+m+4)*hr(n,:).';
           else
            ERCA(i,j)=ERA(i,j+m);
           end
       end
 end
end

function eRA=azimuthC(ERCA,v,fD,r,lambda0)
% eRA=azimuthC(ERCA,v,fD,r,lambda0) compresses 
% azimuth lines (columns) from SAR raw data matrix ERCA. Updated on 
% October 24, 2023. Prof. A. Zozaya. 

mn=size(ERCA);
nfftA=mn(1);
n=mn(2);

for i=1:n
    r0=r(i);
    K_A=(2*v^2)/(r0*lambda0);                   % croos-range (azimuth) space-chirp rate
    HA=exp(-1j*pi*(fD.^2)./K_A);                % azimuth matched filter frequency response
    eRA(:,i)=ifft(HA.*ERCA(:,i),nfftA);         % compression of the i-th column and returning back to cross-range domain
end
end

