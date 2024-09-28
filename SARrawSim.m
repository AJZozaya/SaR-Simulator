% SARrawSim V. 2.0
% SARrawSim generates raw data from a general pulsed linear frequency-modulated 
% synthetic aperture radar (SAR) operating in stripmap mode. In the first part 
% of the code, the probing area is defined. Then, radar parameters, including the 
% transmitted pulse and time span for the A-scope visualization, are established. 
% Subsequently, a set of point targets is randomly generated. A loop is included 
% in which the radar travel is shown in the form of an animation while checking 
% for targets inside the half-power beam-width antenna pattern (this is equivalent 
% to the detection of targets inside the antenna footprint). Then, indexes of 
% targets inside the antenna footprint are extracted (this is equivalent to 
% target ranging). Later on, the echo is created, Gaussian noise is added, 
% and the result is visualized in A-scope type. Finally, a raw data image is 
% rendered, and metadata and raw data are consolidated for further processing.
% 09/28/2024. A. J. Zozaya

clear all
close all
clc

%% Probing scene definition and drawing
% Scene parameters
Dr=200;                     % [m] r-swath
Du=800;                     % [m] u-swath
Rm=2e3;                     % [m] average range from scene to radar
% Probing scene coordinates
rm=Rm-Dr/2;
rM=Rm+Dr/2;
um=-Du/2;
uM=Du/2;
r_scene=[rm, rm, rm, rM, rM, rM, rM, rm];
u_scene=[um, uM, uM, uM, uM, um, um, um];

%% Radar parameters definition

c=3e8;                          % [m/s] speed of light

% Transmitted pulse definition
f0=2.4e9;                       % [Hz] carrier frequency
dr=1;                           % [m] range resolution
du=0.7;                         % [m] cross-range resolution
B=c/(2*dr);                     % [Hz] bandwidth required for such a range resolution
tau=100/B;                      % [s] pulse duration based on the criterion: B x tau = 100
k=B/tau;                        % [Hz/s] chirp rate
osf=1.4;                        % over sampling factor
Fs=osf*B;                       % [Hz] fast time sampling frequency
Ts=1/Fs;                        % [s] fast-time sampling time
t=-tau/2:Ts:tau/2;              % [s] time support for the base-band pulse transmitted
p=exp(1j*pi*k*t.^2);            % base-band transmitted pulse 
%  Time span of a range-line (= to the time span for A-scope visualization purposes)
ti=2*rm/c;                      % [s] time of arrival of the echo rise edge from the nearest target 
tf=2*rM/c;                      % [s] time of arrival of the echo rise edge from the farthest target 
NofS=ceil((tf+1*tau-ti)/Ts);    % minimun number of samples in a range-line required
nfft=2^ceil(log2(NofS));        % number of samples in a range-line for processing purposes and visualization
t=ti:Ts:ti+(nfft-1)*Ts;         % time span in a range-line for processing purposes and visualization
P=fftshift(fft(p,nfft));        % spectrum of the transmitted pulse for processing purposes
df=Fs/nfft;                     % [Hz] frequency spectrum resolution
f=-Fs/2:df:Fs/2-df;             % [Hz] frequency support for spectrum of the transmitted pulse

% UAV and Radar parameters for trajectory drawing

lambda0=c/f0;                   % [m] carrier wavelength 
Dtheta=lambda0/(2*du);          % [rad] 3 dB horizontal beamwidth
v=10;                           % [m/s] uav speed (from [Zozaya2016])
BD=2*v*Dtheta/lambda0;          % [Hz] Doppler bandwidth
tR=1/(osf*BD);                  % [s] pulse repetition time = slow-time samplig time


% Parameters for antenna pattern ground footprint drawing
m=tan(Dtheta/2);                % [ ] curve parameter

%% Random targets generation

% NofT=10;                        % number of targets
% u_target=um+Du*rand(1,NofT);    % target u-coordinates
% r_target=rm+Dr*rand(1,NofT);    % target r-coordinates
% sigma_target=ones(1,NofT).';    % target sigma  

%% Testing targets

% Four (4) targets
% u_target=[0 350 0 -350];
% r_target=[rm rM rM rm+Dr/2];
% sigma_target=ones(1,4).';      

% One target
% u_target=[0 ];
% r_target=[rm+Dr/2];
% sigma_target=1;      

% Square target plus random point targets
NofT=10;                                        % number of random point targets
u_target_r=um+Du*rand(1,NofT);                  % point target u-coordinates
r_target_r=rm+Dr*rand(1,NofT);                  % point target r-coordinates
sigma_target_r=ones(1,NofT);                    % point target sigmas  
ut=[-20:20];                                    % square target u-coordinates
rt=[1980:2020];                                 % square target r-coordinates
[Ut,Rt]=meshgrid(ut,rt);                        % square target (u,r) mesh coordinates 
u_target=Ut(:).';                               % square target u-coordinates
r_target=Rt(:).';                               % square target r-coordinates
sigma_target=0.5*ones(size(u_target));          % square target sigma
u_target=[u_target,u_target_r];                 % square plus point target u-coordinates
r_target=[r_target,r_target_r];                 % square plus point target r-coordinates
sigma_target=[sigma_target,sigma_target_r].';   % sigma of square plus point targets

%% Radar travel animation and raw data generation
u=-Du/2;                        % [m] initial uav position
n=0;                            % row index dummy initial value
figure(1)
while u<=Du/2;
    n=n+1;

    %% Antenna footprint definition for purposes of drawing
    uam=u+m*rm;
    uaM=u+m*rM;
    ubm=u-m*rm;
    ubM=u-m*rM;
    r_footprint=[rm,rm,rm,rM,rM,rM,rM,rm];
    u_footprint=[ubm,uam,uam,uaM,uaM,ubM,ubM,ubm];

    %% Drawing of probing area, antenna footprint and radar position
    subplot(2,1,1)
    plot(r_target,u_target,'.',0,u,'+',r_scene,u_scene,r_footprint,u_footprint,'LineWidth',1.5)
    xlabel('$r$','Interpreter','latex')
    ylabel('$u$','Interpreter','latex','rotation',0)
    xlim([-10 Rm+Dr])
    ylim([-Du/2-10 Du/2+10])
    % axis tight equal
    
    %% Checking for targets inside the antenna footprint and extraction of their indexes
    ioft=find(u_target(:) < u+m*r_target(:) & u_target(:) > u-m*r_target(:));
    l_ioft=length(ioft);
    if l_ioft==0
        e(n,:)=zeros(1,nfft);
    else
        % Delay insertion
        R=sqrt((u_target(ioft)-u).^2+(r_target(ioft)).^2).';
        D=(2/c).*R-ti;
        
        % Echoes without geometric attenuation and without speckle noise
        E=P.*(sigma_target(ioft)).*(exp(-1j*2*pi*D*(f0+f)));
        
        % Echoes with geometric attenuation
        % E=P.*(sigma_target(ioft)).*(exp(-1j*2*pi*D*(f0+f)))./R;
            
        % Echoes with multiplicative Speckle noise
        % E=P.*(sigma_target(ioft)).*(raylrnd(pi/4,l_ioft,nfft).*exp(1j*(-1*pi+1*pi*rand(l_ioft,nfft))).*exp(-1j*2*pi*D*(f0+f)));

        % Echo formation
        e(n,:)=sum(ifft(fftshift(E),nfft,2),1);
    end

    %% Addition of thermal noise
    % N0=0.001.^2;                              % noise power density with geometric attenuation
    N0=0.1.^2;                                % noise power density without geometric attenuation
    e(n,:)=e(n,:)+sqrt(N0)*(randn(1,nfft)+1j*randn(1,nfft)); 
    
    %% Range-line visualization (A-scope)
    subplot(212)
    plot(t*c/(2e3),real(e(n,:)),'LineWidth', 1.5)
    xlim([rm rM+c*tau/2]./(1e3))
    ylim('padded')
    % ylim(target,limitmode,'tickaligned', 'padded')
    xlabel('$r$','Interpreter','LaTeX')
    ylabel('$e(t)$','Interpreter','LaTeX','rotation',0)
    % box off
    grid on
    grid minor
    ax=gca;
    ax.MinorGridAlpha = 1;                  % Make grid lines less transparent.
    ax.MinorGridColor = [0.1, 0.7, 0.2];    % Dark Green.
    if n==500
        print('simulator.svg','-dsvg')
    end
    %% Radar position update
    u=u+v*tR;
    pause(0.1)
end

%% Raw image 
raw_image=e(:,1:NofS);          % extraction of samples of interest for visualization purposes
r=t(1:NofS)*c/2;                % creating the range support vector for visualization purposes
u=um:v*tR:um+(n-1)*v*tR;        % creating the cross-range support support for visualization purposes
h2=figure(2);
set(gcf, 'WindowState', 'maximized');
subplot(141)
pcolor(r,u,abs(raw_image))
xlabel('$r$','Interpreter','LaTeX')
ylabel('$u$','Interpreter','LaTeX')
shading interp
axis equal tight
subplot(142)
pcolor(r,u,angle(raw_image))
xlabel('$r$','Interpreter','LaTeX')
ylabel('$u$','Interpreter','LaTeX')
shading interp
axis equal tight
subplot(143)
pcolor(r,u,real(raw_image))
xlabel('$r$','Interpreter','LaTeX')
ylabel('$u$','Interpreter','LaTeX')
shading interp
axis equal tight
subplot(144)
pcolor(r,u,imag(raw_image))
xlabel('$r$','Interpreter','LaTeX')
ylabel('$u$','Interpreter','LaTeX')
shading interp
axis equal tight
exportgraphics(h2,'rw.jpg','Resolution',300)

%% Completion of metadata data
P=ifftshift(P);
NofC=ceil((tf-ti)/Ts);          % number of range bins for for processing purposes
NofR=n;                         % number of cross-range bins for processing purposes
targets=[r_target; u_target];
scene=[r_scene; u_scene];

% metadata requiered in general: r, u, raw_image, ti, Ts, NofC
% metadata required for range compression: P,e,NofS
% matadata requiered for CRMC purposes: eR,NofR,dr,tR,r,f0,v
% metadata required for azimuth compression: Dtheta

save rawdata.mat r u raw_image P e ti Ts NofC NofS NofR tR f0 v Dtheta targets scene;

