clear all


tmax=3000000000; %max length of simulation
depth=1; %what it says, meters

mortB=0.01;         %mortality of mpb per timestep. Should be dynamic. 1/day (M&F, JGR 2012)
growthrateB=10*0.0062;     %1/day, max growth of mpb, should depend at least on %mud and nutrients, units: 1/day from Niathany 
grazingB=0.0002;       %1/day, ?????????? 
K=1200;            %mgCha/m2carrying capacity, max mpb biomass. %Ecolog. Modelling 2016 Sinha 
Xb=0;              %minimum background mpb biomass, mgCha/m2

dfine=0.000012;	% diameter of fine sed, in m */
dcoarse=0.000125;   % diameter of coarse sed, in m */
rho=1000;             	% water density, kg/m^3 */
rhoS=2600;           	% sediment density, kg/m^3 */
%Shieldscriticalcoarse has been evaluated from the shields diagram
%for dcoarse=0.000125, see also page 106 of Dynamics of marine sand, soulsby 1997*/
Shieldscriticalcoarse=0.065;
g=9.81;  %gravity m/s2

growthB_mac_max=0.7;       %1/day, max growth of macomona. 
K_mac=1200;            %mgC/m2carrying capacity, max number of bugs possible, should depend at least on d50 (simon's paper). 
Xb_mac=0;              %minimum number of bugs, mgCha/m2
wetW2Carbon=0.05*(10^-3)*10^6; %mg C per m2 per individual (from individuals to carbon)
%assumig a shell lenth of 1cm, gives a volume of 0.01*0.01*0.01 m3
%(10^-3 factor in eq below). We assume the specific weight of clams
%is 1000kg/m^3 and 0.05 is the wet weight (AFDW in Ricciardi and
%Bourget, 1998 meps) to carbon, convert from Kg to mg
%**********CHECK DAN PRATT THESIS, PAGE 62

mortB_mac=0.01;         %mortality of macomona per timestep. Should be dynamic. ind/day from judi spreadsheet ...
%We actually realized that Judi's approach results in too much of a large
%contribution and makes the population go negative so ...
% we decided to have just a percentage of the macomona
%mortB_mac=mortB_mac*wetW2Carbon; %mgC/m2ind

%initial and boundary conditions
% B=400*ones(tmax/100,1);           %biomass mgC/m2
% B_mac=300*ones(tmax/100,1);     %biomass mgC/m2
% D_mac=10.0*ones(tmax/100,1);     %detritus 
% D_MPB=10.0*ones(tmax/100,1);     %detritus 
% N=1000000.0*ones(tmax/100,1);     %nutrient 
B(1)=400;           %biomass mgC/m2
B_mac(1)=300;     %biomass mgC/m2
D_mac(1)=10.0;     %detritus 
D_MPB(1)=100.0;     %detritus 
N(1)=100.0;     %nutrient
IN=1000;                   %external source of nutrient. It needs to have the same units as N & Kn
density_mac=0.1; %this can only change between 0 and 1

%Light effect on growth rate (it would be nice to link this to tidal range)
I0=300;          %I0=300miuE/m2/s from M&F
kd=0.4;    %changes depnding on turbidity (mariana) units: 1/meters
I=I0*exp(-kd*depth);
Ki=20;          %Ki=20miuE/m2/s from Bouletreau et al. 2008..probably not for local species
light=I/(I+Ki); %light reduction term

%Mud content effect on growth rate (set to 1 for now. The growth should
%actually depend on nutrient availability which is in turn a function of
%the mud content)
% mudpc=[0.1:0.1:1];
%mudpcf=0.9;
%mud_content=2*(1-1./(1+exp(-mudpcf))); %controls the flux of nutrients out
%of the sediment... not sure how we came up with that
mud_content=0.1;
mud4plot(1)=0.0;
%This is a pore water efficiency factor. More flushing with sand, less
%flushing if it's all mud. Based on hydraulic conductivity?

Kn=10;      %micrograms of Nitrogen/l, from Hydrobiologia, Naithani 2016. Biological uptake, species-dependent
coeff_BD=0.16;  %original value Naithani et al 2016, 1/day
conversion_C2N=5.88; %Naithani et al.
depthbenthiclayer=0.01; %Ecolog. Modelling 2016 Sinha 
%K_uptakenutrient=1.6; %same unit as nutrents mmol*m^-2 day^-1. Data from Sundback & Granelli
K_uptakenutrient=50; %same unit as nutrents mmol*m^-2 day^-1. Data from Sundback & Granelli
coefflushing=0.00001;

%initialization
B_new=B(1);B_mac_new=B_mac(1);
D_MPB_new=D_MPB(1); D_mac_new=D_mac(1); N_new=N(1);N_new_mac=N(1)/2;N_new_MPB=N(1)/2;
dt=0.00001;
dtsave=1000;
isave=0;

 for i=2:tmax
%for i=2:500

N_old=N_new; 
N_old_mac=N_new_mac;
N_old_MPB=N_new_MPB;
B_old=B_new;B_mac_old=B_mac_new;
D_MPB_old=D_MPB_new; D_mac_old=D_mac_new;

%        depthfine=mud_content; depthcoarse=1-mud_content;
%        if mud_content>=1; display('100% muddy'); mud_content=1; end           
%        d50 = dfine*(mud_content) + dcoarse*(1-mud_content); 	%in meters
%         taocriticalcoarse=Shieldscriticalcoarse*g*(rhoS-rho)*dcoarse; %in N/m2
%         taocriticalmax=5*Shieldscriticalcoarse;
%         %Shieldscriticalfine assumes unconsolidated mud, see whitehouse and soulsby p.152 (N/m2)
%         taocriticalfine=0.1; %Shieldscritical is evaluated using Chesher and Ockenden, see whitehouse and soulsby p.152
%         if mud_content<0.2
%             taocritical=taocriticalcoarse+(taocriticalmax-taocriticalcoarse)*mud_content/0.2;
%         else
%             taocritical=taocriticalmax-(taocriticalmax-taocriticalfine)*(mud_content-0.2)/0.8;
%         end
%         tao=0.1;
%         %below we assume that chla content is proportional to our biomass
%         %B. The proportionality is 0.1 at the moment
%         stabchlafactorcoeff=1-0.0018*(0.1*B_old); %from Paarlberg et al ECSS 2005 (eq. 13)  */
%         %Waeles et al (Hydrobiologia 2007, p.74). See also whitehouse and
%         %soulsby p.152; FineErosion is in kg
%         FineErosion=stabchlafactorcoeff*mud_content*0.002*(2*tao/taocritical-1);
%         %we assume coarse sediment does not change. coarse sediment is not eroded or deposited
%         %below we set the depsotion of fine
%         FineDeposition=FineErosion;
%         %now we calculate the change in d50
%         timestep=1; 
%         CELL_WIDTH=1;
%         depth=(FineDeposition-FineErosion)*timestep*CELL_WIDTH*CELL_WIDTH/rhoS;
%         depthactivelayer=0.1; %m
%         newdepthfine=depthfine+depth;
%         newdepthcoarse=depthcoarse-depth;
%         if newdepthcoarse<0
%             display('woooooo  no more coarse sediment')
%         end
%        mud_content=newdepthfine/(newdepthfine+newdepthcoarse);  
       mud_content=0.1;      
       % mud_content=0.1+0.7*i/(tmax); 
       %IN=100*(0.1+816*i/(tmax));% this is the exact condition below which we get a change with 
        
        %biomass dynamics
        %the reduced uptake when sediment is muddier. Look at 14.70, p. 479
        %of the big modelling book
         if IN<N_old_MPB
            growthB=growthrateB*light*(N_old_MPB)/(K_uptakenutrient+(N_old_MPB)); %we need to check the uptake value
         else
             growthB=growthrateB*light*(IN)/(K_uptakenutrient+(IN)); %we need to check the uptake value  
         end    
        growthB_mac=growthB_mac_max*B_old/K;
        %with carrying capacity
        
        %MACOMONA DYNAMICS
K_mac_new=(21.6+0.146*(mud_content*100)-0.005*((mud_content*100)^2))*10000/353.4; % sIMON ET AL 2003 meps. uNITS ARE INDIVIDUALS/QUADRAT (0.25^2)

        B_new=B_old+dt*(growthB*B_old*(1-B_old/K)-mortB*(B_old-Xb)-growthB_mac*B_mac_old*(1-B_mac_old/K_mac_new));
                
        %without carrying capacity
        %B(i)=B(i-1)+dt*(growthB*B(i-1)-mortB*(B(i-1)-Xb)-growthB_mac*B_mac(i-1)*(1-B_mac(i-1)/K_mac(i-1)));
        
        if mud_content>=0.8; K_mac=1; end
        
        %we need to fix this bit
       % K_mac_carb(i)=wetW2Carbon*K_mac(i-1); %mg C per m2 
       if N_old_mac>333333
           mortB_mac=0.1; 
       end
    B_mac_new=B_mac_old+dt*(growthB_mac*B_mac_old*(1-B_mac_old/K_mac_new)-mortB_mac*B_mac_old);
        %D(i)=coeff_BD*mortB*B(i-1)-lostD*D(i-1)-convertedD*D(i-1);
        %coeff_BD could be different for mac and mpb
        
    %coefflushing=c1*mud_content+c2*B_mac_old;
    
    D_MPB_new=D_MPB_old+dt*(mortB*B_old-coeff_BD*D_MPB_old);
    D_mac_new=D_mac_old+dt*(mortB_mac*B_mac_old-coeff_BD*D_mac_old);
    N_new_MPB=N_old_MPB+dt*(coeff_BD*D_MPB_old-growthB*B_old*(1-B_old/K))/depthbenthiclayer/conversion_C2N+dt*N_old_mac*coefflushing; 
    N_new_mac=N_old_mac+dt*(coeff_BD*D_mac_old)/depthbenthiclayer/conversion_C2N-dt*N_old_mac*coefflushing; 
    N_new=N_new_MPB+N_new_mac;
    if rem(i,dtsave*1000)==0
        i
    end
    if rem(i,dtsave)==0
    isave=isave+1;
    mud4plot(isave)=mud_content;
    K_mac(isave)=K_mac_new;
    B(isave)=B_new;
    B_mac(isave)=B_mac_new;
    D_MPB(isave)=D_MPB_new;
    D_mac(isave)=D_mac_new;
    N(isave)=N_new; 
    N_MPB(isave)=N_new_MPB;
    N_mac(isave)=N_new_mac;
    end
% [N(i-1) coeff_BD*D_MPB(i-1) coeff_BD*D_mac(i-1) -growthB*B(i-1)]
% return
        if N_new<0 N_new=0; N(isave)=0; display('no more N'); return; end
        if N_new_MPB<0 N_new_MPB=0; N_MPB(isave)=0; end %display('no more N_new_MPB'); end
        if N_new_mac<0 N_new_mac=0; N_mac(isave)=0; end %display('no more N_new_mac'); end
        
        %****WARNING !!!!! WE NEED TO CHECK
        %depthbenthiclayer/conversion_C2N (28/8/2018)
end
figure
subplot(3,2,1),plot(B/conversion_C2N/depthbenthiclayer,'r'),ylabel('MPB')
subplot(3,2,2),plot(D_mac/conversion_C2N/depthbenthiclayer,'g'),ylabel('Detritus');hold on;plot(D_MPB/conversion_C2N/depthbenthiclayer,'b');
legend('MAC','MPBs')
subplot(3,2,3),plot(N,'b'),ylabel('Nutrient')
subplot(3,2,4),plot(B_mac,'m'),ylabel('Macomona')
subplot(3,2,5),plot(mud4plot,'m'),ylabel('Mud')
subplot(3,2,6),plot(N_MPB,'b'),ylabel('N_{MPB}')
hold on
subplot(3,2,6),plot(N_mac,'r'),ylabel('N_{mac}')
legend('MPB','mac')

%conference figure1
n=length(B_mac);
figure
p=plot([1:n]*dt*1000,B_mac,'m');set(p,'linewidth',3.5)
hold on
p=plot([1:n]*dt*1000,B,'g');set(p,'linewidth',3.5,'color',[0 0.8 0])
legend('Macomona Biomass','MPB Biomass')
xlabel('Time (days)')
ylabel('Biomass (mgC/{m^2})')
figure
p3=plot([1:n]*dt*1000,N_mac,'m');set(p,'linewidth',3.5)

