%%%%%%2U sat heat flux analysis code for particular value of beat , Q_sun,Q_IR%%%%%%%%
%%%%%%%%Different cases can be solved by substituting
%%%%%%%%values%%%%%%%%%%%%%%%%
%%%%%%heat flux an dtemperature calculation%%%%
clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Re=6400; %radius of earth
z=500;% distance of satellitre from earth
rho=0.3; % albedo coefficient
Q_sun=1367;%solar constant 
Q_ir=244; % earth Ir
alpha=1; %absorbtivity
emi=0.7;%emmisivity 
sigma=5.67e-8;
beta = 30;%beta angle
theta = 1:360 ;%orbit angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%view factor%%%%
VF_nadir =(Re/(Re+z))^2; % Earth Facing
a=sqrt(1-VF_nadir);
b=2*asin(a);
c=sin(b);
VF_walls = (pi-b-c)/(2*pi); % Perpendicular to nadir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_entry = 180-asind(Re/(Re+z));
theta_exit = 180+asind(Re/(Re+z));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%area%%%%%%%%%%%%%%%%%
AZ = 0.0227;
AN = 0.0227;
AP=0.01;
AF=0.0227;
AS=0.01;
AA=0.0227;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Q_solar%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta=1:360
Q_solarZ=Q_sun*AZ*cosd(beta)*cosd(theta);%%%%zenith face
for theta =90:270
    Q_solarZ(theta)=0;
end
Q_solarNA= -Q_sun*AN*cosd(beta)*cosd(theta);
for theta=1:360
 if (theta >=90&&theta<theta_entry)
Q_solarNA(theta)= -Q_sun*AN*cosd(beta)*cosd(theta);
else
    Q_solarNA(theta)=0;
 end
end
Q_solarNB= -Q_sun*AN*cosd(beta)*cosd(theta);
for theta=1:360
 if (theta >theta_exit&&theta<=270   )
Q_solarNB(theta)= -Q_sun*AN*cosd(beta)*cosd(theta);
else
    Q_solarNB(theta)=0;
 end
end
Q_solarN=Q_solarNA+Q_solarNB;%%%%%%%%%%%%%%nadir face
Q_solarF=-Q_sun*AF*cosd(beta)*sind(theta);%%%%%%%%%%%%%%forward
for theta=1:360
    if (theta >theta_exit&&theta<=360)
        Q_solarF(theta)=-Q_sun*AF*cosd(beta)*sind(theta);
    else
    Q_solarF(theta)=0;
    end
end
theta=1:360
Q_solarA= Q_sun*AA*cosd(beta)*sind(theta);%%%%%%%%%%%%%%%%%%aft face 
for theta=1:360
    if (theta >=1&&theta<theta_entry)
        Q_solarA(theta)=Q_sun*AA*cosd(beta)*sind(theta);
    else
    Q_solarA(theta)=0;
    end
end

for theta=1:360
if beta>=0
 if (theta <theta_entry||theta>theta_exit)  
Q_solarP(theta)=Q_sun*AP*sind(beta);%%%%%%%%%%%%%%%port face
 else
   Q_solarP(theta)=0;
end
Q_solarS(theta)=0;
elseif beta<0
    Q_solarS=Q_sun*AS*sind(beta);%%%%%%%%%%%%%%%starboard face
 if (theta < theta_entry||theta> theta_exit)  
Q_solarS(theta)=Q_sun*AS*sind(beta);
 else
   Q_solarS(theta)=0;
end
Q_solarP(theta)=0 ;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%albedo%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_albedoZ =zeros(1,360);
theta=1:360
Q_albedoN =Q_sun*rho*VF_nadir*alpha*AN*cosd(beta)*cosd(theta);
for  theta= 90:270
    Q_albedoN(theta)=0.0;
end
theta=1:360
Q_albedoF(theta)=Q_sun*rho*VF_walls*alpha*AF*cosd(beta)*cosd(theta);
for  theta= 90:270
    Q_albedoF(theta)=0.0;
end
theta=1:360
Q_albedoA(theta)=Q_sun*rho*VF_walls*alpha*AA*cosd(beta)*cosd(theta);
for  theta= 90:270
    Q_albedoA(theta)=0.0;
end
theta=1:360
Q_albedoP(theta)=Q_sun*rho*VF_walls*alpha*AP*cosd(beta)*cosd(theta);
for  theta= 90:270
    Q_albedoP(theta)=0.0;
end
theta=1:360
Q_albedoS(theta)=Q_sun*rho*VF_walls*alpha*AP*cosd(beta)*cosd(theta);
for  theta= 90:270
    Q_albedoS(theta)=0.0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%IR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Q_IRZ =zeros(1,360);
Q_IRN= Q_ir*VF_nadir*emi*AN*ones(1,360);
Q_IRF= Q_ir*VF_walls*emi*AF*ones(1,360);
Q_IRA=Q_ir*VF_walls*emi*AA*ones(1,360);
Q_IRP=Q_ir*VF_walls*emi*AP*ones(1,360);
Q_IRS=Q_ir*VF_walls*emi*AS*ones(1,360);

%%%%%%%%%%%%%%%%%%%%%%%%%net flux%%%%%%%%%%%%%%%%%%%%%%%%
theta=1:360
zenith=(Q_solarZ+Q_albedoZ+Q_IRZ); 
nadir=Q_solarN+Q_albedoN+Q_IRN; 
forward=(Q_solarF+Q_albedoF+Q_IRF);
aft=(Q_solarA+Q_albedoA+Q_IRA); 
port=(Q_solarP+Q_albedoP+Q_IRP);
starboard=(Q_solarS+Q_albedoS+Q_IRS);
%%%%%%%%%%%%%%%%%%%%%%%%%Graphs for fluxs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot (theta,Q_solarZ,theta,Q_albedoZ,theta,Q_IRZ);
hleg1 = legend('zenith','Location','NorthEastoutside');
xlabel('Orbital Angle')
ylabel('Environmental Heat Flux (W)')
axis([0 365 -1 50])
figure
plot (theta,Q_solarN,theta,Q_albedoN,theta,Q_IRN);
hleg2 = legend('nadir','Location','NorthEastoutside');
xlabel('Orbital Angle')
ylabel('Environmental Heat Flux (W)')
axis([0 365 -1 50])
figure
plot (theta,Q_solarF,theta,Q_albedoF,theta,Q_IRF);
hleg3 = legend('forward','Location','NorthEastoutside');
xlabel('Orbital Angle')
ylabel('Environmental Heat Flux (W)')
axis([0 365 -1 50])
figure
plot (theta,Q_solarA,theta,Q_albedoA,theta,Q_IRA);
hleg4 = legend('aft','Location','NorthEastoutside');
xlabel('Orbital Angle')
ylabel('Environmental Heat Flux (W)')
axis([0 365 -1 50])
figure
plot (theta,Q_solarP,theta,Q_albedoP,theta,Q_IRP);
hleg5 = legend('port','Location','NorthEastoutside');
xlabel('Orbital Angle')
ylabel('Environmental Heat Flux (W)')
axis([0 365 -1 50])
figure 
plot (theta,Q_solarS,theta,Q_albedoS,theta,Q_IRS);
hleg6 = legend('starboard','Location','NorthEastoutside');
xlabel('Orbital Angle')
ylabel('Environmental Heat Flux (W)')
axis([0 365 -1 50])

figure 
plot(theta,zenith,theta,nadir,theta,forward,theta,aft,theta,port,theta,starboard);
hleg7 = legend('zenith','nadir','forward','aft','port','starboard','Location','NorthEastoutside');
title('Thermal Heat Inputs of 2U sat')
xlabel('Orbital Angle')
ylabel('Environmental Heat Flux (W)')
axis([0 365 -1 50])
%%%%%%%%%%%%%%%%%%%%%%%%temperature calculation for each surface%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta=1:361;
TZ= ((zenith/(emi*AZ*sigma)).^(1/4))-273;
TN= ((nadir/(emi*AN*sigma)).^(1/4))-273;
TF= ((forward/(emi*AF*sigma)).^(1/4))-273;
TA= ((aft/(emi*AA*sigma)).^(1/4))-273;
TP= ((port/(emi*AP*sigma)).^(1/4))-273;
TS= ((starboard/(emi*AS*sigma)).^(1/4))-273;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%graphs for variation of temperature in each face
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
plot(theta,TZ,theta,TN,theta,TF,theta,TA,theta,TP,theta,TS);
hleg8 = legend('zenith','nadir','forward','aft','port','starboard','Location','NorthEastoutside');
title('temperature of 2Usat')
xlabel('Orbital Angle')
ylabel('Environmental Heat Flux (W)')
axis([0 365 -1 50])