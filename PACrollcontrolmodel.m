clear all
clc
close all
global hdotr
% 6-DOF model for quadcopter for PID simulation with direct thrust control
% Roll control
% (Gofly HMB-1 model)
%global Va R climbrate

%inputfile
% thrust inputs T_1 to T_4
m=(2.000); % mass of UAV in Kg (UAV+payload)
rho=1.1; % air density kg/m^3

f1=0.00981; % scaling for gram to Newtons

g=9.81; % acceleration due to gravity  m/s^2
Jxx =0.039525793 ; 
Jyy =0.044924728 ;
Jzz = 0.058153811;
Jxz = 0.005078564;
 
 Jxyz = Jxx*Jzz-Jxz*Jxz;
 j1 = Jxz*(Jxx-Jyy+Jzz)/Jxyz;
 j2 = (Jzz*(Jzz-Jyy)+Jxz*Jxz)/Jxyz;
 j3 = Jzz/Jxyz;
 j4 = Jxz/Jxyz;
 j5 = (Jzz-Jxx)/Jyy;
 j6 = Jxz/Jyy;
 j7 = ((Jxx-Jyy)*Jxx+Jxz*Jxz)/Jxyz;
 j8 = Jxx/Jxyz;
 


 
 % Props T1(left ccw) and T2(right cw) (small rotors)
 % Props T3 (Bottom big rotor cw)
 % Props T4 (Top big rotor ccw)
 
  % origin --- big rotors

 xcg = 0.04; 
 ycg=0.0;
 zcg= 8/100; % down ward positive
 
 yr=12/100; % small right motor moment arm
 yl=12/100; % small left motor moment arm
 
 xr=(0.195-xcg); % small right motor moment arm
 xl=(0.195-xcg); % small left motor moment arm








dt=0.02;
tend=100;

tolval=0.05;
% initial conditions
u=0; v=0; w=0; p=0; q=0; r=0; 
phi=0; theta=0; psi=0;  pn=0; pe=0; h=0.0;

% d1 = diameter fpr prop 1 and 2 in inches
% d2 = diameter fpr prop 3 and 4 in inches
d1=5; d2=13;
% integrator initialconditions
sumh2=0.0;
sumh4=0.0;
sumphi2=0.0;
sumtheta2=0.0;
sumpsi2=0.0;

% thrust initial conditions

cfac=9.81/1000; % thrust grams to Newton

T1=200*cfac; T2=200*cfac; T3=800*cfac; T4=800*cfac;



T1h=200*cfac; T2h=200*cfac; T3h=800*cfac; T4h=800*cfac;  % hover thrust

for i=1:tend/dt





%----------------- force inputs---------------------------


Ftx=0;

%-------------

Fty=0.0;

%-----------------

Ftz=-(T1+T2+T3+T4);

%---------------------

% torque calculation

c1=-T1;
c2=-T2;
c3=-T3;
c4=-T4;


if(T1<=40*cfac)
    rpm1=(7120)/40;
else
    rpm1=-0.037153619683*T1*T1*(1/cfac)*(1/cfac)+54.289145434814*T1*(1/cfac)+5242.025646344952;
end 



    if(T2<=40*cfac)
    rpm2=(7120)/40;
else
    rpm2=-0.037153619683*T2*T2*(1/cfac)*(1/cfac)+54.289145434814*T2*(1/cfac)+5242.025646344952;
    end
    
    if(T3<=176*cfac)
    rpm3=(2385)/176;
else
    rpm3=-0.001537393831*T3*T3*(1/cfac)*(1/cfac)+5.279878995903*T3*(1/cfac)+1490.271882142766;
end 



    if(T4<=176*cfac)
    rpm4=(2385)/176;
else
    rpm4=-0.001537393831*T4*T4*(1/cfac)*(1/cfac)+5.279878995903*T4*(1/cfac)+1490.271882142766;
    end

    

            % power required from momentum theory

pow1=sqrt(T1*T1*T1/(2*rho*0.25*pi*(0.01*d1*2.54)*(0.01*d1*2.54)));
pow2=sqrt(T2*T2*T2/(2*rho*0.25*pi*(0.01*d1*2.54)*(0.01*d1*2.54)));

pow3=sqrt(T3*T3*T3/(2*rho*0.25*pi*(0.01*d2*2.54)*(0.01*d2*2.54)));
pow4=sqrt(T4*T4*T4/(2*rho*0.25*pi*(0.01*d2*2.54)*(0.01*d2*2.54)));



 
if(rpm1>0)
tau1=pow1*60/(2*pi*rpm1);
else
    tau1=0;
end
if(rpm2>0)
tau2=pow2*60/(2*pi*rpm2);
else
    tau2=0;
end
if(rpm3>0)
tau3=pow3*60/(2*pi*rpm3);
else
    tau3=0;
end
if(rpm4>0)
tau4=pow4*60/(2*pi*rpm4);
else
    tau4=0;
end
%------ moment inputs---------------
RM=(T1)*yl-(T2)*yr;
YM=tau1-tau2+tau4-tau3;
PM=-T1*xl-T2*xr+T3*xcg+T4*xcg;

% drag term cx,cy,cz   % std deviation 0.25
cx=1.05;cy=1.05;cz=1.0;

udot=r*v-q*w-g*sin(theta)+Ftx/m-cx*u;
vdot=p*w-r*u+g*cos(theta)*sin(phi)+Fty/m-cy*v;
wdot=q*u-p*v+g*cos(theta)*cos(phi)+Ftz/m-cz*w;

azb=Ftz/m-cz*w;

% damping terms cp cq cr
cp=10; cq=6; cr=7;

pdot=j1*p*q-j2*q*r+j3*RM+j4*YM-1*cp*p;
qdot=j5*p*r-j6*(p*p-r*r)+PM/Jyy-1*cq*q;
rdot=j7*p*q-j1*q*r+j4*RM+j8*YM-1*cr*r;


%-- attitude equations------------

 phidot  = (p + sin(phi)*tan(theta)*q + r*cos(phi)*tan(theta));
 thetadot = (q*cos(phi) - r*sin(phi));
 psidot = ((q*sin(phi) + r*cos(phi))/cos(theta));
 
 %------position equations------------
  pndot = ((cos(psi)*cos(theta))*u+ (cos(psi)*sin(theta)*sin(phi) - sin(psi)*cos(phi))*v+ (sin(psi)*sin(phi)+cos(psi)*sin(theta)*cos(phi))*w);
  pedot = ((sin(psi)*cos(theta))*u + (cos(psi)*cos(phi) + sin(psi)*sin(theta)*sin(phi))*v + (sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi))*w);
  hdot = (sin(theta)*u - cos(theta)*sin(phi)*v - cos(theta)*cos(phi)*w);
  

  
  % state propagation
  

%   
u=u+udot*dt;
v=v+vdot*dt;
w=w+wdot*dt;
p=p+pdot*dt;
q=q+qdot*dt;
r=r+rdot*dt;
phi=phi+phidot*dt;
theta=theta+thetadot*dt;
psi=psi+psidot*dt;
pn=pn+pndot*dt;
pe=pe+pedot*dt;
h=h+hdot*dt;

Va=sqrt(u*u+v*v+w*w);



%storing variables
tim(i)=i*dt;
ustore(i)=u;
vstore(i)=v;
wstore(i)=w;
pstore(i)=p*180/pi;
qstore(i)=q*180/pi;
rstore(i)=r*180/pi;
phistore(i)=phi*180/pi;
thetastore(i)=theta*180/pi;
psistore(i)=psi*180/pi;
pnstore(i)=pn;
pestore(i)=pe;
altitude(i)=h;

airspeedstore(i)=Va;

% controller design --- PIXHAWK control structure

%----------------------------------------------------------
% altitude in cm in pixhawk
% altitude control -- PID


% inputs---
href=0.0;  % maintain current altitude

hp1=1.0;
kph2=30;
kih2=0.0;
kdh2=0.0;
kph4=0.75*3;
kih4=1.50*3;
kdh4=0.0;
% 


%-------------------------

% error in position
h1=((href)-(h))*hp1;
errh1=href-h;
errh1s(i)=errh1;
if(i>1)
   errh1old=errh1s(i-1);
else
    errh1old=0.0;
end
derrh1=(errh1-errh1old)/dt;


if(abs(h1)>250)
    h1=250*sign(h1);
end


% error in velocity
h2=h1-(hdot);
h2s(i)=h2;
sumh2=sumh2+h2*dt;

if(i>1)
    h2old=h2s(i-1);
else
    h2old=h2;
end

dh2=(h2-h2old)/dt;

h3=kph2*h2+kih2*sumh2+kdh2*dh2;

% error in acceleration
% h3 into body frame
h3b=h3+9.81*(cos(theta)*cos(phi));
hzb=-azb;
h4=h3b-hzb; % error in body acceleration
h4s(i)=h4;

if(abs(h4)>250)
    h4=250*sign(h4);
end

sumh4=sumh4+h4*dt;

if(i>1)
    h4old=h4s(i-1);
else
    h4old=h4;
end

dh4=(h4-h4old)/dt;

h5=kph4*h4+kih4*sumh4+kdh4*dh4;

Tcount=h5;

%--------------------------------

% roll angle control -- PAC


%-- inputs----
if(tim(i)<=50)
phiref=5*pi/180; 
elseif((tim(i)>50)&&(tim(i)<=100))
    phiref=0*pi/180;
    elseif((tim(i)>100)&&(tim(i)<=150))
    phiref=5*pi/180;
else
    phiref=0*pi/180;
end




kphi1=4.5;
kpphi2=14*2/kphi1;
kiphi2=8*0.05/kphi1;
kdphi2=0/kphi1;

%-------------------




errphi=phiref-phi;

errphis(i)=errphi;

if(i>1)
    errphiold=errphis(i-1);
else
    errphiold=errphi;
end

derrphi=(errphi-errphiold)/dt;

% PAC computation

   inputpac=[errphi,derrphi,phi,i,phiref];
      phi1=1.0*pacgoflyroll(inputpac);

phi2=phi1-(p);

phi2s(i)=phi2;

sumphi2=sumphi2+phi2*dt;

if(i>1)
    phi2old=phi2s(i-1);
else
    phi2old=phi2;
end

dphi2=(phi2-phi2old)/dt;

phi3=kpphi2*phi2+kiphi2*sumphi2+kdphi2*dphi2;

Pcount=phi3; % total roll count (-4500 to +4500)



%------------------------------------------

% pitch angle control --- PID
%--inputs----
thetaref=0;


ktheta1=4.5;
kptheta2=14/ktheta1;
kitheta2=8*0.1/ktheta1;
kdtheta2=0.0/ktheta1;



%-------------
theta1=sign((thetaref-theta))*(abs(thetaref-theta))*ktheta1;

errtheta=thetaref-theta;

errthetas(i)=errtheta;

if(i>1)
    errthetaold=errthetas(i-1);
else
    errthetaold=errtheta;
end

derrtheta=(errtheta-errthetaold)/dt;


   
theta2=theta1-(q);
theta2s(i)=theta2;

sumtheta2=sumtheta2+theta2*dt;

if(i>1)
    theta2old=theta2s(i-1);
else
   theta2old=theta2;
end

dtheta2=(theta2-theta2old)/dt;

theta3=kptheta2*theta2+kitheta2*sumtheta2+kdtheta2*dtheta2;

Qcount=theta3; % total pitch count (-4500 to +4500)




%------------------------------------------

% Yaw angle control --- PID
%-- inputs----
psiref=0*pi/180;  % 





kpsi1=4.5;
kppsi2=14*5/kpsi1;
kipsi2=8*0.1/kpsi1;
kdpsi2=0/kpsi1;

%-------------------

psi1=(psiref-psi)*kpsi1;

errpsi=psiref-psi;

errpsis(i)=errpsi;

if(i>1)
    errpsiold=errpsis(i-1);
else
    errpsiold=errpsi;
end

derrpsi=(errpsi-errpsiold)/dt;


psi2=psi1-(r);

psi2s(i)=psi2;

sumpsi2=sumpsi2+psi2*dt;


if(i>1)
    psi2old=psi2s(i-1);
else
    psi2old=psi2;
end

dpsi2=(psi2-psi2old)/dt;

psi3=kppsi2*psi2+kipsi2*sumpsi2+kdpsi2*dpsi2;

Rcount=psi3; % total yaw count (-4500 to +4500)



%------------------------------------------

scale=sqrt(Tcount*Tcount+Pcount*Pcount+Qcount*Qcount+Rcount*Rcount);

Tc=(Tcount)/scale;

Pc=(Pcount)/scale;

Qc=(Qcount)/scale;

Rc=(Rcount)/scale;

% storing controls and references

thstore(i)=Tc;
rollstore(i)=Pc;
pitchstore(i)=Qc;
yawstore(i)=Rc;

thetarefstore(i)=thetaref*180/pi;
phirefstore(i)=phiref*180/pi;
psirefstore(i)=psiref*180/pi;
hrefstore(i)=href;

TM=[1,1,1,1; 1,-1,1,-1; -1,-1,1,1; 1,-1,-1,1];
Tv=inv(TM)*[Tc;Pc;Qc;Rc];


T1max=520*cfac*1.0;
T2max=520*cfac*1.0;
T3max=0.5*2370*cfac*1.0;
T4max=0.5*2370*cfac*1.0;

T1=T1h+Tv(1)*abs(h5)*1;
T2=T2h+Tv(2)*abs(h5)*1;
T3=T3h+Tv(3)*abs(h5)*1;
T4=T4h+Tv(4)*abs(h5)*1;
% 
if(T1>T1max)
    T1=T1max;
end
if(T2>T2max)
    T2=T2max;
end
if(T3>T3max)
    T3=T3max;
end
if(T4>T4max)
    T4=T4max;
end
    
    
end



figure(1)
subplot(3,1,1)
plot(tim,ustore)
grid on
xlabel('time (s)')
ylabel('u (m/s)')
subplot(3,1,2)
plot(tim,vstore)
grid on
xlabel('time (s)')
ylabel('v (m/s)')
subplot(3,1,3)
plot(tim,wstore)
grid on
xlabel('time (s)')
ylabel('w (m/s)')

figure(2)
subplot(3,1,1)
plot(tim,pstore)
grid on
xlabel('time (s)')
ylabel('p (deg/s)')
subplot(3,1,2)
plot(tim,qstore)
grid on
xlabel('time (s)')
ylabel('q (deg/s)')
subplot(3,1,3)
plot(tim,rstore)
grid on
xlabel('time (s)')
ylabel('r (deg/s)')



figure(3)
plot(tim,airspeedstore)
grid on
xlabel('time (s)')
ylabel('Va (m/s)')



figure(4)
plot(tim,thstore)
hold on
plot(tim,rollstore,'r')
plot(tim,pitchstore,'k')
plot(tim,yawstore,'g')
grid on
hold on
xlabel('time (s)')
ylabel('thrust, roll, pitch, yaw inputs')




figure(5)
plot(tim,phirefstore,'r')
hold on
plot(tim,phistore)
grid on
xlabel('time (s)')
ylabel('\phi (deg)')

