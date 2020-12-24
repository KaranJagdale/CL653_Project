%Written in Matlab
%please uncomment relevant part of code to obtain various plots
%avoid plotting all the plots at once to avoide plots being override
%Currently the code plots estimation error with bound and NEES
clf
randn('state',4)   %4 is by far the best seed.
%disp('Karan')
%Definitions
Ts = 0.1; N = 1000;

sigma1 = sqrt(0.2); sigma2 = 0.007*pi/180; sigmar = 50; sigmatheta = 0.1*pi/180; sigmaphi = 0.1*pi/180; 
%unit of sigmatheta and sigmaphi is radian

Qd = diag([0 sigma1^2 0 sigma1^2 0 sigma1^2 sigma2^2]);
R = diag([sigmar^2 sigmatheta^2 sigmaphi^2]);

x0 = [1000; 0; 2650; 150; 200; 0; 3*pi/180];  %turn rate is in rad/sec
xeststart = [1050; 10; 2650*1.05; 150*1.05; 200*1.05; 10; 3*1.05*pi/180];  %tis value will be used for the starting estimate for all filters
%Definitions

%True State generation
Xt = x0;   %Xt is used to store true states theta is stored as radian
for i = 0:1:N-1
    [t,X] = ode45(@(t,x)dyn(t,x), [i*Ts, (i+1)*Ts], Xt(:,size(Xt,2)));
    prossnoise = mvnrnd([0;0;0;0;0;0;0],Qd,1);
%     if i<3
%         disp(prossnoise')
%         disp((X(size(X,1),:)'))
%     end
    XtR = X(size(X,1),:)' + prossnoise';
    Xt = [Xt XtR];
    %Xt = [Xt X(size(X,1),:)'];
end



%Measurement generations, Va are the actual measurements
Va = [];
for i = 1:1:N+1
    xc = Xt(:,i);
    v = [sqrt(xc(1)^2 + xc(3)^2 + xc(5)^2);
        atan2(xc(3),xc(1));
        atan2(xc(5),sqrt(xc(1)^2 + xc(3)^2))];
    if v(2)<0
        v(2) = v(2) + 2*pi;
    end
    va = v + mvnrnd([0;0;0],R,1)';
    Va = [Va va];
end

x500 = Xt(:,500);     %this is the point if linearization
xdot500 = dyn(1,Xt(:,500)); % this term is equal to B

%uncomment following plots to obtain measurements
% figure(2)
% plot(TKF, Va(2,:)*180/pi)
% title('y2 (theta) vs time')
% xlabel('time (S)')
%legend('y2','y3')

%disp(x500)
% figure(3)
% plot(TKF,Va(1,:))
% title('y1 (range) vs time')
% xlabel('time (S)')
% 
% 
% figure(4)
% plot(TKF,Va(3,:)*180/pi)
% title('y3 (phi) vs time')
% xlabel('time (S)')

%After descretizing, we get
A = [0 1 0 0 0 0 0;
    0 0 0 -0.0524 0 0 129.5;
    0 0 0 1 0 0 0;
    0 0.0524 0 0 0 0 -75.7;
    0 0 0 0 0 1 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];


Phi = expm(A*Ts); %for discretised system
[tgamma, Gamma] = ode45(@(t,Gamma)expm(A*t)*xdot500, [0, 0.1], [0 0 0 0 0 0 0]');
Gammatilde = Gamma(size(Gamma,1),:); %Gain of control for discretised system
%disp(Gammatilde)

%Now applying kalman filter
%xcap0_0 = x0*1.05-x500;  %remember that x is a deviation variable0
xcap0_0 = xeststart-x500;
pcap0_0 = eye(7)*0.5;
C = Clin(x500);         %linearizing measurements
KFstates = x500 + xcap0_0;
KFvar = pcap0_0;
TKF = 0;
xcap = xcap0_0;
pcap = pcap0_0;
KFE = [];
KFSRP = [];
KFSRU = [];
KFesterr = [];  %estimation error of the filter
KFestbound = [];   %3*std values for plotting
for i = 1:1:size(Va,2)-1
    
    devY = Va(:,i) - Va(:,500);
    xpred = Phi*xcap + Gammatilde';  %u = 1
    ppred = Phi*pcap*Phi' + Qd;
    KFSRP = [KFSRP max(abs(eig(ppred)))];   %predicted covarience spectral radius
    if det(C*ppred*C' + R) == 0
        disp(C*ppred*C' + R)
        break
    end
%     disp(det(C*ppred*C' + R))
%     disp(C*ppred*C' + R)
    KFG = ppred*C'/((C*ppred*C' + R));  %Kalman gain
    KFe = devY - C*xpred;
    %KFe = Va(:,i) - mes(xpred + x500);
    xestKF = xpred + KFG*KFe;
    xcap = xestKF;
    pestKF = (eye(7) - KFG*C)*ppred;
    KFSRU = [KFSRU max(abs(eig(pestKF)))];  %updated covariance spectral radius
    pcap = pestKF;
    eststateKF = xcap + x500;
    KFstates = [KFstates eststateKF];
    KFvar = [KFvar pcap];
    TKF = [TKF i*Ts];
    KFE = [KFE KFe];
    kfesterr = Xt(:,i) - eststateKF;
    KFesterr = [KFesterr kfesterr];
    kfestbound = diag(pcap);
    KFestbound = [KFestbound kfestbound];
end

% figure(1)
% plotrange = 1000;
% plotstate = 7;
% plot(TKF(1:plotrange),Xt(plotstate,1:plotrange))
% hold on
% plot(TKF(1:plotrange),KFstates(plotstate,1:plotrange),'+-')
% legend('x2','x2KF')
% xlabel('Time (Seconds)') 

%------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%Now applying EKF
%xcap0_0E = x0*1.05;
xcap0_0E = xeststart;
pcap0_0E = eye(7)*1;
%pcap0_0E = diag([30 1 2 1 2 1 0.0001])*5;
EKFstates = xcap0_0E;
EKFvar = pcap0_0E;
TEKF = 0;
xcapE = xcap0_0E;
pcapE = pcap0_0E;
EKFE = [];
EKFSRP = [];
EKFSRU = [];
EKFesterr = [];
EKFestbound = [];
for i = 1:1:size(Va,2)-1
    [tEKF,XpredE] = ode45(@(t,x)dyn(t,x), [i*Ts (i+1)*Ts], xcapE);  %integration for the prediction step
    xpredE = XpredE(size(XpredE,1),:)';
    A = jacobian(xcapE);
    PhiE = expm(A*Ts);
    BE = dyn(1,xcapE);  %passting t=1 as it doesn not matter what t is passed
%     [tgammaE, GammaE] = ode45(@(t,Gamma)expm(A*t)*BE, [0, Ts], [0 0 0 0 0 0 0]');
%     GammatildeE = GammaE(size(GammaE,1),:);
%     if i < 5
%         disp(Va(:,i) - mes(xcapE))
%         disp('Karan')
%         disp(Va(:,i) - mes(xpredE))
%     end
    PpredE = PhiE*pcapE*PhiE' + Qd;
    EKFSRP = [EKFSRP max(abs(eig(PpredE)))];
    CE = Clin(xpredE);
    EKFG = PpredE*CE'/((CE*PpredE*CE' + R));
    EKFe = Va(:,i) - mes(xpredE);
    xestEKF = xpredE + EKFG*EKFe;
    xcapE = xestEKF;
    pestEKF = (eye(7) - EKFG*CE)*PpredE;
    EKFSRU = [EKFSRU max(abs(eig(pestEKF)))];
    pcapE = pestEKF;
    EKFstates = [EKFstates xcapE];
    EKFvar = [EKFvar pcapE];
    EKFE = [EKFE EKFe];
    ekfesterr = Xt(:,i) - xcapE;
    EKFesterr = [EKFesterr ekfesterr];
    ekfestbound = diag(pcapE);
    EKFestbound = [EKFestbound ekfestbound];
end


% figure(3)
% plotrange = 1000;
% plotstate = 7;
% plot(TKF(1:plotrange),Xt(plotstate,1:plotrange))
% hold on
% plot(TKF(1:plotrange),EKFstates(plotstate,1:plotrange),'+-')
% legend('x2','x2EKF')
% xlabel('Time (Seconds)') 

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%Now applying UKF

%xcap0_0U = x0*1.05;
xcap0_0U = xeststart;
xcapU = xcap0_0U;
pcap0_0U = eye(7)*1;
pcapU = pcap0_0U;
ns = 7; kappa = 1;   %ns is number of states
rho = sqrt(ns + kappa);
UKFstates = xcapU;
UKFvar = pcapU;
UKFE = [];
UKFSRP = [];
UKFSRU = [];
UKFesterr = [];
UKFestbound = [];
for i = 1:1:size(Va,2)-1
    %First, the prediction step-----------------------------------------------
    Apred = sqrtm(pcapU);
    chipred = zeros(ns,2*ns+1);  %this contains all the smaples for predict
    chipred(:,1) = xcapU;
    wpred = zeros(1,2*ns+1);    %contains weights of the samples
    wpred(1) = kappa/(ns+kappa);
    
    for j = 1:1:ns
        chipred(:,j+1) = xcapU + rho*Apred(:,j);
        wpred(j+1) = 1/(2*(ns + kappa));
    end
    
    for k = 1:1:ns
        chipred(:,k+1+ns) = xcapU - rho*Apred(:,k);
        wpred(k+ns+1) = 1/(2*(ns + kappa));
    end
    
    %Propogating the samples through the model
    chipredprop = zeros(ns,2*ns+1);
    for l = 1:1:2*ns+1
        chipredprop(:,l) = modelprop(chipred(:,l),Ts*(i-1),Ts*i);   %Confused about which i to use
    end
    xpredU = chipredprop*wpred';    %x(k|k-1) using the propogated samples
    xpredU = zeros(ns,1);
    for n= 1:1:2*ns+1
        xpredU = xpredU + wpred(n)*chipredprop(:,n);
    end
    ppredUs = zeros(ns);              %predicted covarience from samples
    for m = 1:1:2*ns+1
        ppredUs = ppredUs + wpred(m)*(chipredprop(:,m)-xpredU)*(chipredprop(:,m)-xpredU)';
    end
    ppredU = ppredUs + Qd;  %actual predicted covarience
    UKFSRP = [UKFSRP max(abs(eig(ppredU)))];
    %Now, the update steps-------------------------------------------------
    Aup = sqrtm(ppredU);
    
    chiup = zeros(ns,2*ns+1);  %this contains all the smaples for update
    chiup(:,1) = xpredU;
    wup = zeros(1,2*ns+1);    %contains weights of the samples
    wup(1) = kappa/(ns+kappa);
    
    for j = 1:1:ns
        chiup(:,j+1) = xpredU + rho*Aup(:,j);
        wup(j+1) = 1/(2*(ns + kappa));
    end
    
    for k = 1:1:ns
        chiup(:,k+1+ns) = xpredU - rho*Aup(:,k);
        wup(k+ns+1) = 1/(2*(ns + kappa));
    end
    
    Ys = zeros(3,2*ns+1);   %measurement samples obtained by predicted states
    for l = 1:1:2*ns+1
        Ys(:,l) = mes(chiup(:,l));
    end
    ypredU = Ys*wup';   %mean of the measurement samples
    ypredU = zeros(3,1);
    for n = 1:1:2*ns+1
        ypredU = ypredU + wup(n)*Ys(:,n);
    end
    %Now computing Pee and Pepseps
    pmesmes = zeros(3);
    pmesstate = zeros(ns,3);
    for m = 1:1:2*ns+1
        pmesstate = pmesstate + wup(m)*(chipredprop(:,m) - xpredU)*(Ys(:,m)-ypredU)';
        pmesmes = pmesmes + wup(m)*(Ys(:,m)-ypredU)*(Ys(:,m)-ypredU)';
    end
    
    pmesmes = pmesmes + R;
    UKFG = pmesstate/pmesmes;
    UKFe = Va(:,i) - ypredU;
    
    xcapU = xpredU + UKFG*UKFe;
    pcapU = ppredU - UKFG*pmesmes*UKFG';
    UKFSRU = [UKFSRU max(abs(eig(pcapU)))];
    
    UKFstates = [UKFstates xcapU];
    UKFvar = [UKFvar pcapU];
    UKFE = [UKFE UKFe];
    ukfesterr = Xt(:,i) - xcapU;
    UKFesterr = [UKFesterr ukfesterr];
    ukfestbound = diag(pcapU);
    UKFestbound = [UKFestbound ukfestbound];
end

% figure(4)
% plotrangeS = 1;
% plotrangeE = 1000;
% plotstate = 7;
% plot(TKF(plotrangeS:plotrangeE),Xt(plotstate,plotrangeS:plotrangeE))
% hold on
% plot(TKF(plotrangeS:plotrangeE),UKFstates(plotstate,plotrangeS:plotrangeE),'+-')
% legend('x2','x2UKF')
% xlabel('Time (Seconds)') 

%Estimated states vs time----------------------------------
%Uncomment figure 5,6,7 to get state estimate and innovation plots
% figure(5)    
% subplot(2,1,1)
% plotrangeS = 1;
% plotrangeE = 1000;
% plotstate = 5;
% plot(TKF(plotrangeS:plotrangeE),Xt(plotstate,plotrangeS:plotrangeE))
% hold on
% plot(TKF(plotrangeS:plotrangeE),UKFstates(plotstate,plotrangeS:plotrangeE),'-')
% plot(TKF(plotrangeS:plotrangeE),KFstates(plotstate,plotrangeS:plotrangeE),'-')
% plot(TKF(plotrangeS:plotrangeE),EKFstates(plotstate,plotrangeS:plotrangeE),'-')
% legend('true','UKF','KF', 'EKF')
% xlabel('Time (S)') 
% ylabel('x5 (m)')
% title('x5 state')
% grid on
% 
% subplot(2,1,2)
% plotrangeS = 1;
% plotrangeE = 200;
% plotstate = 5;
% plot(TKF(plotrangeS:plotrangeE),Xt(plotstate,plotrangeS:plotrangeE))
% hold on
% plot(TKF(plotrangeS:plotrangeE),UKFstates(plotstate,plotrangeS:plotrangeE),'+-')
% plot(TKF(plotrangeS:plotrangeE),KFstates(plotstate,plotrangeS:plotrangeE),'o-')
% plot(TKF(plotrangeS:plotrangeE),EKFstates(plotstate,plotrangeS:plotrangeE),'*-')
% legend('true','UKF','KF', 'EKF')
% xlabel('Time (S)') 
% ylabel('x5 (m)')
% title('x5 state (magnified)')
% grid on
% %-----------------------------------------------------------------
% 
% %another state------------------------------------------
% figure(6) 
% subplot(2,1,1)
% plotrangeS = 1;
% plotrangeE = 1000;
% plotstate = 3;
% plot(TKF(plotrangeS:plotrangeE),Xt(plotstate,plotrangeS:plotrangeE))
% hold on
% plot(TKF(plotrangeS:plotrangeE),UKFstates(plotstate,plotrangeS:plotrangeE),'-')
% plot(TKF(plotrangeS:plotrangeE),KFstates(plotstate,plotrangeS:plotrangeE),'-')
% plot(TKF(plotrangeS:plotrangeE),EKFstates(plotstate,plotrangeS:plotrangeE),'-')
% legend('true','UKF','KF', 'EKF')
% xlabel('Time (S)') 
% ylabel('x3 (m)')
% title('x3 state')
% grid on
% 
% subplot(2,1,2)
% plotrangeS = 1;
% plotrangeE = 150;
% plotstate = 3;
% plot(TKF(plotrangeS:plotrangeE),Xt(plotstate,plotrangeS:plotrangeE))
% hold on
% plot(TKF(plotrangeS:plotrangeE),UKFstates(plotstate,plotrangeS:plotrangeE),'+-')
% plot(TKF(plotrangeS:plotrangeE),KFstates(plotstate,plotrangeS:plotrangeE),'o-')
% plot(TKF(plotrangeS:plotrangeE),EKFstates(plotstate,plotrangeS:plotrangeE),'*-')
% legend('true','UKF','KF', 'EKF')
% xlabel('Time (S)') 
% ylabel('x3 (m)')
% title('x3 state (magnified)')
% grid on
% 
% %-------------------------------------------------------------------------
% 
% %Innovations with time------------------------------------------------
% figure(7)
% subplot(2,1,1)
% plotrangeS = 1;
% plotrangeE = 1000;
% plotstate = 3;
% plot(TKF(plotrangeS:plotrangeE),UKFE(plotstate,plotrangeS:plotrangeE),'+')
% hold on
% plot(TKF(plotrangeS:plotrangeE),KFE(plotstate,plotrangeS:plotrangeE),'o')
% plot(TKF(plotrangeS:plotrangeE),EKFE(plotstate,plotrangeS:plotrangeE),'*')
% legend('UKF','KF','EKF')
% xlabel('Time (S)') 
% ylabel('innovation')
% title('y3 innovation')
% 
% subplot(2,1,2)
% plotrangeS = 1;
% plotrangeE = 1000;
% plotstate = 2;
% plot(TKF(plotrangeS:plotrangeE),UKFE(plotstate,plotrangeS:plotrangeE),'+')
% hold on
% plot(TKF(plotrangeS:plotrangeE),KFE(plotstate,plotrangeS:plotrangeE),'o')
% plot(TKF(plotrangeS:plotrangeE),EKFE(plotstate,plotrangeS:plotrangeE),'*')
% legend('UKF','KF','EKF')
% xlabel('Time (S)') 
% ylabel('innovation')
% title('y2 innovation')

%Uncomment to plot spectral radius
% figure(2)
% plot(TKF(1:10:1000), KFSRP(1:10:1000), '*-')
% hold on
% plot(TKF(1:10:1000), KFSRU(1:10:1000), '+-')
% legend('Predicted', 'Updated')
% xlabel('Time (S)')
% ylabel('Spectral radius')
% title('KF covarience matrix spectral radius')
% 
% figure(3)
% plot(TKF(1:10:1000), EKFSRP(1:10:1000), '*-')
% hold on
% plot(TKF(1:10:1000), EKFSRU(1:10:1000), '+-')
% legend('Predicted', 'Updated')
% xlabel('Time (S)')
% ylabel('Spectral radius')
% title('EKF covarience matrix spectral radius')
% 
% figure(4)
% plot(TKF(1:10:1000), UKFSRP(1:10:1000), '*-')
% hold on
% plot(TKF(1:10:1000), UKFSRU(1:10:1000), '+-')
% legend('Predicted', 'Updated')
% xlabel('Time (S)')
% ylabel('Spectral radius')
% title('UKF covarience matrix spectral radius')

%plotting estimation error with bounds
figure(9)
subplot(2,1,1)
stateplot = 1;
errorbar(TKF(stateplot:20:1000),KFesterr(stateplot,1:20:1000),3*sqrt(KFestbound(stateplot,1:20:1000)))
hold on
errorbar(TKF(stateplot:20:1000),EKFesterr(stateplot,1:20:1000),3*sqrt(EKFestbound(stateplot,1:20:1000)))
errorbar(TKF(stateplot:20:1000),UKFesterr(stateplot,1:20:1000),3*sqrt(UKFestbound(stateplot,1:20:1000)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Estimation error for state 1')
 

subplot(2,1,2)
stateplot = 1;
range = 1:6:300;
errorbar(TKF(range),KFesterr(stateplot,range),3*sqrt(KFestbound(stateplot,range)))
hold on
errorbar(TKF(range),EKFesterr(stateplot,range),3*sqrt(EKFestbound(stateplot,range)))
errorbar(TKF(range),UKFesterr(stateplot,range),3*sqrt(UKFestbound(stateplot,range)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Magnified estimation error state 1')

figure(10)
subplot(2,1,1)
stateplot = 2;
errorbar(TKF(1:20:1000),KFesterr(stateplot,1:20:1000),3*sqrt(KFestbound(stateplot,1:20:1000)))
hold on
errorbar(TKF(1:20:1000),EKFesterr(stateplot,1:20:1000),3*sqrt(EKFestbound(stateplot,1:20:1000)))
errorbar(TKF(1:20:1000),UKFesterr(stateplot,1:20:1000),3*sqrt(UKFestbound(stateplot,1:20:1000)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Estimation error state 2')

subplot(2,1,2)
stateplot = 2;
range = 1:4:200;
errorbar(TKF(range),KFesterr(stateplot,range),3*sqrt(KFestbound(stateplot,range)))
hold on
errorbar(TKF(range),EKFesterr(stateplot,range),3*sqrt(EKFestbound(stateplot,range)))
errorbar(TKF(range),UKFesterr(stateplot,range),3*sqrt(UKFestbound(stateplot,range)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Magnified estimation error state 2')

figure(11)
subplot(2,1,1)
stateplot = 3;
errorbar(TKF(1:20:1000),KFesterr(stateplot,1:20:1000),3*sqrt(KFestbound(stateplot,1:20:1000)))
hold on
errorbar(TKF(1:20:1000),EKFesterr(stateplot,1:20:1000),3*sqrt(EKFestbound(stateplot,1:20:1000)))
errorbar(TKF(1:20:1000),UKFesterr(stateplot,1:20:1000),3*sqrt(UKFestbound(stateplot,1:20:1000)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Estimation error state 3')

subplot(2,1,2)
stateplot = 3;
range = 1:8:400;
errorbar(TKF(range),KFesterr(stateplot,range),3*sqrt(KFestbound(stateplot,range)))
hold on
errorbar(TKF(range),EKFesterr(stateplot,range),3*sqrt(EKFestbound(stateplot,range)))
errorbar(TKF(range),UKFesterr(stateplot,range),3*sqrt(UKFestbound(stateplot,range)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Magnified estimation error state 3')


figure(12)
subplot(2,1,1)
stateplot = 4;
errorbar(TKF(1:20:1000),KFesterr(stateplot,1:20:1000),3*sqrt(KFestbound(stateplot,1:20:1000)))
hold on
errorbar(TKF(1:20:1000),EKFesterr(stateplot,1:20:1000),3*sqrt(EKFestbound(stateplot,1:20:1000)))
errorbar(TKF(1:20:1000),UKFesterr(stateplot,1:20:1000),3*sqrt(UKFestbound(stateplot,1:20:1000)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Estimation error state 4')

subplot(2,1,2)
stateplot = 4;
range = 1:8:400;
errorbar(TKF(range),KFesterr(stateplot,range),3*sqrt(KFestbound(stateplot,range)))
hold on
errorbar(TKF(range),EKFesterr(stateplot,range),3*sqrt(EKFestbound(stateplot,range)))
errorbar(TKF(range),UKFesterr(stateplot,range),3*sqrt(UKFestbound(stateplot,range)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Magnified estimation error state 4')


figure(13)
subplot(2,1,1)
stateplot = 5;
errorbar(TKF(1:20:1000),KFesterr(stateplot,1:20:1000),3*sqrt(KFestbound(stateplot,1:20:1000)))
hold on
errorbar(TKF(1:20:1000),EKFesterr(stateplot,1:20:1000),3*sqrt(EKFestbound(stateplot,1:20:1000)))
errorbar(TKF(1:20:1000),UKFesterr(stateplot,1:20:1000),3*sqrt(UKFestbound(stateplot,1:20:1000)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Estimation error state 5')

subplot(2,1,2)
stateplot = 5;
range = 1:8:400;
errorbar(TKF(range),KFesterr(stateplot,range),3*sqrt(KFestbound(stateplot,range)))
hold on
errorbar(TKF(range),EKFesterr(stateplot,range),3*sqrt(EKFestbound(stateplot,range)))
errorbar(TKF(range),UKFesterr(stateplot,range),3*sqrt(UKFestbound(stateplot,range)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Magnified estimation error state 5')


figure(14)
subplot(2,1,1)
stateplot = 6;
errorbar(TKF(1:20:1000),KFesterr(stateplot,1:20:1000),3*sqrt(KFestbound(stateplot,1:20:1000)))
hold on
errorbar(TKF(1:20:1000),EKFesterr(stateplot,1:20:1000),3*sqrt(EKFestbound(stateplot,1:20:1000)))
errorbar(TKF(1:20:1000),UKFesterr(stateplot,1:20:1000),3*sqrt(UKFestbound(stateplot,1:20:1000)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Estimation error state 6')

subplot(2,1,2)
stateplot = 6;
range = 1:8:400;
errorbar(TKF(range),KFesterr(stateplot,range),3*sqrt(KFestbound(stateplot,range)))
hold on
errorbar(TKF(range),EKFesterr(stateplot,range),3*sqrt(EKFestbound(stateplot,range)))
errorbar(TKF(range),UKFesterr(stateplot,range),3*sqrt(UKFestbound(stateplot,range)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Magnified estimation error state 6')


figure(15)
subplot(2,1,1)
stateplot = 7;
errorbar(TKF(1:20:1000),KFesterr(stateplot,1:20:1000),3*sqrt(KFestbound(stateplot,1:20:1000)))
hold on
errorbar(TKF(1:20:1000),EKFesterr(stateplot,1:20:1000),3*sqrt(EKFestbound(stateplot,1:20:1000)))
errorbar(TKF(1:20:1000),UKFesterr(stateplot,1:20:1000),3*sqrt(UKFestbound(stateplot,1:20:1000)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Estimation error state 7')

subplot(2,1,2)
stateplot = 7;
range = 1:2:100;
errorbar(TKF(range),KFesterr(stateplot,range),3*sqrt(KFestbound(stateplot,range)))
hold on
errorbar(TKF(range),EKFesterr(stateplot,range),3*sqrt(EKFestbound(stateplot,range)))
errorbar(TKF(range),UKFesterr(stateplot,range),3*sqrt(UKFestbound(stateplot,range)))
legend('KF', 'EKF', 'UKF')
xlabel('Time (S)')
title('Magnified estimation error state 7')


%Calculating mena and varience in each innovation state with time
KFinnomean = [mean(KFE(1,:));mean(KFE(2,:));mean(KFE(2,:))];
KFinnovar = [var(KFE(1,:)); var(KFE(2,:)); var(KFE(3,:))];

EKFinnomean = [mean(EKFE(1,:));mean(EKFE(2,:));mean(EKFE(2,:))];
EKFinnovar = [var(EKFE(1,:)); var(EKFE(2,:)); var(EKFE(3,:))];

UKFinnomean = [mean(UKFE(1,:));mean(UKFE(2,:));mean(UKFE(2,:))];
UKFinnovar = [var(UKFE(1,:)); var(UKFE(2,:)); var(UKFE(3,:))];

%Calculating RMSE
KFRMSE = mean((KFesterr.^2)');
EKFRMSE = mean((EKFesterr.^2)');
UKFRMSE = mean((UKFesterr.^2)');

%Computing beta
alfa = 0.05;   %significance level
zhai1 = chi2inv(alfa, ns); zhai2 = chi2inv(1-alfa, ns); 
KFBeta = [];  %vector which will contain values wrt time
KFbcount = 0;
for i = 1:1:N
    KFbeta = KFesterr(:,i)'/KFvar(:,7*i+1:7*i+7)*KFesterr(:,i);
    KFBeta = [KFBeta KFbeta];
    if KFbeta> zhai2 || KFbeta <zhai1
        KFbcount = KFbcount +1;
    end
end

EKFBeta = [];  %vector which will contain values wrt time
EKFbcount = 0;
for i = 1:1:N
    EKFbeta = EKFesterr(:,i)'/EKFvar(:,7*i+1:7*i+7)*EKFesterr(:,i);
    EKFBeta = [EKFBeta EKFbeta];
    if EKFbeta > zhai2 || EKFbeta<zhai1
        EKFbcount = EKFbcount + 1;
    end
end

UKFBeta = [];  %vector which will contain values wrt time
UKFbcount = 0;
for i = 1:1:N
    UKFbeta = UKFesterr(:,i)'/UKFvar(:,7*i+1:7*i+7)*UKFesterr(:,i);
    UKFBeta = [UKFBeta UKFbeta];
    if UKFbeta > zhai2 || UKFbeta <zhai1
        UKFbcount = UKFbcount +1;
    end
    
end


% disp(KFbcount/1000)
% disp(EKFbcount/1000)
% disp(UKFbcount/1000)
% figure(15)
% range = 1:1000;
% plot(TKF(range),KFBeta(range))
% hold on
% plot(TKF(range),EKFBeta(range))
% plot(TKF(range),UKFBeta(range))
% yline(zhai1)
% yline(zhai2)
% % plot(TKF(range),ones(1,size(TKF(range),1))*zhai1)
% % plot(TKF(range),ones(1,size(TKF(range),1))*zhai2)
% legend('KF','EKF','UKF','zhai1','zhai2')
% title('Normalised estimation error squared')
% xlabel('time (S)')



