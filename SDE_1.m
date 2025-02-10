% Stochastic Differential Equations (reproducing figure 3.1)
% Kevin Roberts
% February 2025

clear all
close all
clc

%%
randn('state',12);

dt=0.001;
finaltime=1; 

n=finaltime/dt;
numberofrealisations=6;
X=zeros(numberofrealisations,n);
time=[0:dt:finaltime];

for i=1:numberofrealisations
    dX = sqrt(dt)*randn(1,n);   
    X(i,:)=cumsum(dX(1,:));             
end

figure(1);
set(gca,'Fontsize',18);
plot(time,[0,X(1,:)],'r');   
hold on;
plot(time,[0,X(2,:)],'g');   
plot(time,[0,X(3,:)],'c');   
plot(time,[0,X(4,:)],'b');   
plot(time,[0,X(5,:)],'m');   
plot(time,[0,X(6,:)],'Color',[0.9 0.9 0]);  
line([0 1],[0 0],'Linewidth',2,'Linestyle','--','Color','k');
xlabel('$t$','interpreter','latex');
ylabel('$X$','interpreter','latex');
set(get(gca,'ylabel'),'rotation',0);
axis([0 1 -1.8 1.8]);
set(gca,'Fontsize',18);

%%

randn('state',100);

dt=0.001;
finaltime=5; 

n=finaltime/dt;
numberofrealisations=6;
X=zeros(numberofrealisations,n);
time=[0:dt:finaltime];

for i=1:numberofrealisations
    dX = dt*ones(1,n)+sqrt(dt)*randn(1,n);   
    X(i,:)=cumsum(dX(1,:));             
end

figure(1);
set(gca,'Fontsize',18);
plot(time,[0,X(1,:)],'c');   
hold on;
plot(time,[0,X(2,:)],'m');   
plot(time,[0,X(3,:)],'g');   
plot(time,[0,X(4,:)],'Color',[0.9 0.9 0]);   
plot(time,[0,X(5,:)],'b');   
plot(time,[0,X(6,:)],'r');
line([0 5],[0 5],'Linewidth',2,'Linestyle','--','Color','k');
xlabel('$t$','interpreter','latex');
ylabel('$X$','interpreter','latex');
set(get(gca,'ylabel'),'rotation',0);
axis([0 5 -1.2 7]);
set(gca,'Fontsize',18);

%%

randn('state',5);

dt=0.1;
finaltime=10*60;  % 10 minutes
D=0.0001;         % 10^{-4} mm^2 sec^{-1}

n=finaltime/dt;
numberofrealisations=6;
X=zeros(numberofrealisations,n);
Y=zeros(numberofrealisations,n);

for i=1:numberofrealisations
    dX = sqrt(2*D*dt)*randn(1,n);   
    X(i,:)=cumsum(dX(1,:));             
    dY = sqrt(2*D*dt)*randn(1,n);   
    Y(i,:) = cumsum(dY(1,:));             
end;

Xplot(:,:)=X(:,5:5:n);
Yplot(:,:)=Y(:,5:5:n);

figure(1);
set(gca,'Fontsize',18);
plot([0,Xplot(4,:)],[0,Yplot(4,:)],'Color',[0.9 0.9 0]);   
hold on;
plot([0,Xplot(1,:)],[0,Yplot(1,:)],'m');   
plot([0,Xplot(2,:)],[0,Yplot(2,:)],'c');   
plot([0,Xplot(3,:)],[0,Yplot(3,:)],'g');   
plot([0,Xplot(5,:)],[0,Yplot(5,:)],'b');   
plot([0,Xplot(6,:)],[0,Yplot(6,:)],'r');  
line([-2.5 2.5],[0 0],'Color','k');
line([0 0],[-2.5 2.5],'Color','k');
plot([X(1,n)],[Y(1,n)],'ko','Linewidth',6);   
plot([X(2,n)],[Y(2,n)],'ko','Linewidth',6);   
plot([X(3,n)],[Y(3,n)],'ko','Linewidth',6);   
plot([X(4,n)],[Y(4,n)],'ko','Linewidth',6);   
plot([X(5,n)],[Y(5,n)],'ko','Linewidth',6);   
plot([X(6,n)],[Y(6,n)],'ko','Linewidth',6);  
xlabel('$x$ [mm]','interpreter','latex');
ylabel('$y$ [mm]','interpreter','latex');
set(gca,'XTick',[-0.6 -0.4 -0.2 0 0.2 0.4 0.6]);
set(gca,'YTick',[-0.6 -0.4 -0.2 0 0.2 0.4 0.6]);
axis([-0.7 0.7 -0.7 0.7]);
set(gca,'Fontsize',18);

%%

randn('state',20);

k1=0.001;
k2=0.75;
k3=165;
k4=10000;
k5=200;

dt=0.01;
finaltime=800; 
n=finaltime/dt;

X=zeros(n+1,1);
time=[0:dt:finaltime];

for i=1:n
    X(i+1) = X(i)+dt*(-k1*X(i)*X(i)*X(i)+k2*X(i)*X(i)-k3*X(i)+k4)+k5*sqrt(dt)*randn(1,1);   
end

[t,z] = ode15s(@myode,[0 5],[0]);

figure(1);
set(gca,'Fontsize',18);
line([5000 5000],[0 0],'Color','b','Linewidth',4);
hold on;
line([5000 5000],[100 100],'Color','r','Linewidth',4);
plot(time,X,'b');   
hold on;
plot(t,z,'r','Linewidth',3);
line([5 800],[100 100],'Color','r','Linewidth',3);
xlabel('$t$','interpreter','latex');
ylabel('$X$','interpreter','latex');
set(get(gca,'ylabel'),'rotation',0);
hh=legend('SDE','ODE');
set(hh,'interpreter','latex','Fontsize',18,'location','northwest');
axis([0 800 0 500]);
box on;
set(gca,'Fontsize',18);

function dydt = myode(t,z)
k1=0.001;
k2=0.75;
k3=165;
k4=10000;
dydt = [-k1*z(1)*z(1)*z(1)+k2*z(1)*z(1)-k3*z(1)+k4];
end