% Stochastic vs Deterministic modeling (reproducing figure 2.2)
% Kevin Roberts
% February 2025

clear all
close all
clc

rand('state',44);

k1=0.00025;
k2=0.18;
k3=37.5;
k4=2200;

X=0;
time=0;
timeSSA=0;
kk=0;

while (time<100)
   timeSSAprev=timeSSA; 
   timefinSSA=500;
   kk=kk+1;
   while (timeSSA<timeSSAprev+timefinSSA)
      rr=rand(2,1);
      timeSSA=timeSSA+1;
      a0=k1*X*(X-1)*(X-2)+k2*X*(X-1)+k3*X+k4;
      tau=(1/a0)*log(1/rr(1));
      if (rr(2)*a0<(k2*X*(X-1)+k4))
          X=X+1;
      else
          X=X-1;
      end
      time=time+tau;
   end;     
XX(kk)=X;
tt(kk)=time;
end

[t,z] = ode45(@myode,[0 100],[0]);

figure(1);
set(gca,'Fontsize',18);
line([5000 5000],[0 0],'Color','b','Linewidth',4);
hold on;
line([5000 5000],[0 0],'Color','r','Linewidth',4);
h=plot(tt,XX);
set(h,'Color','b','Linewidth',1);
plot(t,z,'r','Linewidth',3);
xlabel('time [min]','interpreter','latex');
ylabel('number of $A$ molecules','interpreter','latex');
axis([0 100 0 550]);
box on;
set(gca,'Fontsize',18);


function dydt = myode(t,z)
    k1=0.00025;
    k2=0.18;
    k3=37.5;
    k4=2200;
    dydt = [-k1*z(1)*z(1)*z(1)+k2*z(1)*z(1)-k3*z(1)+k4];
end