% Comparing deterministic vs stochastic (reproducing figure 2.5)
% Kevin Roberts
% February 2025

clear all
close all
clc

rand('state',5);

k1=0.00004;
k2=50;
k3=10;
k4=25;


A=10;
B=10;

time=0;
timeSSA=0;
kk=0;

while (time<5000)
   timeSSAprev=timeSSA; 
   timefinSSA=500;
   kk=kk+1;
   while (timeSSA<timeSSAprev+timefinSSA)
      if (A>100) 
          timefinSSA=10;
      end
      rr=rand(2,1);
      timeSSA=timeSSA+1;
      a0=k1*A*(A-1)*B+k2+k3*A+k4;
      tau=(1/a0)*log(1/rr(1));
      ss=k1*A*(A-1)*B;     
      if (ss>rr(2)*a0)
          A=A+1;
          B=B-1;
      else
          ss=ss+k2;
          if (ss>rr(2)*a0)
             A=A+1;
          else
             ss=ss+k3*A;
             if (ss>rr(2)*a0)
                 A=A-1;
             else
                 ss=ss+k4;
                 if (ss>rr(2)*a0)
                    B=B+1;
                 end
             end
          end
      end
      time=time+tau;
   end     
AA(kk)=A;
BB(kk)=B;
tt(kk)=time;
end

[t,z] = ode45(@myode,[0 5000],[10; 10]);

figure(1);
set(gca,'Fontsize',18);
h=semilogy([100 100],[20000 20000]);
set(h,'Color','b','Linewidth',4);
hold on;
h=semilogy([100 100],[20000 20000]);
set(h,'Color','r','Linewidth',4);
h=semilogy(tt/60,AA);
set(h,'Color','b','Linewidth',2);
h=semilogy(t/60,z(:,1));
set(h,'Color','r','Linewidth',4);
xlabel('time [min]','interpreter','latex');
ylabel('number of $A$ molecules','interpreter','latex');
hh=legend('stochastic','deterministic');
set(hh,'interpreter','latex','Fontsize',18);
axis([0 82 1 10000]);
set(gca,'YTick',[1 10 100 1000 10000]);
set(gca,'Fontsize',18);

figure(2);
set(gca,'Fontsize',18);
plot(tt/60,AA,'g','Linewidth',4);
hold on;
plot([100 100],[20000 20000],'b--','Linewidth',4);
plot(tt/60,BB,'b--','Linewidth',2);
xlabel('time [min]','interpreter','latex');
ylabel('number of molecules','interpreter','latex');
hh=legend('$A(t)$','$B(t)$');
set(hh,'interpreter','latex','Fontsize',18);
axis([0 82 1 10000]);
set(gca,'Fontsize',18);

figure(3);
set(gca,'Fontsize',18);
h=semilogx([20000 20000],[20000 20000]);
set(h,'Color','b','Linewidth',4);
hold on;
h=semilogx([20000 20000],[20000 20000]);
set(h,'Color','r','Linewidth',4);
xval1=[2:0.15:24.8 25:1:100 110:10:1000 1100:100:6200];
xval2=[2:0.15:24.8 25:1:160];
ynul1=(k3*xval1-k2)./(k1*xval1.*xval1);
ynul2=k4./(k1*xval2.*xval2);
semilogx(AA,BB,'b');
h=semilogx(z(1,1),z(1,2));
set(h,'Color','r','Linewidth',2);
xlabel('number of $A$ molecules','interpreter','latex');
ylabel('number of $B$ molecules','interpreter','latex');
h=semilogx(xval1,ynul1);
set(h,'Color','g','Linewidth',4);
h=semilogx(xval2,ynul2);
set(h,'Color','g','Linewidth',4);
h=semilogx([(k4+k2)/k3],[k4*k3*k3/(k1*(k4+k2)*(k4+k2))],'ro');
set(h,'Markersize',4,'Linewidth',6);
axis([1 10000 0 13000]);
h=semilogx(z(:,1),z(:,2));
set(h,'Color','r','Linewidth',2);
hh=legend('stochastic','deterministic');
set(hh,'interpreter','latex','Fontsize',18);
set(gca,'XTick',[1 10 100 1000 10000]);
set(gca,'Fontsize',18);

function dydt = myode(t,z)
k1=0.00004;
k2=50;
k3=10;
k4=25;
dydt =[k1*z(1)*z(1)*z(2)+k2-k3*z(1); -k1*z(1)*z(1)*z(2)+k4];
end