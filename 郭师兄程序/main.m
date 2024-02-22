% Rosen-Zener quantum battery
% 鍙傛暟鍙樺寲
tic;
lmada0=2;
m=length(lmada0);
N0=10;
p=length(N0);
E1n=zeros(m,p);
E2n=zeros(m,p);
E3n=zeros(m,p);
E4n=zeros(m,p);
E5n=zeros(m,p);
E6n=zeros(m,p);
E7n=zeros(m,p);
for j=1:1:m
    lmada=lmada0(j);
for k=1:1:p
    N=N0(k);
h=0.001;S=N/2;T=0.1*pi;v0=10;v1=20;v2=40;v3=60;
Sz=zeros(N+1,N+1);
H=zeros(N+1,N+1);
H0=zeros(N+1,N+1); 
H1=zeros(N+1,N+1); 
H11=zeros(N+1,N+1); 
H2=zeros(N+1,N+1); 
t0=0:h:T; 
n=length(t0);
E111=zeros(1,n);
E222=zeros(1,n);
E333=zeros(1,n);
E3=zeros(1,n);
P1=zeros(1,n);
P2=zeros(1,n);
P3=zeros(1,n);
b=zeros(n,1);
y1=zeros(N+1,1); y1(1)=1;
y2=zeros(N+1,1); y2(1)=1;
y3=zeros(N+1,1); y3(1)=1;
u=zeros(1,n);c1=zeros(1,n);c2=zeros(1,n);c3=zeros(1,n);
b1=zeros(1,n);b2=zeros(1,n);b3=zeros(1,n);b4=zeros(1,n);
b5=zeros(1,n);b6=zeros(1,n);b7=zeros(1,n);b8=zeros(1,n);
b9=zeros(1,n);b10=zeros(1,n);b11=zeros(1,n);b12=zeros(1,n);
E1=zeros(1,n);E2=zeros(1,n);E11=zeros(1,n);E12=zeros(1,n);Em=zeros(1,n);
F1=zeros(1,n);F2=zeros(1,n);E21=zeros(1,n);E22=zeros(1,n);E=zeros(1,n);pm=zeros(1,n);
for i=1:n
    t=t0(i);
w0=1;t00=10*T;w=1/T;
%     v=v0*cos(t/T);
v=v0*sin(pi*t/T).^2;
% v=v0*sech((t-t00)/T);
    for i1=1:N+1 
        for i2=1:N+1 
            if i1==i2
                m=-S+i2-1; 
                H0(i1,i2)=w0*m+lmada*(m.^2)/N;
            elseif i1+1==i2
                m=-S+i2-1;
                H0(i1,i2)=0.5*v*sqrt(S*(S+1)-m*(m-1));
            elseif i1-1==i2
                m=-S+i2-1;
                H0(i1,i2)=0.5*v*sqrt(S*(S+1)-m*(m+1));
            else
                H0(i1,i2)=0;
            end
        end
    end
    y1=grkt(t,y1,h,H0);
    Et=0;
    Et1=0;
    Et2=0;
    for i2=1:N+1
        m=-S+i2-1;
        Et1=Et1+abs(y1(i2))^2*(w0*m);           %哈密顿量期望
        Et2=Et2+abs(y1(i2))^2*((w0*m).^2);    %哈密顿量平方的期望
        Et=sqrt(Et2-Et1.^2);    
    end
    E111(i)=(Et);      %能量的涨落
    E222(i)=(Et1+N/2);      %能量没除N
%      E222(i)=(Et1+N/2)/N;      %能量除N
    P1(i)=E222(i)/(i*h);         %平均功率
    EE1(i)=N/2*(1+cos(pi-v0*t/2+v0*T/(4*pi)*sin(2*pi*t/T)));       %能量的解析式（2-36）
    EE2(i)=N/2*(1+cos(pi-v1*t/2+v1*T/(4*pi)*sin(2*pi*t/T)));       %能量的解析式（2-36）
    EE3(i)=N/2*(1+cos(pi-v2*t/2+v2*T/(4*pi)*sin(2*pi*t/T)));       %能量的解析式（2-36）
    EE4(i)=N/2*(1+cos(pi-v3*t/2+v3*T/(4*pi)*sin(2*pi*t/T)));       %能量的解析式（2-36）
    
    Ep1(i)=N*(1+cos(pi-v0*t/2+v0*T/(4*pi)*sin(2*pi*t/T)))/(2*t);       %功率的解析式（2-36）
    Ep2(i)=N*(1+cos(pi-v1*t/2+v1*T/(4*pi)*sin(2*pi*t/T)))/(2*t);       %功率的解析式（2-36）
    Ep3(i)=N*(1+cos(pi-v2*t/2+v2*T/(4*pi)*sin(2*pi*t/T)))/(2*t);       %功率的解析式（2-36）
    Ep4(i)=N*(1+cos(pi-v3*t/2+v3*T/(4*pi)*sin(2*pi*t/T)))/(2*t);       %功率的解析式（2-36）
    
    Eee1(i)=log2(sqrt((10-1)*exp(pi/2)))*abs(sin(pi-v0*t/2+v0*T/(4*pi)*sin(2*pi*t/T)));   %熵的解析式
    Eee2(i)=log2(sqrt((10-1)*exp(pi/2)))*abs(sin(pi-v1*t/2+v1*T/(4*pi)*sin(2*pi*t/T)));   %熵的解析式
    Eee3(i)=log2(sqrt((10-1)*exp(pi/2)))*abs(sin(pi-v2*t/2+v2*T/(4*pi)*sin(2*pi*t/T)));   %熵的解析式
    Eee4(i)=log2(sqrt((10-1)*exp(pi/2)))*abs(sin(pi-v3*t/2+v3*T/(4*pi)*sin(2*pi*t/T)));   %熵的解析式
    
    Et3=0;   Et4=0;
    Et4=Et4+sqrt((N)/4-(E222(i)-N/2)^2/N);       %能量涨落解析解，除以了N

    for i2=1:N+1      
    Et3=Et3-abs(y1(i2))^2*log(abs(y1(i2))^2)/log(2);        %冯诺依曼熵
   
   
    end
    E333(i)=Et3;
    E444(i)=Et4;
% u(i)=-v0/w*sin(w*t);
u(i)=pi-v0*(t/2-(T*sin(2*pi*t/T))/(4*pi));
% u(i)=pi-2*v0*T*(atan(tanh((t-t00)/(2*T))))-pi/2*v0*T;
v=pi;
c1(i)=-w0*cos(u(i))*sin(v);
c2(i)=-w0*sin(u(i));
c3(i)=w0*cos(u(i))*cos(v);
b1(i)=0;b4(i)=0;b5(i)=0;b6(i)=0;b7(i)=0;b10(i)=0;
% b2(i)=-4*sin(v0*t/4)^2*besselj(0,-v0*T/(4*pi))/v0-2*besselj(1,-v0*T/(4*pi))*T*((8*pi+(v0*T-4*pi))*cos(v0*t/2+2*pi*t/T)-(v0*T+4*pi)*cos(-v0*t/2+2*pi*t/T))/(16*pi^2-v0^2*T^2);
% b3(i)=-2*besselj(1,-v0*T/(4*pi))*(T*sin(-v0*t/2+2*pi*t/T)/(4*pi-v0*T)-T*sin(v0*t/2+2*pi*t/T)/(4*pi+v0*T))+2*sin(v0*t/2)/v0*besselj(0,-v0*T/(4*pi));
% b8(i)=lmada/N*((1-cos(v0*t))/v0*besselj(0,-v0*T/(2*pi))+2*besselj(1,-v0*T/(2*pi))*T*(-2*pi+2*pi*cos(2*pi*t/T)*cos(v0*t)+T*v0*sin(v0*t)*sin(2*pi*t/T))/(-4*pi^2+v0^2*T^2));
b2(i)=-4*sin(v0*T/4)^2*besselj(0,-v0*T/(4*pi))/v0-2*besselj(1,-v0*T/(4*pi))*16*pi*T*(sin(v0*T/4))^2/(16*pi^2-v0^2*T^2);
b3(i)=-2*besselj(1,-v0*T/(4*pi))*(-8*pi*T*(sin(v0*T/2))^2)/(16*pi^2-v0^2*T^2)+2*sin(v0*T/2)*besselj(0,-v0*T/(4*pi))/v0;
% b8(i)=lmada/N*(1-cos(v0*T))*(besselj(0,-v0*T/(2*pi))/v0+4*pi*T*besselj(1,-v0*T/(2*pi))/(-4*pi^2+v0^2*T^2));
% % % % % b2(i)=2*w0*besselj(1,-v0/w)*((1-cos(w*t))/w);
% % b2(i)=2*w0*besselj(1,-v0*T/(4*pi))*(pi*t-v0/4*(t^2-T^2*sin(pi*t/T)^2/(4*pi)));
% % b3(i)=w0*besselj(0,-v0*T/(4*pi))*t;
% % b8(i)=lmada/N*2*besselj(1,-2*v0*T/(4*pi))*((1-cos(w*t))/w);
% % b8(i)=lmada/N*2*besselj(1,-2*v0/w)*(pi*t-v0/4*(t^2-T^2*sin(pi*t/T)^2/(4*pi)));
% b9(i)=b8(i);
% % b11(i)=8*lmada/N*(besselj(1,-v0/w)).^2*(t/2-(sin(2*w*t))/(4*w));
% b12(i)=lmada/N*sin(v0*T)*(T+besselj(0,-v0*T/(4*pi))/v0-4*pi*T*besselj(1,-v0*T/(2*pi))/(-4*pi^2+v0^2*T^2));
% E11(i)=c3(i)+(b2(i)*c1(i)-b1(i)*c2(i))*sin(sqrt(b1(i).^2+b2(i).^2+b3(i).^2))/sqrt(b1(i).^2+b2(i).^2+b3(i).^2);
% E12(i)=((b1(i).^2+b2(i).^2)*c3(i)-b3(i)*(b1(i)*c1(i)+b2(i)*c2(i)))*(cos(sqrt(b1(i).^2+b2(i).^2+b3(i).^2))-1)/(b1(i).^2+b2(i).^2+b3(i).^2);
% E1(i)=E11(i)+E12(i);
% F1(i)=(b1(i)*b2(i)*(b6(i)+b7(i))+(b2(i).^2+b3(i).^2)*(b8(i)+b9(i)))*c1(i)-(b1(i)*b2(i)*(b8(i)+b9(i))+(b1(i).^2+b3(i).^2)*(b6(i)+b7(i)))*c2(i)+(b2(i)*b3(i)*(b6(i)+b7(i))-b1(i)*b3(i)*(b8(i)+b9(i)))*c3(i);
% F2(i)=(b8(i)+b9(i))*(b3(i)*c2(i)-b2(i)*c3(i));
% E21(i)=F1(i)*(sin(sqrt(b1(i).^2+b2(i).^2+b3(i).^2))-sqrt(b1(i).^2+b2(i).^2+b3(i).^2))/(b1(i).^2+b2(i).^2+b3(i).^2).^(3/2);
% E22(i)=F2(i)*(cos(sqrt(b1(i).^2+b2(i).^2+b3(i).^2))-1)/(b1(i).^2+b2(i).^2+b3(i).^2);
% E2(i)=E21(i)+E22(i);
% Es(i)=log2(sqrt((N-1)*exp(v/2))*(1-b8(i)*b9(i)*(b2(i)*c3(i)-b3(i)*c2(i))))*abs(sin(u(i)));
% E(i)=(c3(i)*N/2)^2;
% Em(i)=sqrt((E2(i)-E(i))/N);
% E(i)=(u(i).^2+(b2(i)+v).^2)*(cos(u(i).^2+(b2(i)+v).^2+b3(i).^2)-1)/(u(i).^2+(b2(i)+v).^2+b3(i).^2);
% E(i)=N/2-N/2*E1(i);
% P2(i)=E(i)/(i*h);
% E(i)=2.7*sin(u(i));
% E(i)=2*exp(-2*b2(i)-2*pi)*((b2(i)+pi)+exp(b2(i)+pi)*(2.3*(b2(i)+pi)*cos(c3(i)+b3(i))-4.6*(c3(i)+b3(i))*sin(c3(i)+b3(i))));
end
E1n(j,k)=max(E222);
E2n(j,k)=-max(E)+2;
% E3n(j,k)=max(E111);
% E4n(j,k)=max(E333);
% E5n(j,k)=E111(n);
% E6n(j,k)=max(E111);
% E7n(j,k)=max(E333);
end
end
figure(1)
% plot(t0,E,'o')
% hold on
% plot(t0,    Ep1,'linewidth',2)
% hold on
% plot(t0,    Ep2,'linewidth',2)
% hold on
% plot(t0,    Ep3,'linewidth',2)
% hold on
% plot(t0,    Ep4,'linewidth',2)

% figure(2)
% plot(t0,    EE1,'linewidth',2)
% hold on
% plot(t0,    EE2,'linewidth',2)
% hold on
% plot(t0,    EE3,'linewidth',2)
% hold on
% plot(t0,    EE4,'linewidth',2)

plot(t0,    Eee1,'linewidth',2)
hold on
plot(t0,    Eee2,'linewidth',2)
hold on
plot(t0,    Eee3,'linewidth',2)
hold on
plot(t0,    Eee4,'linewidth',2)

% xlabel('t');ylabel('E(t)')
% figure(2)
% plot(lmada0,E2n,'-','linewidth',4)
% hold on
% plot(lmada0,E1n,'--','linewidth',4)
axis([0 0.3 0 2.8])

 
xlabel('\lambda');ylabel('Emax')
% axis([-20,20,0.99,1]) hold on plot(lmada0,E2n,'linewidth',2) hold on
% plot(lmada0,E5n,'linewidth',2) xlabel('t');ylabel('E(\tau)') figure(2)
% plot(T0,E3n,'linewidth',2) hold on plot(T0,E4n,'linewidth',2) hold on
% plot(T0,E5n,'linewidth',2) xlabel('T');ylabel('E(\tau)') figure(1)
% plot(N0,E1n,'linewidth',2) figure(2) plot(N0,E2n,'linewidth',2) figure(3)
% plot(N0,E3n,'linewidth',2) figure(4) plot(N0,E4n,'linewidth',2)
% xlabel('t');ylabel('E(t)') f1=fopen('E_2.txt','wt+'); fprintf(f1,'%f %f
% \n',[t0;E2]); fclose(f1); f1=fopen('sigma.txt','wt+'); fprintf(f1,'%f %f
% \n',[t0;E1]); fclose(f1); f1=fopen('von.txt','wt+'); fprintf(f1,'%f %f
% \n',[t0;E3]); fclose(f1); f1=fopen('P.txt','wt+'); fprintf(f1,'%f %f
% \n',[t0;P2]); fclose(f1); f1=fopen('P_A.txt','wt+'); fprintf(f1,'%f %f %f
% %f\n',[A0;PP]); fclose(f1); plot(lmada0,e2,'linewidth',2)
% plot(N0,e1,'linewidth',2) hold on t=t0/pi;
toc;


% clear
% clc
% c0=0:0.1:5;
% d=length(c0);
% w0=0.1;
% m=length(w0);
% pn=zeros(d,m);
% for kk=1:1:d
%     c=c0(kk);
% 
% for i=1:1:m
%     w=w0(i);
% y=zeros(2,1); y(1)=1;
% h=0.001;
% k=1;
% g=0.01;
% T=1;
% t0=-T:h:T;
% n=length(t0);
% p=zeros(n,1);
% q=zeros(n,1);
%   for j=1:n
%     t=t0(j);
%     %tau=1000;
%     alpha=0.1;
%     Gamma_t=alpha*t;
%     %Gamma_t=(4/tau)*t-2; 
%     %omega=0.1;
%     omega_t=0.2; 
%     H0 = [Gamma_t/2 + c/2*(abs(y(2))^2 - abs(y(1))^2), omega_t/2;
%              omega_t/2, -Gamma_t/2 - c/2*(abs(y(2))^2 - abs(y(1))^2)];
%     dtheat=2*g*k./(4*g.^2+(w-k*t).^2);
%     Hcd=(1/2).*dtheat.*[0, -1i;1i, 0];
%     %H1 = 1i*omega/2*omega_t/(omega_t^2 + Gamma_t^2)*[0,-1i;1i,0];
%     H = H0+Hcd;
%     y=grkt(t,y,h,H);
%     %y=y./sqrt(sum(abs(y).^2));
%     y=y./sqrt(sum(abs(y).^2));
%     p(j)=abs(y(2)).^2;
%     q(j)=abs(y(1)).^2;
%     %E(j)=p(j);
%     %P(j)=E(j)/t;   
%   end 
%  pn(kk,i)=p(n); 
%  end
% 
% end
% 
% figure(1) 
% plot(c0,pn,'-','linewidth',2);
% hold on;
% plot(t0,q,'-','linewidth',2);









