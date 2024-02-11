clear all;
close all;
clc;

%%
%naca 4 digit aerofoil (naca uvw) nomenclature
u=input('first digit of naca = ');
v=input('second digit of naca = ');
w=input('last two digit of naca = ');
% u=2;v=4;w=12;
ymc_c=u*0.01;
xmc_c=v*0.1;
tm=w*0.01;
U=100;
%chord
c=1;

%%
%cosine clustring
%choose even number of points on circle
n=160;
theta0=2*pi/(n-1);
i=1:n/2;
x_c=0.5*c*(1-cos((i-0.5)*theta0));

%%
%camber line formation
for j=1:length(x_c)
    if(x_c(j)>=0) && (x_c(j)<=xmc_c)
        ycx_c(j)=ymc_c*(2*(x_c(j)/xmc_c)-(x_c(j)/xmc_c)^2);
        dyc_dx(j)=ymc_c*(2*(1/xmc_c)-(2*(x_c(j))/(xmc_c)^2));
    else
        ycx_c(j)=ymc_c*(2*((1-x_c(j))/(1-xmc_c))-((1-x_c(j))/(1-xmc_c))^2);
        dyc_dx(j)=ymc_c*(2*(-1/(1-xmc_c))+(2*(1-x_c(j))/(1-xmc_c)^2));
    end
end

%%
%
theta0=0:theta0:pi;
syms alpha
a0=alpha-(1/pi)*trapz(theta0,dyc_dx,2);
e=trapz(theta0,1+cos(theta0));

a1=(2/pi)*trapz(theta0,dyc_dx.*cos(theta0));
b=pi/2;
h=a1*b;
syms alpha
eqn=a0*e+h==0;
f=solve(eqn,alpha);
a_l0=double(f);
a=(-10:1:15)*pi/180;
RA=1:1:8;

N=100;
i=1:1:100;
theta=(i-1)*pi/(N-1);

%elliptic wing
%cl and cd
for o=1:length(RA)
    cl_e=2*pi.*(a-a_l0)./(1+2./RA(o));
    cd_e=cl_e.^2/(pi*RA(o));
%     cl_e_tilda=cl_e;
    figure(1)
    plot(a*180/pi,cl_e);
    hold on;
    xlabel('--alpha (in degree) ---->');
    ylabel('--- lift coefficient elliptic---->');
    legend('RA=1','RA=2','RA=3','RA=4','RA=5','RA=6','RA=7','RA=8');
    figure(2)
    plot(cl_e,cd_e);
    hold on;
    xlabel('-- lift coefficient elliptic ---->');
    ylabel('--- induced drag coefficient elliptic---->');
    legend('RA=1','RA=2','RA=3','RA=4','RA=5','RA=6','RA=7','RA=8');
end

for o=1:length(RA)
%     fix angle of attack 10 degree
    cl_e_tilda=2*pi.*(10*pi/180-a_l0)./(1+2./RA(o));
%     fix angle of attack 10 degree
    a_i_e=cl_e_tilda/(pi*RA(o));
    figure(3)
    plot([theta(1) theta(end)],[cl_e_tilda cl_e_tilda],'LineWidth',1);
    hold on;
    ylabel('-- sectional lift coefficient elliptic ---->');
    xlabel('---  spanwise location ---->');
    legend('RA=1','RA=2','RA=3','RA=4','RA=5','RA=6','RA=7','RA=8');
    figure(4)
    plot([theta(1) theta(end)],[a_i_e a_i_e],'LineWidth',1);
    hold on;
    ylabel('-- downwash elliptic ---->');
    xlabel('---  spanwise location ---->');
    legend('RA=1','RA=2','RA=3','RA=4','RA=5','RA=6','RA=7','RA=8');
end



%rectangular wing
%cl
C=zeros(N,N);

for o=1:length(RA)
    for j=1:N
        C(1,j)=j^2;
    end
    for j=1:N
        for i=2:(N-1)
            C(i,j)=(2*RA(o)./pi+j./sin(theta(i))).*sin(j*theta(i));
        end
    end
    for j=1:N
        C(N,j)=(j^2)*(-1)^(j+1);
    end
    B=ones(N,1);
    A=C\B;
    cl_r=pi*RA(o)*A(1)*(a-a_l0);
    sigma=3*(A(3)/A(1))^2+5*(A(5)/A(1))^2+7*(A(7)/A(1))^2;
    cd_r=cl_r.^2.*(1+sigma)/(pi*RA(o));
%     fix angle of attack 10 degree
    cl_r_tilda=4*RA(o)*(10*pi/180-a_l0)*(A(1)*sin(theta)+A(3)*sin(3*theta)+A(5)*sin(5*theta)+A(7)*sin(7*theta));
%     fix angle of attack 10 degree
    a_i_r=(10*pi/180-a_l0)*(1*A(1)*sin(theta)+3*A(3)*sin(3*theta)+5*A(5)*sin(5*theta)+7*A(7)*sin(7*theta))./sin(theta);
    figure(5)
    plot(a*180/pi,cl_r);
    hold on;
    xlabel('--alpha (in degree) ---->');
    ylabel('--- lift coefficient rectangular ---->');
    legend('RA=1','RA=2','RA=3','RA=4','RA=5','RA=6','RA=7','RA=8');
    figure(6)
    plot(cl_r,cd_r);
    hold on;
    xlabel('-- lift coefficient rectangular ---->');
    ylabel('--- induced drag coefficient rectangular ---->');
    legend('RA=1','RA=2','RA=3','RA=4','RA=5','RA=6','RA=7','RA=8');
    figure(7)
    plot(-0.5*cos(theta),cl_r_tilda);
    hold on;
    ylabel('-- sectional lift coefficient rectangular ---->');
    xlabel('---  spanwise location ---->');
    legend('RA=1','RA=2','RA=3','RA=4','RA=5','RA=6','RA=7','RA=8');
    figure(8)
    plot(-0.5*cos(theta),a_i_r);
    hold on;
    ylabel('-- downwash rectangular ---->');
    xlabel('---  spanwise location ---->');
    legend('RA=1','RA=2','RA=3','RA=4','RA=5','RA=6','RA=7','RA=8');
end



