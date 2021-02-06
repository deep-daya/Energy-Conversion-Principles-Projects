function [overall_des_R,overall_des_w_tur, overall_des_sigma, overall_des_eff]=task_3( rou, v1, alpha, Cl, n, omega)
%setting R, sigma, and other variables
R1=linspace(1,22,51); 
sigma1=linspace(0,1,51);
Kh1 = 0.085*R1;
rh1 = 0.1*R1; 
%Using task_1 to find W_tur and Efficiency for all different R and sigma
for i = 1:length(R1)
    disp(i);
    for j = 1:length(sigma1)
        [w_tur(i,j),eta(i,j)]=task_1(rou, v1, alpha, Cl, n, Kh1(i),sigma1(j), rh1(i), R1(i),omega );
    end
end
%finding a value closest to 1500W
[a,b]= find(w_tur<1550 & w_tur>=1500);
for l = 1:length(a)
    des_2_w_tur(l) = w_tur(a(l),b(l));
    des_2_eff(l) = eta(a(l),b(l));
end
[p,q] = max(des_2_eff);
des_2_eff = p;
des_2_w_tur = des_2_w_tur(q);
%finding a value of efficiency closest to 1 
[c,d]= find(eta<1.00 &eta>0.98);
for l = 1:length(c)
    des_1_w_tur(l) = w_tur(c(l),d(l));
    des_1_eff(l) = eta(c(l),d(l));
end
[e,f] = max(des_1_eff);
des_1_eff = e;
des_1_w_tur = des_1_w_tur(f);
des_1_R= R1(c(f));
des_1_sigma = sigma1(d(f));
des_2_R = R1(a(q));
des_2_sigma = sigma1(b(q));
des_R = [des_1_R,des_2_R];
des_eff = [des_1_eff,des_2_eff];
des_w_tur = [des_1_w_tur,des_2_w_tur];
%choosing between two alternatives for min sigma
[overall_des_sigma,index] = min([des_1_sigma,des_2_sigma]);
overall_des_R= des_R(index);
overall_des_w_tur = des_w_tur(index);
overall_des_eff = des_eff(index);
%plotting change in setup angle vs R
des_ov_rh = overall_des_R*0.1;
lambda = omega*overall_des_R/v1;
des_ov_c = des_ov_rh:0.001:overall_des_R;
des_ov_xi = alpha*pi/180 - atan(2*overall_des_R./3./lambda./des_ov_c);
des_ov_xi1=des_ov_xi*180/pi;
plot(des_ov_c, des_ov_xi1)
xlabel('Radius (m)')
ylabel('Setup Angle (Degrees)')
title('Setup Angle vs Radius')

end

function [w_tur, eta]=task_1(rou, v1, alpha, Cl, n, Kh,sigma, rh, R,omega )

lambda = omega*R/v1;
rh_bar= rh/R;
a=2*R/3/lambda;
b=1-rh_bar;

syms r;
w = Cl/3*(n*rou*lambda*Kh*(v1^2)/R) * ((r^2+a^2)^(0.5) *((r^2+a^2)*(lambda*v1/3/R*(1+sigma*rh_bar/b) - sigma*r*lambda*v1/4/(R^2)/b) + log(r+(r^2+a^2)^(0.5))*(2*sigma*R^2*v1)/(81*lambda^3*b)));
w_b = 16/27*0.5*rou*v1^3*pi *r^2;

%turbine work
ans = subs(w,r, [rh R]);
w_tur =double(ans(2)-ans(1));

%betz work 
ans1 = subs(w_b,r, [rh R]);
w_betz = double(ans1(2)-ans1(1));

eta=w_tur/w_betz;

%c=rh:0.001:R;
%xi = alpha*pi/180 - atan(2*R./3./lambda./c);
%xi1=xi*180/pi;
%plot(c, xi1)
 

end