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

c=rh:0.001:R;
xi = alpha*pi/180 - atan(2*R./3./lambda./c);
xi1=xi*180/pi;
plot(c, xi1)
title ('Setup Angle vs. Radius');
xlabel('Radius (m)', 'FontWeight','bold');
ylabel('Setup Angle (degrees)', 'FontWeight','bold');
 

end
