% Setting up variables and values
t= 13;
lambda=37.9;
d= 120;
zeta= 200;
epsilon = 36;
% Using equations 2, 3, 6, 7, 8, 9 and solving them to get different
% parameters
%Finding Hour angle
alpha = 15*(t-12);
%Finding Declination Angle
del = 23.44*sind((360/365.25)*(d-80));
%Finding Zenith angle
cos_zen_ang= sind(lambda)*sind(del)+ cosd(lambda)*cosd(del)*cosd(alpha);
chi = acosd(cos_zen_ang);
%finding Solar Azimuth angle
tan_xi = sind(alpha) /((sind(lambda)*cosd(alpha)-cosd(lambda)*tand(del)));
if alpha >= 0 && tan_xi >= 0
    xi = 180 + atand(tan_xi);
elseif alpha >= 0 &&  tan_xi < 0
    xi = 360 + atand(tan_xi);
elseif alpha < 0 && tan_xi >= 0
    xi = atand(tan_xi);
elseif alpha<0 && tan_xi < 0
    xi = 180+ atand(tan_xi);
end 
A = 1310;
B = 0.18;
%Finding Direct Normal Radiation
I_dn = A*exp(-B/(sind(90-chi)));
%Finding Direct Solar Flux per unit area
Id = I_dn*(cos_zen_ang*cosd(epsilon)+sind(epsilon)*sind(chi)*cosd(xi - zeta));
Id