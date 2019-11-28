%Setting up parameters
t_0= 10;
t_i=16;
lambda=37.9;
d= 120;
zeta= 200;
epsilon = 36;
%Setting up time as a per second change
t = t_0:(1/3600):t_i;
%Setting up Parameters from task 1
for i = 1:length(t) 
    alpha(i) = 15*(t(i)-12);
end
del = 23.44*sind((360/365.25)*(d-80));
cos_zen_ang = zeros(0,length(alpha));
chi = zeros(0,length(alpha));
tan_xi = zeros(0,length(alpha));
for i = 1:length(t)
cos_zen_ang(i)= sind(lambda)*sind(del)+ cosd(lambda)*cosd(del)*cosd(alpha(i));
chi(i) = acosd(cos_zen_ang(i));
end
for i = 1:length(t)
    tan_xi(i) = (sind(alpha(i))) /((sind(lambda)*cosd(alpha(i))-cosd(lambda)*tand(del)));
end
for i = 1:length(t)
    if alpha(i) > 0 && tan_xi(i) > 0
        xi(i) = 180 + atand(tan_xi(i));
    elseif alpha(i) > 0 &&  tan_xi(i) < 0
        xi(i) = 360 + atand(tan_xi(i));
    elseif alpha(i) < 0 && tan_xi(i) >0
        xi(i) = atand(tan_xi(i));
    elseif alpha(i)<0 && tan_xi(i) < 0
        xi(i) = 180+ atand(tan_xi(i));
    end
end 
A = 1310;
B = 0.18;
for i = 1:length(chi)
    I_dn(i) = A*exp(-B/(sind(90-chi(i))));
end
%Finding Direct solar flux per unit area per second
for i = 1:length(t)
     Id_i(i) = I_dn(i)*(cos_zen_ang(i)*cosd(epsilon)+sind(epsilon).*sind(chi(i))*cosd(xi(i) - zeta));
end
%Using the trapezoidal rule for finding the total Incident radiation
Id_f = 0;
for i = 1:length(Id_i)-1
Id_f = ((((t_i-t_0)*3600)/length(t))*(((Id_i(i)) + Id_i(i+1))/2))+ Id_f;
end
Id_f