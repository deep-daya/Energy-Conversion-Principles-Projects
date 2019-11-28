function [T_out, Id_ncol] = task_4(t_0,t_i,d, zeta,epsilon, T_i,m)
%Setting up parameters as in previous tasks
lambda = 37.9;
t = t_0:1/3600:t_i;
for i = 1:length(t)
    T_a(i)= (7/6)*t(i)-(8/3);
end
Ac = 3.25;
alpha_c =0.85;
tau_g = 0.89;
h_convo = 7;
del_g = 0.007;
kg = 1.3;
h_convi = 3.1;
del_ins = .06;
k_ins = 0.045;
cp = 4186;
T_in = [];
for i = 1:length(t) 
    alpha(i) = 15*(t(i)-12);
%Converting Tin to the same length as T_ambient for task 4
end
if length(T_i)==1
    for i =1:length(t)
        T_in(end+1) = T_i;
    end
else
    T_in = T_i;
end;
del = 23.44*sind((360/365.25)*(d-80));
cos_zen_ang = zeros(0,length(alpha));
for i = 1:length(t)
cos_zen_ang(i)= sind(lambda)*sind(del)+ cosd(lambda)*cosd(del)*cosd(alpha(i));
chi(i) = acosd(cos_zen_ang(i));
end
for i = 1:length(t)
    tan_xi(i) = sind(alpha(i)) /((sind(lambda)*cosd(alpha(i))-cosd(lambda)*tand(del)));
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
    Id_i(i) = I_dn(i)*(cos_zen_ang(i)*cosd(epsilon)+sind(epsilon)*sind(chi(i))*cosd(xi(i) - zeta));
end
U_loss = ((1/(h_convo*Ac) + del_g/(kg*Ac) + 1/(h_convi*Ac))^-1 + (1/(h_convo*Ac) + (del_ins/(k_ins*Ac)))^-1)/Ac;
F_r = (1- exp(-Ac*U_loss/(m*cp)))/(Ac*U_loss/(m*cp));
Id = 0;
for i = 1:length(t)-1
Id = (((t_i-t_0)*3600/length(t))*(((Id_i(i)) + Id_i(i+1))/2))+Id;
end
for i =1: length(t)
n_colli(i) = F_r*(tau_g *alpha_c - (U_loss./Id_i(i)).*(T_in(i)-T_a(i)));
end
Id_ncol = 0;
for i = 1:length(t)-1
Id_ncol = (((t_i-t_0)*3600/length(t))*(((Id_i(i))*n_colli(i) + n_colli(i+1)*Id_i(i+1))/2))+Id_ncol;
end
n_colll = Id_ncol/Id;
%Finding T_out per second using equation 13
for i = 1:length(t)
T_out(i)= (Id_i(i)*n_colli(i)*Ac)/(cp*m)+ T_in(i);
end
%Finding T_avg
T_out_avg = sum(T_out(:))/length(T_out);
