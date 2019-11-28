function [n_colll,Id] = task_3(Ac,cp,m,t_0,t_i,lambda,d, zeta,epsilon, alpha_c, tau_g, h_convo,del_g,kg,h_convi, del_ins, k_ins, T_i,T_a)
%Setting up time and given variables
t = t_0:1/3600:t_i;
for i = 1:length(t) 
    alpha(i) = 15*(t(i)-12);
end
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
%Finding the conductance, heat removal factors and collector efficiencies
%per second
U_loss = ((1/(h_convo*Ac) + del_g/(kg*Ac) + 1/(h_convi*Ac))^-1 + (1/(h_convo*Ac) + (del_ins/(k_ins*Ac)))^-1)/Ac;
F_r = (1- exp(-Ac*U_loss/(m*cp)))/(Ac*U_loss/(m*cp));
Id = 0;
for i = 1:length(t)-1
Id = (((t_i-t_0)*3600/length(t))*(((Id_i(i)) + Id_i(i+1))/2))+Id;
end
n_coll = F_r*(tau_g *alpha_c - (U_loss./Id_i)*(T_i-T_a));
%Finding total energy absorbed in watts
Id_ncol = 0;
for i = 1:length(t)-1
Id_ncol = (((t_i-t_0)*3600/length(t))*(((Id_i(i))*n_coll(i) + n_coll(i)*Id_i(i+1))/2))+Id_ncol;
end
n_colll = Id_ncol/Id;