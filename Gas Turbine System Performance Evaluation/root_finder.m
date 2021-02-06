function [T3_f, cp23] = root_finder(Qs, T2_r, n_dot_air)
syms cp_23 T3
    eqn = [cp_23 == 28.11+ 0.1967*10^(-2)*(0.5*(T3+T2_r)) + 0.4802*10^(-5)*(0.5*(T3+T2_r))^2 - 1.966*10^(-9)*(0.5*(T3+T2_r))^3,T3 == Qs/(n_dot_air*cp_23)+T2_r];
v = solve(eqn, [cp_23, T3]);
cp_23= vpa(v.cp_23);
T3 = vpa(v.T3);
p = find(T3 >800& T3 < 1500);
T3_f = double(T3(p));
cp23 = double(cp_23(p));
end
