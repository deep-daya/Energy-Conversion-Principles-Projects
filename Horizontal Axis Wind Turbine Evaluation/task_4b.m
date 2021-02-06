function [fract]=task_4b(I_L, I_cell, V_rev)
V_L = V_rev - 0.11 - 0.041*I_L;
V_cell = V_rev + 0.11 + 0.041*I_cell;
fract = V_L/V_cell;
end