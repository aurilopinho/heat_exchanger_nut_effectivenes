%clear all

% NUT, Cr
N = 51;
NUT = linspace (0, 10, N);

M = 5;
Cr = linspace (0, 1, M)';

% legend to plot
s_leg = num2str (Cr, 'Cr = %0.2f');

[Cr, NUT] = meshgrid (Cr, NUT);

% Parallel flow
eps = heat_exchanger_eps (NUT, Cr, 'parallel flow');

plot (NUT, eps, '-*');
title ('Parallel Flow Heat Exchange');
xlabel ('NUT');
ylabel ('eps');
legend (s_leg, 'Location','southeast');
grid

pause

% Counter flow
eps = heat_exchanger_eps (NUT, Cr, 'counter flow');

plot (NUT, eps, '-*');
title ('Counter Flow Heat Exchange');
xlabel ('NUT');
ylabel ('eps');
legend (s_leg, 'Location','southeast');
grid

pause

% Shell and tube, one shell pass (2, 4, 6, ... tubes passes)
eps = heat_exchanger_eps (NUT, Cr, 'single shell pass');

plot (NUT, eps, '-*');
title ('Shell and Tube Heat Exchange: one shell pass');
xlabel ('NUT');
ylabel ('eps');
legend (s_leg, 'Location','southeast');
grid

pause

% Shell and tube, 2 shell passes (2n, 4n, 6n, ... tubes passes)
eps = heat_exchanger_eps (NUT, Cr, 'multiple shell passes', 3);

plot (NUT, eps, '-*');
title ('Shell and Tube Heat Exchange: multiple shell passes');
xlabel ('NUT');
ylabel ('eps');
legend (s_leg, 'Location','southeast');
grid

pause

% Cross flow with Single Pass and Both Fluid Unmixed
eps = heat_exchanger_eps (NUT, Cr, 'cross flow both unmixed');

plot (NUT, eps, '-*');
title ('Cross Flow Heat Exchange: both unmixed');
xlabel ('NUT');
ylabel ('eps');
legend (s_leg, 'Location','southeast');
grid

pause

%   Cross flow with Single Pass and C_min Fluid Unmixed
eps = heat_exchanger_eps (NUT, Cr, 'cross flow Cmin unmixed');

plot (NUT, eps, '-*');
title ('Cross Flow Heat Exchange: C_{min} unmixed');
xlabel ('NUT');
ylabel ('eps');
legend (s_leg, 'Location','southeast');
grid

pause

%   Cross Flow with Single Pass and C_min Fluid Unmixed
eps = heat_exchanger_eps (NUT, Cr, 'cross flow Cmax unmixed');

plot (NUT, eps, '-*');
title ('Cross Flow Heat Exchange: C_{max} unmixed');
xlabel ('NUT');
ylabel ('eps');
legend (s_leg, 'Location','southeast');
grid

pause

close all
