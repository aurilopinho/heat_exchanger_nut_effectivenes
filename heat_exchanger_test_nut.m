%clear all

% NUT, Cr
N = 26;
eps = linspace (0, 1, N);

M = 5;
Cr = linspace (0, 1, M)';

% legend to plot
s_leg = num2str (Cr, 'Cr = %0.2f');

[Cr, eps] = meshgrid (Cr, eps);

% Parallel flow
NUT = heat_exchanger_nut (eps, Cr, 'parallel flow');

plot (eps, NUT, '-*');
title ('Parallel Flow Heat Exchange');
xlabel ('eps');
ylabel ('NUT');
legend (s_leg, 'Location','northwest');
grid

pause

% Counter flow
NUT = heat_exchanger_nut (eps, Cr, 'counter flow');

plot (eps, NUT, '-*');
title ('Counter Flow Heat Exchange');
xlabel ('eps');
ylabel ('NUT');
legend (s_leg, 'Location','northwest');
grid

pause

% Shell and tube, one shell pass (2, 4, 6, ... tubes passes)
NUT = heat_exchanger_nut (eps, Cr, 'single shell pass');

plot (eps, NUT, '-*');
title ('Shell and Tube Heat Exchange: one shell pass');
xlabel ('eps');
ylabel ('NUT');
legend (s_leg, 'Location','northwest');
grid

pause

% Shell and tube, 2 shell passes (2n, 4n, 6n, ... tubes passes)
NUT = heat_exchanger_nut (eps, Cr, 'multiple shell passes', 2);

plot (eps, NUT, '-*');
title ('Shell and Tube Heat Exchange: multiple shell passes');
xlabel ('eps');
ylabel ('NUT');
legend (s_leg, 'Location','northwest');
grid

pause

% Cross flow with Single Pass and Both Fluid Unmixed
NUT = heat_exchanger_nut (eps, Cr, 'cross flow both unmixed');

plot (eps, NUT, '-*');
title ('Cross Flow Heat Exchange: both unmixed');
xlabel ('eps');
ylabel ('NUT');
legend (s_leg, 'Location','northwest');
grid

pause

%   Cross flow with Single Pass and C_min Fluid Unmixed
NUT = heat_exchanger_nut (eps, Cr, 'cross flow Cmin unmixed');

plot (eps, NUT, '-*');
title ('Cross Flow Heat Exchange: C_{min} unmixed');
xlabel ('eps');
ylabel ('NUT');
legend (s_leg, 'Location','northwest');
grid

pause

%   Cross Flow with Single Pass and C_min Fluid Unmixed
NUT = heat_exchanger_nut (eps, Cr, 'cross flow Cmax unmixed');

plot (eps, NUT, '-*');
title ('Cross Flow Heat Exchange: C_{max} unmixed');
xlabel ('eps');
ylabel ('NUT');
legend (s_leg, 'Location','northwest');
grid

pause

close all
