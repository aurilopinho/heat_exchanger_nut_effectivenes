%clear all

M = 3;
N = 3;

NUT = linspace (0, 10, N);
Cr = linspace (0, 1, M);

[NUT, Cr] = meshgrid (NUT, Cr);

n = 2;
err_max = 1.0E-05;

eps = zeros (M, N);
NUT_ = zeros (M, N);

%   Parallel flow
% % ss
for I = 1: M
    for J = 1: N
        eps(I,J) = heat_exchanger_eps (NUT(I,J), Cr(I,J), 'parallel flow');
        NUT_(I,J) = heat_exchanger_nut (eps(I,J), Cr(I,J), 'parallel flow');
    end
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Parallel flow ss: %.2E\n', err_m);
end

% % sv
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,1), 'parallel flow');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,1), 'parallel flow');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Parallel flow sv: %.2E\n', err_m);
end

% % vs
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(1,J), Cr(:,J), 'parallel flow');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'parallel flow');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Parallel flow vs: %.2E\n', err_m);
end

% % rr
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,:), 'parallel flow');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,:), 'parallel flow');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Parallel flow vv: %.2E\n', err_m);
end

% % cc
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(:,J), Cr(:,J), 'parallel flow');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'parallel flow');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Parallel flow vv: %.2E\n', err_m);
end

% % mm
eps = heat_exchanger_eps (NUT, Cr, 'parallel flow');
NUT_ = heat_exchanger_nut (eps, Cr, 'parallel flow');

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Parallel flow mm: %.2E\n', err_m);
end

fprintf ('Parallel Flow: error max: %0.4E\n', err_m);

%   Counter flow
% % ss
for I = 1: M
    for J = 1: N
        eps(I,J) = heat_exchanger_eps (NUT(I,J), Cr(I,J), 'counter flow');
        NUT_(I,J) = heat_exchanger_nut (eps(I,J), Cr(I,J), 'counter flow');
    end
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Counter flow ss: %.2E\n', err_m);
end

% % sv
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,1), 'counter flow');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,1), 'counter flow');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Counter flow sv: %.2E\n', err_m);
end

% % vs
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(1,J), Cr(:,J), 'counter flow');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'counter flow');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Counter flow vs: %.2E\n', err_m);
end

% % rr
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,:), 'counter flow');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,:), 'counter flow');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Counter flow vv: %.2E\n', err_m);
end

% % cc
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(:,J), Cr(:,J), 'counter flow');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'counter flow');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Counter flow vv: %.2E\n', err_m);
end

% % mm
eps = heat_exchanger_eps (NUT, Cr, 'counter flow');
NUT_ = heat_exchanger_nut (eps, Cr, 'counter flow');

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Counter flow mm: %.2E\n', err_m);
end

fprintf ('Counter Flow: error max: %0.4E\n', err_m);

%   Shell and tube, one shell pass (2, 4, 6, ... tubes passes)
% % ss
for I = 1: M
    for J = 1: N
        eps(I,J) = heat_exchanger_eps (NUT(I,J), Cr(I,J), 'single shell pass');
        NUT_(I,J) = heat_exchanger_nut (eps(I,J), Cr(I,J), 'single shell pass');
    end
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, one shell pass ss: %.2E\n', err_m);
end

% % sv
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,1), 'single shell pass');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,1), 'single shell pass');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, one shell pass sv: %.2E\n', err_m);
end

% % vs
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(1,J), Cr(:,J), 'single shell pass');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'single shell pass');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, one shell pass vs: %.2E\n', err_m);
end

% % rr
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,:), 'single shell pass');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,:), 'single shell pass');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, one shell pass vv: %.2E\n', err_m);
end

% % cc
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(:,J), Cr(:,J), 'single shell pass');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'single shell pass');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, one shell pass vv: %.2E\n', err_m);
end

% % mm
eps = heat_exchanger_eps (NUT, Cr, 'single shell pass');
NUT_ = heat_exchanger_nut (eps, Cr, 'single shell pass');

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, one shell pass mm: %.2E\n', err_m);
end

fprintf ('Shell and tube, 1 shell pass: error max: %0.4E\n', err_m);

%   Shell and tube, 2 shell passes (2n, 4n, 6n, ... tubes passes)
% % ss
for I = 1: M
    for J = 1: N
        eps(I,J) = heat_exchanger_eps (NUT(I,J), Cr(I,J), 'multiple shell passes', n);
        NUT_(I,J) = heat_exchanger_nut (eps(I,J), Cr(I,J), 'multiple shell passes', n);
    end
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, 2 shell passes ss: %.2E\n', err_m);
end

% % sv
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,1), 'multiple shell passes', n);
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,1), 'multiple shell passes', n);
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, 2 shell passes sv: %.2E\n', err_m);
end

% % vs
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(1,J), Cr(:,J), 'multiple shell passes', n);
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'multiple shell passes', n);
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, 2 shell passes vs: %.2E\n', err_m);
end

% % rr
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,:), n);
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,:), n);
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, 2 shell passes vv: %.2E\n', err_m);
end

% % cc
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(:,J), Cr(:,J), 'multiple shell passes', n);
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'multiple shell passes', n);
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, 2 shell passes vv: %.2E\n', err_m);
end

% % mm
eps = heat_exchanger_eps (NUT, Cr, 'multiple shell passes', n);
NUT_ = heat_exchanger_nut (eps, Cr, 'multiple shell passes', n);

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Shell and tube, 2 shell passes mm: %.2E\n', err_m);
end

fprintf ('Shell and tube, 2 shell passes: error max: %0.4E\n', err_m);

%   Cross Flow with Single Pass and both fluid unmixed
% % ss
for I = 1: M
    for J = 1: N
        eps(I,J) = heat_exchanger_eps (NUT(I,J), Cr(I,J), 'cross flow both unmixed');
        NUT_(I,J) = heat_exchanger_nut (eps(I,J), Cr(I,J), 'cross flow both unmixed');
    end
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and both fluid unmixed ss: %.2E\n', err_m);
end

% % sv
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,1), 'cross flow both unmixed');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,1), 'cross flow both unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and both fluid unmixed sv: %.2E\n', err_m);
end

% % vs
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(1,J), Cr(:,J), 'cross flow both unmixed');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'cross flow both unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and both fluid unmixed vs: %.2E\n', err_m);
end

% % rr
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,:), 'cross flow both unmixed');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,:), 'cross flow both unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and both fluid unmixed rr: %.2E\n', err_m);
end

% % cc
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(:,J), Cr(:,J), 'cross flow both unmixed');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'cross flow both unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and both fluid unmixed cc: %.2E\n', err_m);
end

% % mm
eps = heat_exchanger_eps (NUT, Cr, 'cross flow both unmixed');
NUT_ = heat_exchanger_nut (eps, Cr, 'cross flow both unmixed');

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and both fluid unmixed mm: %.2E\n', err_m);
end

fprintf ('Cross flow with and both fluid unmixed: error max: %0.4E\n', err_m);

%   Cross flow with single pass and C_max fluid unmixed
% % ss
for I = 1: M
    for J = 1: N
        eps(I,J) = heat_exchanger_eps (NUT(I,J), Cr(I,J), 'cross flow Cmax unmixed');
        NUT_(I,J) = heat_exchanger_nut (eps(I,J), Cr(I,J), 'cross flow Cmax unmixed');
    end
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_max fluid unmixed ss: %.2E\n', err_m);
end

% % sv
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,1), 'cross flow Cmax unmixed');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,1), 'cross flow Cmax unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_max fluid unmixed sv: %.2E\n', err_m);
end

% % vs
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(1,J), Cr(:,J), 'cross flow Cmax unmixed');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'cross flow Cmax unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_max fluid unmixed vs: %.2E\n', err_m);
end

% % rr
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,:), 'cross flow Cmax unmixed');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,:), 'cross flow Cmax unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_max fluid unmixed rr: %.2E\n', err_m);
end

% % cc
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(:,J), Cr(:,J), 'cross flow Cmax unmixed');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'cross flow Cmax unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_max fluid unmixed cc: %.2E\n', err_m);
end

% % mm
eps = heat_exchanger_eps (NUT, Cr, 'cross flow Cmax unmixed');
NUT_ = heat_exchanger_nut (eps, Cr, 'cross flow Cmax unmixed');

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_max fluid unmixed mm: %.2E\n', err_m);
end

fprintf ('Cross flow with and C_max fluid unmixed: error max: %0.4E\n', err_m);

%   Cross flow with single pass and C_min fluid unmixed
% % ss
for I = 1: M
    for J = 1: N
        eps(I,J) = heat_exchanger_eps (NUT(I,J), Cr(I,J), 'cross flow Cmin unmixed');
        NUT_(I,J) = heat_exchanger_nut (eps(I,J), Cr(I,J), 'cross flow Cmin unmixed');
    end
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_min fluid unmixed ss: %.2E\n', err_m);
end

% % sv
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,1), 'cross flow Cmin unmixed');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,1), 'cross flow Cmin unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_min fluid unmixed sv: %.2E\n', err_m);
end

% % vs
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(1,J), Cr(:,J), 'cross flow Cmin unmixed');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'cross flow Cmin unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_min fluid unmixed vs: %.2E\n', err_m);
end

% % rr
for I = 1: M
    eps(I,:) = heat_exchanger_eps (NUT(I,:), Cr(I,:), 'cross flow Cmin unmixed');
    NUT_(I,:) = heat_exchanger_nut (eps(I,:), Cr(I,:), 'cross flow Cmin unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_min fluid unmixed rr: %.2E\n', err_m);
end

% % cc
for J = 1: N
    eps(:,J) = heat_exchanger_eps (NUT(:,J), Cr(:,J), 'cross flow Cmin unmixed');
    NUT_(:,J) = heat_exchanger_nut (eps(:,J), Cr(:,J), 'cross flow Cmin unmixed');
end

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_min fluid unmixed cc: %.2E\n', err_m);
end

% % mm
eps = heat_exchanger_eps (NUT, Cr, 'cross flow Cmin unmixed');
NUT_ = heat_exchanger_nut (eps, Cr, 'cross flow Cmin unmixed');

err = 1 - NUT ./ NUT_;
err_m = max (max (abs (err)));
if (err_m > err_max)
    error ('Cross flow with single pass and C_min fluid unmixed mm: %.2E\n', err_m);
end

fprintf ('Cross flow with and C_min fluid unmixed: error max: %0.4E\n', err_m);
