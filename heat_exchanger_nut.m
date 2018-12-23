function [NUT] = heat_exchanger_nut (eps, Cr, type, n)

    % This program calculate heat exchange NUT from
    %  effectiveness and heat capacity ratio folowing textbook
    %  Fundamentals of Heat and Mass Transfer (7th ed.) 
    %  by BERGMAN, T. L., LAVINE, A. S., INCROPERA, F. P. 
    %  and DEWITT, D. P.. Table 11.4 pp 724.
    %
    % Input parameters explained:
    %
    %   eps: heat exchanger effectiveness (q / qmax)
    %        q: actual heat transfer rate for a heat exchanger
    %        qmax: maximum possible heat transfer rate
    %
    %   obs: eps must be >= 0 and <= 1 (error) or < effectiveness limit from 
    %        type and Cr (warning).
    %
    %   Cr: heat capacity ratio (Cmax / Cmin)
    %       Cmax: greater heat capacity between both fluids
    %       Cmin: smaller heat capacity between both fluids
    %
    %   obs: Cr must be >= 0 and <= 1
    %
    %   obs: Both, eps and Cr can be scallar, vector or matrix.
    %
    %   type: heat exchange type. Valid entries
    %
    %         1 - parallel flow
    %         2 - counter flow
    %         3 - single shell pass
    %         4 - multiple shell passes
    %         5 - cross flow both unmixed
    %         6 - cross flow Cmax unmixed
    %         7 - cross flow Cmin unmixed
    %
    %   obs: type can use number or case sensitive string.
    %
    %   n: number of shell passes on multiple shell passes heat exchange type.
    %
    % Output parameter explained:
    %
    %   eps: number of transfer units (U A / Cmin)
    %        U: overall heat transfer coefficient
    %        A: heat transfer area
    %        Cmin: smaller heat capacity between both fluids
    %
    % Examples of the program in use:
    %
    %   Input:  eps = 0.5;
    %   Input:  Cr = 0.5;
    %   Input:  NUT = heat_exchanger_nut (eps, Cr, 'parallel flow')
    %   Output: NUT =
    %
    %               0.9242
    %
    %   Input:  eps = [0 0.5 0.99];
    %   Input:  Cr = 0.5;
    %   Input:  NUT = heat_exchanger_nut (eps, Cr, 2)
    %   Output: NUT =
    %
    %               0    0.8109    7.8439
    %
    %   Input:  eps = 0.5;
    %   Input:  Cr = [0 0.5 1];
    %   Input:  NUT = heat_exchanger_nut (eps, Cr, 'single shell pass')
    %   Output: NUT =
    %
    %               0.6931    0.8608    1.2465
    %
    %   Input:  eps = [0 0.5 0.99];
    %   Input:  Cr = [0 0.5 1];
    %   Input:  NUT = heat_exchanger_nut (eps, Cr, 4, 2)
    %   Output  Warning: There is any eps >= eps_lim to Cr in multiple shell 
    %           passes heat exchanger. The result is NaN. 
    %           > In heat_exchanger_nut (line 185)
    %
    %   Output: NUT =
    %
    %               0    0.4112       NaN
    %
    %   Input:  eps = [0 0.5 0.99];
    %   Input:  Cr = [0 0.5 1];
    %   Input:  [eps Cr] = meshgrid (eps, Cr);
    %   Input:  NUT = heat_exchanger_nut (eps, Cr, 'cross flow both unmixed')
    %   Output: NUT =
    %
    %               0    0.6931    4.6052
    %               0    0.8584   42.6732
    %               0    1.1498  512.0000
    %

    persistent type_names

    if isempty (type_names)
        type_names = {'parallel flow', 'counter flow', 'single shell pass', ...
                      'multiple shell passes', 'cross flow both unmixed', ...
                      'cross flow Cmax unmixed', 'cross flow Cmin unmixed'};
    end

    % default type: parallel flow
    if nargin == 2
        N = 1;
    elseif isscalar (type)
        N = type;
    else
        [~, N] = ismember (type, type_names);
    end

    % number of pass in multiple shell passes: n == 1
    if N == 4 && (nargin < 4 || n < 2)
        N = 3;
    end

    % test for 0 <= eps <= 1
    if any (eps(:) > 1) || any (eps(:) < 0)
        error ('/teps must be >= 0 and <= 1');
    end

    if any (Cr(:) > 1) || any (Cr(:) < 0)
        error ('/tCr must be <= 1 and <= 0');
    end

    % Set the same size to eps and Cr
    if isscalar (Cr)
        Cr = zeros (size (eps)) + Cr;
    elseif isscalar (eps)
        eps = zeros (size (Cr)) + eps;
    end

    NUT = NaN (size (eps));

    switch N

        % Parallel flow
        case 1

            Cr = 1 + Cr;
            eps_lim = 1 ./ Cr;
            L = eps < eps_lim;
            if any (~L(:))
                warning ('There is any eps >= eps_lim to Cr in %s %s', ...
                         type_names {1}, 'heat exchanger. The result is NaN.'); 
            end

            NUT(L) = - log (1 - eps(L) .* Cr(L)) ./ Cr(L);

        % Counter flow
        case 2

            Le = eps < 1;
            if any (~Le(:))
                warning ('There is any eps >= eps_lim to Cr in %s %s', ...
                         type_names {2}, 'heat exchanger. The result is NaN.'); 
            end

            L = Le & Cr < 1;
            NUT(L) = log ((eps(L) - 1) ./ (eps(L) .* Cr(L) - 1)) ./ (Cr(L) - 1);

            L = Le & Cr == 1;
            NUT(L) = eps(L) ./ (1 - eps(L));

        % Single shell pass
        case 3

            Cr2 = sqrt (1 + Cr .* Cr);
            eps_lim = 2 ./ (1 + Cr + Cr2);
            Le = eps < eps_lim;
            if any (~Le(:))
                warning ('There is any eps >= eps_lim to Cr in %s %s', ...
                         type_names {3}, 'heat exchanger. The result is NaN.'); 
            end

            L = Le & eps > 0;
            eps(L) = (2 ./ eps(L) - 1 - Cr(L)) ./ Cr2(L);
            NUT(L) = log ((eps(L) + 1) ./ (eps(L) - 1)) ./ Cr2(L);

            L = Le & eps == 0;
            NUT(L) = 0;

        % Multiple shell pass
        case 4

            L = eps > 0 & Cr < 1;
            eps(L) = ((eps(L) .* Cr(L) - 1) ./ (eps(L) - 1)) .^ (1 / n);
            eps(L) = (eps(L) - 1) ./ (eps(L) - Cr(L));

            L = eps > 0 & Cr == 1;
            eps(L) = 1 ./ (1 + n * (1 ./ eps(L) - 1));

            Cr2 = sqrt (1 + Cr .* Cr);
            eps_lim = 2 ./ (1 + Cr + Cr2);

            Le = eps < eps_lim;
            if any (~Le(:))
                warning ('There is any eps >= eps_lim to Cr in %s %s', ...
                         type_names {4}, 'heat exchanger. The result is NaN.'); 
            end

            L = Le & eps > 0;
            eps(L) = (2 ./ eps(L) - 1 - Cr(L)) ./ Cr2(L);
            NUT(L) = log ((eps(L) + 1) ./ (eps(L) - 1)) ./ Cr2(L);
            L = Le & eps == 0;
            NUT(L) = 0;

        % Cross flow both unmixed
        case 5

            Le = eps < 1;
            if any (~Le(:))
                warning ('There is any eps >= eps_lim to Cr in %s %s', ...
                         type_names {5}, 'heat exchanger. The result is NaN.'); 
            end

            L = eps > 0 & Le & Cr > 0;

            yo = - Cr(L) .^ (9/7) .* log (1 - eps(L));
            x1 = zeros (size (yo));
            x2 = zeros (size (yo)) + 128;
            x = zeros (size (yo)) + 64;

            err_max = 1.0E0-10;
            i_max = 50;
            err = 1.0;
            i = 0;

            while max (max (err)) > err_max & i < i_max

                y = x .^ (2/7) .* (1 - exp (- x));

                L1 = y < yo;
                x1(L1) = x(L1);
                L1 = y > yo;
                x2(L1) = x(L1);
                x = (x1 + x2) / 2;

                err = (x2 -  x1) ./ x;
                i = i + 1;

            end

            NUT(L) = (x ./ Cr(L)) .^ (9/7);
            L = eps == 0;
            NUT (L) = 0;
            L = Le & Cr == 0;
            NUT (L) = - log (1 - eps(L));

        % Cross flow Cmax unmixed
        case 6

            eps_lim = ones (size (Cr));
            L = Cr > 0;
            eps_lim(L) = 1 - exp (- 1 ./ Cr(L));
            Le = eps < eps_lim;
            if any (~Le(:))
                warning ('There is any eps >= eps_lim to Cr in %s %s', ...
                         type_names {6}, 'heat exchanger. The result is NaN.'); 
            end

            L = Le & L;
            NUT(L) = - log (1 + Cr(L) .* log (1 - eps(L))) ./ Cr(L);

            L = Le & Cr == 0;
            NUT(L) = - log (1 - eps(L));

        % Cross flow Cmin unmixed
        case 7

            eps_lim = ones (size (Cr));
            L = Cr > 0;
            eps_lim(L) = (1 - exp (- Cr(L))) ./ Cr(L);
            Le = eps < eps_lim;
            if any (~Le(:))
                warning ('There is any eps >= eps_lim to Cr in %s %s', ...
                         type_names {7}, 'heat exchanger. The result is NaN.'); 
            end

            L = Le & L;
            NUT(L) = - log (1 + log (1 - eps(L) .* Cr(L)) ./ Cr(L));

            L = Le & Cr == 0;
            NUT(L) = - log (1 - eps(L));

        % No valid type
        otherwise

            str = sprintf ('\n\t''%s''', type_names{:});
            error ('There is not ''%s'' heat exchange. Valid types:%s.', ...
                    type, str);

    end

end
