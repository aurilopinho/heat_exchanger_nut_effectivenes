function [eps] = heat_exchanger_eps (NUT, Cr, type, n)

    % This program calculate heat exchange effectiveness 
    %  from NUT and heat capacity ratio folowing textbook
    %  Fundamentals of Heat and Mass Transfer (7th ed.) 
    %  by BERGMAN, T. L., LAVINE, A. S., INCROPERA, F. P. 
    %  and DEWITT, D. P.. Table 11.3 pp 724.
    %
    %  Thanks to Professor Nelson Fernandes Inforzato for providential help
    %
    % Input parameters explained:
    %
    %   NUT: number of transfer units (U A / Cmin)
    %        U: overall heat transfer coefficient
    %        A: heat transfer area
    %        Cmin: smaller heat capacity between both fluids
    %
    %   obs: NUT must be >= 0
    %
    %   Cr: heat capacity ratio (Cmax / Cmin)
    %       Cmax: greater heat capacity between both fluids
    %       Cmin: smaller heat capacity between both fluids
    %
    %   obs: Cr must be >= 0 and <= 1
    %
    %   obs: Both, NUT and Cr can be scallar, vector or matrix.
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
    %   eps: heat exchanger effectiveness (q / qmax)
    %        q: actual heat transfer rate for a heat exchanger
    %        qmax: maximum possible heat transfer rate
    %
    % Examples of the program in use:
    %
    %   Input:  NUT = 5;
    %   Input:  Cr = 0.5;
    %   Input:  eps = heat_exchanger_eps (NUT, Cr, 'parallel flow')
    %   Output: eps = 
    %
    %               0.6663
    %
    %   Input:  NUT = [1 5 10];
    %   Input:  Cr = 0.5;
    %   Input:  eps = heat_exchanger_eps (NUT, Cr, 2)
    %   Output: eps =
    %
    %               0.5647    0.9572    0.9966
    %
    %   Input:  NUT = 5;
    %   Input:  Cr = [0 0.5 1];
    %   Input:  eps = heat_exchanger_eps (NUT, Cr, 'single shell pass')
    %   Output: eps =
    %
    %               0.9933    0.7615    0.5854
    %
    %   Input:  NUT = [1 5 10];
    %   Input:  Cr = [0 0.5 1];
    %   Input:  eps = heat_exchanger_eps (NUT, Cr, 4, 2)
    %   Output: eps =
    % 
    %               0.8647    0.9199    0.7388
    %
    %   Input:  NUT = [1 5 10];
    %   Input:  Cr = [0 0.5 1];
    %   Input:  [NUT Cr] = meshgrid (NUT, Cr);
    %   Input:  eps = heat_exchanger_eps (NUT, Cr, 'cross flow both unmixed')
    %   Output: eps =
    % 
    %               0.6321    0.9933    1.0000
    %               0.5448    0.9058    0.9580
    %               0.4685    0.7501    0.8106
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

    % test for NUT > 0
    if any (NUT(:) < 0)
        error ('/tNUT must be >= 0');
    end

    % test for 0 <= Cr <= 1
    if any (Cr(:) > 1) || any (Cr(:) < 0)
        error ('/tCr must be <= 1 and <= 0');
    end

    % Set the same size to NUT and Cr
    if isscalar (Cr)
        Cr = zeros (size (NUT)) + Cr;
    elseif isscalar (NUT)
        NUT = zeros (size (Cr)) + NUT;
    end

    eps = zeros (size (NUT));

    switch N

        % Parallel flow
        case 1

            Cr = 1 + Cr;
            eps = (1 - exp (- NUT .* Cr)) ./ Cr;

        % Counter flow
        case 2

            L = Cr < 1;
            eps(L) = exp (- NUT(L) .* (1 - Cr(L)));
            eps(L) = (1 - eps(L)) ./ (1 - Cr(L) .* eps(L));
            eps(~L) = 1 ./ (1 + 1 ./ NUT(~L));;

        % Single shell pass
        case 3

            Cr2 = sqrt (1 + Cr .* Cr);
            eps = exp (- NUT .* Cr2);
            eps = 2 ./ (1 + Cr + Cr2 .* (1 + eps) ./ (1 - eps));

        % Multiple shell pass
        case 4

            Cr2 = sqrt (1 + Cr .* Cr);
            eps = exp (- NUT .* Cr2);
            eps = 2 ./ (1 + Cr + Cr2 .* (1 + eps) ./ (1 - eps));

            L = Cr < 1;
            eps(L) = ((1 - eps(L) .* Cr(L)) ./ (1 - eps(L))) .^ n;
            eps(L) = (eps(L) - 1) ./ (eps(L) - Cr(L));
            eps(~L) = 1 ./ (1 + (1 ./ eps(~L) - 1) / n);

        % Cross flow both unmixed
        case 5

            L = Cr > 0;
            NUT_ = NUT(L) .^ (1/9);
            NUT2 = NUT_ .* NUT_;
            NUT7 = NUT2 .* NUT2 .* NUT2 .* NUT_;
            eps(L) = 1 - exp (NUT2 .* (exp (- Cr(L) .* NUT7) - 1) ./ Cr(L));
            eps(~L) = 1 - exp (- NUT(~L));

        % Cross flow Cmax unmixed
        case 6

            L = Cr > 0;
            eps(L) = 1 - exp ((exp (- Cr(L) .* NUT(L)) - 1) ./ Cr(L));
            eps(~L) = 1 - exp (- NUT(~L));

        % Cross flow Cmin unmixed
        case 7

            L = Cr > 0;
            eps(L) = 1 - exp (- NUT(L));
            eps(L) = (1 - exp (- Cr(L) .* eps(L))) ./ Cr(L);
            eps(~L) = 1 - exp (- NUT(~L));

        % No valid type
        otherwise

            str = sprintf ('\n\t''%s''', type_names{:});
            error ('There is not ''%s'' heat exchange.\nValid types:%s.', ...
                    type, str);

    end

end
