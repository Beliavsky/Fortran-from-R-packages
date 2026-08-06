%%%%% bm_lts %%%%%
% Benchmark for LIBRA::ltsreg() on several n and p.
%
%   V.Todorov: 16.08.2004
% 
% For each n and p (specified by the arrays <an> and <ap>) a data set 
%   is generated (see the function gendata() for description of the model)
%   and LTS is computed by calling ltsregress(). All defaults are accepted, 
%   except the trimming proportion alpha - it is set to 0,5
%
% Each test (for given n and p) is performed several times (specified
%   by the parameter <nrep>) and the result is averaged.
%   
% Required input argument: 
%    nrep : Number of times the tests are executed
%
function bm_lts(nrep)
    if(nrep <= 0)
        nrep = 1;
    end
    disp('*** Benchmark for LIBRA: ltsreg() ***')
    disp(sprintf('\nThe results are averaged on nrep = %d runs.\n', nrep))
    disp('       N    P    TIME (sec.)')
    eps = 0.4;
    an = [ 100 500 1000 10000 50000];
    ap = [2 3 5 10];
    ttime = 0;
    for(i = 1:length(an))
        n = an(i);
        for(j = 1:length(ap))
            p = ap(j);
            [x, y] = gendata(n, p, eps);
            tic;
            for(k = 1:nrep)
                ltsregres(x, y, 'alpha', 0.50, 'plots', 0);
            end
            tt=toc;
            ttime = ttime + tt;
            tt=tt/nrep;
            disp(sprintf('%8d %4d %10.2f', n, p, tt))
        end
    end
    disp(sprintf('\nTotal elapsed time (nrep=%d): %8.0f sec.\n', nrep, ttime))

%--------------------------------------------------
% gendata() 
% Generates a data set with bad leverage points (outliers in x-space),
% n observations in p dimensions acording to the model:
%
%   yi = Xi1 + Xi2 + ... + ei
%
% where ei - N(0,1) is the error term, Xi,j for j=1...p-1 - N(0,100) are 
% the non-trivial explanatory variables and Xip is the intercept term.
% The outliers in the x-space are introduced by replacing eps percent of
% Xi1 by values distributed as N(100,100).
%
% Defaults: eps=0
%
function [x, y] = gendata(n, p, eps);

    if(eps < 0 | eps >= 0.5)
        error('eps must be in [0,0.5)')
    end

    p = p-1;
    x = randn(n,p);
    y = zeros(n,1);
    for(i = 1:n)
        y(i) = sum(x(i,:)) + 1 + rand(1);
    end
    
    if(eps > 0)
        nbad = floor(eps * n);
        xind = randperm(n);
        xind = sort(xind(1:nbad));

        %   generate nbad observations distributed as 
        %   N(m, s^2) = N(100,100) => = m + s*rand(nbad)
        xbad = 100 + 10*randn(nbad,1);

        for(i = 1:nbad)
            x(xind(i),1) = xbad(i,1);
        end
    end
    
    result = [x, y];
