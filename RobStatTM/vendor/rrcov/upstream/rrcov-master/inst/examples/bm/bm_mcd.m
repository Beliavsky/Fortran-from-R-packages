%%%%% bm_mcd %%%%%
% Benchmark for LIBRA::mcdcov() on several n and p.
%
%   V.Todorov: 16.08.2004
%
% For each n and p (specified by the arrays <an> and <ap>) a data set 
%   is generated (see the function gendata() for description of the model)
%   and MCD is computed by calling mcdcov(). All defaults are accepted.
%
% Each test (for given n and p) is performed several times (specified
%   by the parameter <nrep>) and the result is averaged.
%   
% Required input argument: 
%    nrep : Number of times the tests are executed
%
function bm_mcd(nrep)
    if(nrep <= 0)
        nrep = 1;
    end
    disp('*** Benchmark for LIBRA: mcdcov() ***')
    disp(sprintf('\nThe results are averaged on nrep = %d runs.\n', nrep))
    disp('       N    P    TIME (sec.)')
    eps = 0.4;
    b = 10;
    an = [ 100 500 1000 10000 50000];
    ap = [2 5 10 20 30];
    ttime = 0;
    for(i = 1:length(an))
        n = an(i);
        for(j = 1:length(ap))
            p = ap(j);
            if(5*p <= n)
                x = gendata(n, p, eps, b);
                tic;
                for(k = 1:nrep)
                    mcdcov(x, 'plots', 0);
                end
                tt = toc;
                ttime = ttime + tt;
                tt=tt/nrep;
                disp(sprintf('%8d %4d %10.2f', n, p, tt))
            end
        end
    end
    disp(sprintf('\nTotal elapsed time (nrep=%d): %8.0f sec.\n', nrep, ttime))
        
%--------------------------------------------------
% gendata() 
% Generates a location contaminated multivariate 
% normal sample of n observations in p dimensions
%    (1-eps)*Np(0,Ip) + eps*Np(m,Ip)
% where 
%    m = (b,b,...,b)
% Defaults: eps=0 and b=10
%
function result = gendata(n, p, eps, b)

    if(eps < 0 | eps >= 0.5)
        error('eps must be in [0,0.5)')
    end

    mu = zeros(p,1);
    sigma = diag(ones(p,1));
    x = mvnrnd(mu, sigma, n);

%   generate now the outliers
    if(eps > 0)
        mubad = b + mu;
        nbad = floor(eps * n);
        xind = randperm(n);
        xind = sort(xind(1:nbad));
        
        xbad = mvnrnd(mubad, sigma, nbad);
        for(i = 1:nbad)
            x(xind(i),:) = xbad(i,:);
        end
    end
    
    result = x;

