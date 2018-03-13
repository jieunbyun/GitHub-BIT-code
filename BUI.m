function PI = BUI( x,y )
%{
This function evaluates the posterior probability of independence using
reference prior and median discretization.

INPUT:
[x y] : 2-D data

Output: array of posterior probability of independence with bins
sequentially being halved.

Modification from v2: Changed algorithm of discretization
%}

N = length( x );

PI_ = 2; nbin = 4; PI = [];
while ( (PI_ > 1e-3) && (N/nbin^2 > 5) ) || ( nbin < 5 )
    edgex = quantile( x, (1:nbin-1)/nbin );
    edgey = quantile( y, (1:nbin-1)/nbin );
    edgex = [ min(x)-1;edgex'; max(x)+1 ];
    edgey = [ min(y)-1;edgey'; max(y)+1 ];
    counts = histcounts2( x,y,edgex,edgey );
    
    I = nbin; IJ = I^2;
    cij = counts(:); ci = sum(counts,1); cj = sum(counts,2);
    PI_ = -(IJ-1)*log( pi ) - log( cij(end)+1 ) + ...
       sum( gammaln( cij(1:end-1)+.5 ) ) + ...
       sum( gammaln( cumsum( cij(end:-1:2) ) + 1.5 ) ) - sum( gammaln( cumsum( cij(end-1:-1:1) ) + 2 + cij(end) ) );
    PI_ = PI_ - ( -(2*I-2)*log( pi ) - log( ci(end)+1 ) - log( cj(end)+1 ) + ...             
        sum( gammaln( ci(1:end-1)+.5 ) ) + sum( gammaln( cj(1:end-1)+.5 ) ) + ...
        sum( gammaln( cumsum( ci(end:-1:2) ) + 1.5 ) ) + sum( gammaln( cumsum( cj(end:-1:2) ) + 1.5 ) ) - ...
        sum( gammaln( cumsum( ci(end-1:-1:1) ) + 2 + ci(end) ) ) - ...
        sum( gammaln( cumsum( cj(end-1:-1:1) ) + 2 + cj(end) ) ) );
    PI_ = exp(PI_); PI_ = 1/(1+PI_);
    PI = [PI;PI_];
    
    nbin = nbin+1;
end

        