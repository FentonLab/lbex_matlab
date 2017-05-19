function [ p0, p1 ] = getLegendre( x, degree )
%
%   Routine to evaluate Legendre functions of order 0 & 1
%
%   Output:
%   p0(i,d), i=1,2,...length(x); d=1,2,..., degree: Legendre order 0
%   p1(i,d), i=1,2,...length(x); d=1,2,..., degree: Legendre order 1
%
%   Notes: 
%   1. Stability of the recursion: Have compared against matlab routine for
%   orders 0,1, and degree 600. The max errors are O(1e-13).
%   2. In MATLAB, these would correspond to un-normalized functions
%   3. Test program in test folder, called testGetLegendre.m
%   3. To get exactly the same output in matlab, the following commands will
%   work:
%   x=[-0.75, -0.25, 0, .25, .75]; degree = 4;
%   [p0,p1]=getLegendre( x, degree ); ours = [p0(:,end)'; p1(:,end)'];
%   tmp = legendre( degree, x, 'norm' ); theirs = tmp(1:2, :);
%   p0MaxError = max( ours(1,:) - theirs(1,:) )
%   p1MaxError = max( ours(2,:) - theirs(2,:) )

x = x(:);
p0 = zeros( length(x), degree ); p1 = zeros( length(x), degree );
x2 = x.*x;
%
% Legendre polynomials
%
p0( :, 1 ) = x;
p0( :, 2 ) = 1.5 * x2 - 0.5;
p0( :, 3 ) = x .* ( 2.5 * x2 - 1.5 );
%
% Associated Legendre - order 1 - G+R pp 1033
%
p1( :, 1 ) = -sqrt( 1 - x2 );
p1( :, 2 ) = 3 * x .* p1( :, 1 );
p1( :, 3 ) = p1(:,1) .* ( 2.5 * x2 - 1.5 );

%
% Use recursion, G+R pp 1021, 8.731.#4(1)
%
if degree > 3
    for n = 3 : degree-1
        p0( :, n+1 ) = ( ( n+n+1 )*p0( :, n ).*x - n*p0( :, n-1 ) ) / (n+1);
    end
    for n = 3 : degree-1
        p1( :, n+1 ) = ( (n+n+1)*p1( :, n ).*x - (n+1)*p1( :, n-1 ) ) / n;
    end
end

%
% Normalize them
%
%ns = [ 1 : degree ];
%p0 = repmat( sqrt( ns + 0.5 ), length(x), 1 ) .*p0; 
%p1 = -repmat( sqrt( (ns + 0.5).*factorial( ns-1 )./factorial( ns+1 ) ), length(x), 1 ) .* p1;

    
    