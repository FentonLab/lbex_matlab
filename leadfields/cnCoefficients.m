function [ cn, omegan ] = cnCoefficients( params )
%
%   Function to evaluate the c_n coefficients in the EEG leadfield
%   expansion. Used 2 times: (a) To evaluate cn(1),cn(2),...cn(N), N finite
%   typically N=3, for the leading terms, and (b) the remaining cn,
%   [N+1,\infty] for evaluation of a_n, coefficients used in approximating
%   the series remainder. Function is used internally, and not to be
%   directly called by user. When called for use in (a), set nMax = 3, and
%   when used for (b), set nMax = 300. 
%
% Inputs:
%   params.nMax: Should be infinity, in theory. In practice, nMax=300 suggested
%
%   params.radius(i), i=1,2,...,4 must be input ordered such that:
%                      r(1) > d=r(2) > c=r(3) > b=r(4)
%
%   params.conductivity(i), i = 1,2,...,4
%
% Outputs:
%   cn(i), i=1,2,...,nMax
%   omegan(i), i=1,2,...,nMax. Used for (b) only. 
%
if ~isfield( params, 'nMax' ), params.nMax = 300; end
nMax = params.nMax;
radius = params.radius;
% Beware... Sun's index is in reverse
% scalp is sigma_4 and innermost brain is sigma_1
% we assumed conductivity index to follow the radius array so we need to
% switch the conductivity
conductivity(1) = params.conductivity(4); % brain
conductivity(2) = params.conductivity(3); % cerebrospinal fluid
conductivity(3) = params.conductivity(2); % skull
conductivity(4) = params.conductivity(1); % scalp

% These must be normalised wrt scalpRadius
b = radius(4)/radius(1); c = radius(3)/radius(1); d = radius(2)/radius(1); 

%
% Conductivity ratio's
%
s12 = conductivity(1) / conductivity(2);
s23 = conductivity(2) / conductivity(3);
s34 = conductivity(3) / conductivity(4);

% 
% To be checked for underflow/overflow...
%
bP = zeros( nMax, 1 ); cP = zeros( nMax, 1 ); dP = zeros( nMax, 1 );
b2 = b*b; bP(1) = b*b2; c2 = c*c; cP(1) = c*c2; d2 = d*d; dP(1) = d*d2;
for n = 2 : nMax
    bP( n ) = bP( n-1 )*b2; 
    cP( n ) = cP( n-1 )*c2; 
    dP( n ) = dP( n-1 )*d2;
end

% cn: eqn(2,3) pp 1244
cn = zeros( nMax, 1 );
for n = 1 : nMax
    np1 = n+1;
    
    lambda = dP(n)*( bP(n)*n*(s12-1)*(s23-1)*np1 + cP(n)*(s12*n + np1)*(s23*n + np1) )* ...
                   ( (s34*n + np1) + dP(n)*np1*(s34-1) ) + ...
         np1*cP(n)*( bP(n)*(s12-1)*(s23*np1 + n) + cP(n)*(s12*n + np1)*(s23-1) ) * ...
                   ( n*(s34-1) + dP(n)*(s34*np1 + n) );
    cn(n) = ((2*n + 1)^4)*dP(n)*cP(n)/lambda;

    if n ==1
        fprintf ( 'cn(n) = %s, lambda = %s \n ', cn(n), lambda);
    end
    
%     x1 = dP(n)*( bP(n)*n*(s12-1)*(s23-1)*np1 + cP(n)*(s12*n + np1)*(s23*n + np1) );
%     x2 = ( (s34*n + np1) + dP(n)*np1*(s34-1) ) ;
%     x3 = np1*cP(n)*( bP(n)*(s12-1)*(s23*np1 + n) + cP(n)*(s12*n + np1)*(s23-1) ) ;         
%     x4 = ( n*(s34-1) + dP(n)*(s34*np1 + n) );
%     x5 = x1*x2 + x3*x4;
%     
%     if n ==1
%         fprintf ( 'dP(n) = %s\n , bP(n) = %s\n , s12 = %s\n , s23 = %s\n, np1 = %s\n , cP(n) = %s\n', dP(n) , bP(n) , s12 , s23 , np1 , cP(n) );        
%         fprintf ( 'x1 = %s\n , x2 = %s\n , x3 = %s\n , x4 = %s\n, x5 = %s \n', x1, x2, x3, x4, x5);
%     end
%    
%     cn(n) = ((2*n + 1)^4)*dP(n)*cP(n)/x5;

%     if n ==1
%         fprintf ( 'after cn(n) = %s\n ', cn(n));
%     end
% 
%     if n ==1
%         fprintf ( 'lambda1 = %s lambda01 = %s lambda02 = %s dP(n) = %s', lambda1, lambda01, lambda02, dP(n));
%     end
%     lambda2 = ( (s34*n + np1) + dP(n)*np1*(s34-1) ) ; 
%     if n ==1
%         fprintf ( 'lambda2 = %s\n', lambda2);
%     end
%     lambda3 = np1*cP(n)*( bP(n)*(s12-1)*(s23*np1 + n) + cP(n)*(s12*n + np1)*(s23-1) ) ; 
%     if n ==1
%         fprintf ( 'lambda3 = %s\n', lambda3);
%     end
%     lambda4 = ( n*(s34-1) + dP(n)*(s34*np1 + n) ) ;
%     if n ==1
%         fprintf ( 'lambda4 = %s\n', lambda4);
%     end
%     ll = lambda1 * lambda2 * lambda3 * lambda4 ;
% 
%     lambda = dP(n)*( bP(n)*n*(s12-1)*(s23-1)*np1 + cP(n)*(s12*n + np1)*(s23*n + np1) )* ...
%                    ( (s34*n + np1) + dP(n)*np1*(s34-1) ) + ...
%          np1*cP(n)*( bP(n)*(s12-1)*(s23*np1 + n) + cP(n)*(s12*n + np1)*(s23-1) ) * ...
%                    ( n*(s34-1) + dP(n)*(s34*np1 + n) );
%     cn(n) = ((2*n + 1)^4)*dP(n)*cP(n)/lambda;
    
end

% Omega_n
omegan = b*([1:nMax]') .* ( 2*([1:nMax]') - 1 ) ./ bP;
fprintf ( 'omegan %s', omegan);
