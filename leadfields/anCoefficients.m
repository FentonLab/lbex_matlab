function [ an ] = anCoefficients( params )
%
%   This routine is an internal one, called by getHeadModel.m. It is not 
%   directly accessed by the user. 
%
%   Function to evaluate the a_n coefficients in the EEG leadfield
%   expansion. Function is used internally, and not to be
%   directly called by user. Typically, set N=K=3
%
% Inputs:
%   params.seriesCutOff: Truncation index of the series. Used to split direct
%        sum & approximation of remainder. Typically, seriesCutOff = 3,4.
%
%   params.tailApproxDegree: Max. degree of polynomial used to approximate the 
%       remainder of the series. Typically, tailApproxDegree = 3,4.
%
%   params.nMax: Should be infinity, in theory. In practice, nMax=300 suggested
%
%   params.radius(i), i=1,2,...,4 must be input normalized such that:
%                      r(1)=1, d=r(2), c=r(3), b=r(4); 1 > d > c > b
%
%   params.conductivity(i), i = 1,2,...,4
%
% Outputs:
%   an(i), i = 1,2,...,K+1


if ~isfield( params, 'tailApproxDegree' ), params.tailApproxDegree = 3; end
if ~isfield( params, 'seriesCutOff' ), params.seriesCutOff = 3; end
if ~isfield( params, 'nMax' ), params.nMax = 300; end
N = params.seriesCutOff; 
K = params.tailApproxDegree;
nMax = params.nMax;

[cn, omegan] = cnCoefficients( params );
if 0, figure, plot(cn(1:50)), end

% For reasons of possible overflow/underflow, don't try to optimize the 
% following. It needs to be looked into...

% Ax = b
nk = [N+1:nMax]' ; %size = nMax-N
cO = cn(nk)./omegan(nk);
rhs=zeros(K+1, 1); tmp=zeros(2*K+1, 1); nkPower=zeros(size(nk));
for k = 1 : K+1
    nkPower = nk.^(k-1);
    rhs(k) = sum( cO .* nkPower ); % K x 1
    tmp(k) = sum( nkPower ./ omegan(nk) );
end

for k = K+2 : K+K+1
    tmp(k) = sum( (nk.^(k-1)) ./ omegan(nk) );
end

A=zeros(K+1);
for k = 1 : K+1
    A(:,k) = tmp(k:k+K);
end
an = A\rhs; % (K+1) x 1


