function b = absorptionFilter(N,Fs,distance,depth,rho,S,T)

%ABSORPTIONFILTER Computes a sound absorption FIR filter for the ocean
%
%  b = absorptionFilter(N,Fs,distance,depth,rho,S,T)
%    N is the order of the filter
%    Fs is the sampling rate
%    distance is in meters (optional)
%    depth is distance in meters (optional)
%    rho is density in kg/m3 (optional)
%    S is salinity in parts per thousand (optional)
%    T is temperature in deg C (optional)
%    b contains FIR filter coefficients
%
% See absorption for default values and computation details.  The function uses
% absorption to compute the frequency response of the filter and then fir2 to
% generate the filter.
%
% Author: Mandar Chitre
% Last modified: Aug 18, 2004
% $Revision: 1.1 $

% optional parameters
if nargin < 3, distance = []; end
if nargin < 4, depth = []; end
if nargin < 5, rho = []; end
if nargin < 6, S = []; end
if nargin < 7, T = []; end

% compute filter response
f = 0:0.1:1;
for j = 1:length(f)
  freq = f(j)*Fs/2;
  a(j) = absorption(freq,distance,depth,rho,S,T);
end

% compute filter
b = fir2(N,f,a);
