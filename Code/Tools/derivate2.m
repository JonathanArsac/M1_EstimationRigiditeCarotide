function [s,b] = derivate2(u, N, option)
%derivate Least-Squares Strain Estimator.
%  s = derivate(u,N), where u is a vector or a matrix of displacement data, 
%  (N>1) is the lentgth of the kernel for the local strain estimation s.
%  If u is a matrix, the strain is estimated for each column vector of the
%  input matrix
%
%  [s,b] = derivate(), return s and the filter coefficients b. The number of  
%  coefficients N is fixed to 10% of the length of the displacement vector(s)
%  u, as default value. 
%
%  s = derivate(u,N,option), option is:
%  - 'same' compute s after zero padding u so that the output has 
%    the same length than the input, 		
%  - 'valid' (this is the default). No zero paddind is done. s is 
%    estimated only from the initial (valid) points of u. 
%

%  Elastography Toolbox 1.0
%  Ph. Delachartre, March 2002
%  revised october 2004
%
%  Reference:
%  [1] F. Kallel and J. Ophir, 
%      A Least-Squares Strain Estimator for Elastography, 
%      Ultrasonic Imaging, 19, 1997, 195-208. 


% Check input parameters:
if nargin < 1 | nargin > 3, eval('help derivate'); error; end;
if length(u) < 2 | ndims(u) > 2, error('u must be a vector or matrix');
end
if nargin < 3, option = 'valid'; end;% set the default option
[L,J] = size(u);
if L < 2, row = 1; u = u'; L = J; J = 1; end;
if nargin < 2, N = 0.1*L; end;% set the value of N   
if N < 2, error('N must be a positive integer > 1'); end;
N = round(N);% N must be integer
if N > L, error('You either need a shorter kernel or more data'); end
 
switch option
	case 'same'
      u = [zeros(fix(N/2),J); u; zeros(fix(N/2)-(rem(N,2)==0),J)];
      M = L;
   case 'valid'
      M = L-N+1;
   otherwise
      error('unknown option')
end

A = [1:N;ones(1,N)];
H = (A*A')\A;
b = -H(1,:); % coefficients of the differentiator filter
s = filter(b,1,u);
s = s(N:end,:);
