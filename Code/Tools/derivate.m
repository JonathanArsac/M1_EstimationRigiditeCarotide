function [s,N] = derivate(u, N, option)
%derivate Least-Squares Strain Estimator.
%  s = derivate(u,N), where u is a vector or a matrix of displacement data, 
%  (N>1) is the lentgth of the kernel for the local strain estimation s.
%  If u is a matrix, the strain is estimated for each column vector of the
%  input matrix
%
%  [s,N] = derivate(), return s and the default value N fixed to  
%  10% of the length of the displacement vector(s) u. 
%
%  s = derivate(u,N,option), option is:
%  - 'same' compute s after zero padding u so that the output has 
%    the same length than the input, 		
%  - 'valid' (this is the default). No zero paddind is done. s is 
%    estimated only from the initial (valid) points of u. 
%

%  Elastography Toolbox 1.0
%  Ph. Delachartre, March 2002
%
%  Reference:
%  [1] F. Kallel and J. Ophir, 
%      A Least-Squares Strain Estimator for Elastography, 
%      Ultrasonic Imaging, 19, 1997, 195-208. 

progress_bar=waitbar(0,'calculating strain');

% Check input parameters:
if nargin < 1 | nargin > 3, eval('help LSstrain'); error; end;
if length(u) < 2 | ndims(u) > 2, error('u must be a vector or matrix');
end
if nargin < 3, option = 'valid'; end;% set the default option
[L,J] = size(u);
if L < 2, row = 1; u = u'; L = J; J = 1; else, row = 0; end;
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

A = [[1:N+M-1];ones(1,N+M-1)]';
s = zeros(M,J);

waitbar(0.1,progress_bar);

for i = 1:M
   Ai = A(i:i+N-1,:);
   H = (Ai'*Ai)\Ai';
   H = H(1,:)'*ones(1,J);
   s(i,:) = sum(H.*u(i:i+N-1,:));
   
   waitbar((0.85*i/M)+0.1,progress_bar);
end

if row
   s = s';
end
   waitbar(1,progress_bar);
   close(progress_bar);