function q = quadvec(f, varargin)
% QUADVEC calls quad but assumes you can't vectorize your input function.
%   QUADVEC has the same calling syntax as QUAD but doesn't require that
%   the integrand be vectorized.  This makes QUADVEC suitable to perform
%   double-, triple-, and higher order integrals over a non-rectangular
%   domain.
%   Q = QUADVEC(FUN,A,B) tries to approximate the integral of scalar-valued
%   function FUN from A to B to within an error of 1.e-6 using recursive
%   adaptive Simpson quadrature. FUN is a function handle. 
%
%   Example:
%   Compute area of ellipse by computing the result in one quadrant and
%   then multiplying by 4.  Note that the inner integrand has integration
%   limits that are a function of the semi-major and semi-minor ellipse
%   axes, a and b.  
%   First we set up the function handle to the inner integrand, fh.
%   a = 4; b = 3;
%   fh = @(x) quadl(@(y) ones(size(y)), 0, b*sqrt(1-(x/a).^2));
%   ea = 4 * quadvec(fh , 0, a);
%   ea-pi*a*b  % compare integration result to closed-form solution
%
%  Note: This functionality has since been added to MATLAB with the function
%  QUADV.

%   Copyright 2006 The MathWorks, Inc.

  q = quadl(@g, varargin{:}); % like quadl, but supplies g as the argument
  function y = g(X) % make f into a "vectorized" function
    y = zeros(size(X));
    for i = 1:numel(X)
      y(i) = f(X(i)); % this f refers to the argument of quadvec
    end
  end
end
