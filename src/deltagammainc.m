function [rho,sigma,strval] = deltagammainc(x,y,mu,p)
%%                  Generalized incomplete Gamma function                  
% 
% Usage: [rho,sigma,strval] = deltagammainc(x,y,mu,p)
% 
% Description : compute (rho,sigma) such as
%
%  rho*exp(sigma) = 
%         I_{x,y}^{mu,p} := integral over [x,y] of s^(p-1)*exp(-mu*s) ds
% 
% for 0 <= x <= y <= infinity, mu ~= 0 and p > 0 (when mu < 0, p must be 
% integer and y must be finite).
%
% Inputs : (x,y,mu,p) can be scalar or higher dimensional arrays
% with same number of elements.
%
%   x  : scalar or array of non-negative double numbers
%   y  : scalar or array of (possibly infinite) numbers such as all(x<=y) 
%   mu : scalar or array of nonzero double numbers
%   p  : scalar or array of positive numbers such as p(k) is integer when
%        mu(k) < 0 (or equivalently, such as all(p(mu<0)==floor(p(mu<0))) 
% 
% Outputs : the outputs (rho,sigma,strval) will have the same size as x. 
% 
%   rho    : scalar or higher dimensional array such as rho(k) is the
%            computed mantissa of the integral I_{x(k),y(k)}^{mu(k),(k)}
%            for 1 <= k <= numel(x)
%   sigma  : scalar or higher dimensional array such as sigma(k) is the
%            computed exponent of the integral I_{x(k),y(k)}^{mu(k),(k)}
%            for 1 <= k <= numel(x)
%   strval : string scalar or array such as strval(k) is the base 10
%            scientific notation of the computed value of
%            I_{x(k),y(k)}^{mu(k),p(k)} for 1 <= k <= numel(x)

%% Check MEX compilation
if(3 ~= exist("deltagammainc_mexinterface","file"))
    dir = fileparts(mfilename('fullpath'));
    cmd = ['mex -R2018a -silent -lm CFLAGS="\$CFLAGS -std=c99" ', sprintf('%s%skernel.c %s%sdeltagammainc_mexinterface.c -output %s%sdeltagammainc_mexinterface',dir,filesep,dir,filesep,dir,filesep)];
    warning("MEX file must be compiled before using this module. Running MEX file compilation now:");
    warning off MATLAB:mex:GccVersion_link
    try
        eval(cmd);
        warning("MEX compilation done.");
    catch
        error("MEX compilation failed");
    end
end

%% Check consistency
if(~all(x>=0))
    error("x must be >= 0");
end
if(~all(x<=y))
    error("x must be less or equal to y");
end
if(~all(mu~=0))
    error("mu = 0 is not allowed")
end
if(~all(p(mu<0) == floor(p(mu<0))))
    error("p must have integer value when mu < 0");
end
if(~all(isfinite(y(mu<0))))
    error("y = +inf is not allowed when mu < 0");
end

%% Core of the module

% compute integral (call the mex interface module)
[rho,sigma] = deltagammainc_mexinterface(x,y,mu,p); 
rho = reshape(rho,size(x)); 
sigma = reshape(sigma,size(x)); 

% format (rho,sigma) into base 10 reprensentation (in string format)
strval = repmat("",size(x)); 
ID = rho==0 & isinf(sigma) & sigma < 0; 
strval(ID) = sprintf("%.17e",0.);
a=rho(~ID).*exp(sigma(~ID));
ID = find(~ID); 
id = a>1e-308 & a<1e308; 
str = split(sprintf("%.17e ",a(id)));
strval(id) = str(1:end-1); 
id = ~id;  
l10 = log(10);
b = sigma(id)/l10;
b0 = floor(b+log10(rho(id)));
a = rho(id).*exp((b-b0)*l10);
c = [a';b0']; 
str = split(sprintf("%.17fe%+03ld ",c(:)));
strval(ID(id)) = str(1:end-1);

end

