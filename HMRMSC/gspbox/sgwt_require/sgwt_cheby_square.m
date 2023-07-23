% sgwt_cheby_square : Chebyshev coefficients for square of polynomial
%
% function d=sgwt_cheby_square(c)
%
% Inputs :
% c - Chebyshev coefficients for p(x) = sum c(1+k) T_k(x) ; 0<=K<=M
%
% Outputs :
% d - Chebyshev coefficients for p(x)^2 = sum d(1+k) T_k(x) ;
%     0<=k<=2*M
%
%   Url: http://lts2research.epfl.ch/gsp/doc/sgwt_require/sgwt_cheby_square.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond. 
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

function d=sgwt_cheby_square(c)
M=numel(c)-1;
cp=c;
cp(1)=.5*c(1);

% adjust cp so that
% p(x) = sum cp(1+k) T_k(x)
% for all k>=0 (rather than with special case for k=0)
%
% Then formula for dp in terms of cp is simpler.
% Ref: my notes, july 20, 2009
dp=zeros(1,2*M+1);
% keep in mind : due to indexing from 1
% c(1+k) is k'th Chebyshev coefficient 

for m=0:(2*M)
    if (m==0)
        dp(1+m)=dp(1+m)+.5*cp(1)^2;
        for i=0:M
            dp(1+m)=dp(1+m)+.5*cp(1+i)^2;
        end
    elseif (m<=M)
        for i=0:m
            dp(1+m)=dp(1+m)+.5*cp(1+i)*cp(1+m-i);
        end
        for i=0:(M-m)
            dp(1+m)=dp(1+m)+.5*cp(1+i)*cp(1+i+m);
        end
        for i=m:M
            dp(1+m)=dp(1+m)+.5*cp(1+i)*cp(1+i-m);
        end
    else % M<m<=2*M
        for i=(m-M):M
            dp(1+m)=dp(1+m)+.5*cp(1+i)*cp(1+m-i);
        end
    end
end
d=dp;
d(1)=2*dp(1);

