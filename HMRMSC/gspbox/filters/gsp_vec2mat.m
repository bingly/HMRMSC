function [ d ] = gsp_vec2mat( d,Nf )
%GSP_VEC2MAT vector to matrix transform
%   Usage: d  = gsp_vec2mat( d,Nf );
%
%   Input parameters:
%       d       : Data
%       Nf      : Number of filter
%
%   Ouput parameter
%       d       : Data
%   
%   Reshape the data from the vector form to the matrix form
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_vec2mat.php

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

% TESTING: test_filter

[M,N] = size(d);

d = reshape(d,M/Nf,Nf,N);


end


