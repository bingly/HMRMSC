function X = gsp_itft(G,Xdot)
% GSP_ITFT Compute the inverse Time Fourier Transform of a time-vertex signal
%   Usage:  X = gsp_itft(G,Xhat)
%
%   Input parameters:
%         G      : Time-Vertex graph structure
%         Xdot   : Time-Vertex Time Fourier Coefficients
%   Output parameters:
%         X      : Time-Vertex signal
%   Compute the inverse Time Fourier Transform of a time-vertex signal
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_itft.php

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

% Author : Francesco Grassi
% Date   : September 2016


NFFT = G.jtv.NFFT;

if isempty(NFFT)
    NFFT = size(Xdot,2);
end;
    

normalize = sqrt(NFFT);

switch G.jtv.transform
    case 'dft'
        X = ifft(Xdot,NFFT,2)*normalize;
        if sum(abs(vec(imag(X)))) < 1e-13 *norm(X(:))
            X = real(X);
        end
    case 'dct'
        X = idct(Xdot.',NFFT).';
    otherwise
        error('Unknown transform');
end


end

