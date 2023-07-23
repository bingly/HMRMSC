function [s] = gsp_filter_synthesis(G, filters, c, param)
%GSP_FILTER_SYNTHESIS Synthesis operator of a gsp filterbank
%   Usage:  s = gsp_filter_synthesis(G, filters, c);
%           s = gsp_filter_synthesis(G, filters, c, param);
%
%   Input parameters:
%         G         : Graph structure.
%         filters   : Set of spectral graph filters.
%         c         : Transform coefficients
%         param     : Optional parameter
%   Output parameters:
%         signal    : sythesis signal
%
%   'gsp_filter_synthesis(G,filters,c)' computes the sythesis
%   operator for coefficient c, where the atoms of the transform 
%   dictionary are generalized translations of each graph spectral filter
%   to each vertex on the graph.
%
%      f = D * c 
%
%   where the columns of D are g_{i,m}=T_i g_m, and T_i is a
%   generalized translation operator applied to each filter 
%   hat{g}_m(cdot).  
%
%   Each column of c is the response of the signal to one filter.
%
%   Example:
%
%         Nf = 4;
%         G = gsp_sensor(30);
%         G = gsp_estimate_lmax(G);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_mexican_hat(G, Nf);  
%         f = zeros(G.N,1);
%         f(1) = 1;
%         f = G.L^2*f;
%         ff = gsp_filter_analysis(G,g,f);
%         f2 = gsp_filter_synthesis(G,g,ff);
%         paramplot.show_edges = 1;
%         figure()
%         subplot(211)
%         gsp_plot_filter(G,g)
%         subplot(223)
%         gsp_plot_signal(G,f,paramplot);
%         subplot(224)
%         gsp_plot_signal(G,f2,paramplot);       
%
%   Additional parameters
%   ---------------------
% 
%    param.method  : Select the method ot be used for the computation. 
%      'exact'     : Exact method using the graph Fourier matrix
%      'cheby'     : Chebyshev polynomial approximation
%      'lanczos'   : Lanczos approximation
%     Default: if the Fourier matrix is present: 'exact' otherwise 'cheby'
%    param.order : Degree of the Chebyshev approximation
%     (default=30). 
%    param.verbose : Verbosity level (0 no log - 1 display warnings)
%     (default 1).   
%
%   See also: gsp_filter_analysis gsp_filter_inverse
% 
%   References:
%     D. K. Hammond, P. Vandergheynst, and R. Gribonval. Wavelets on graphs
%     via spectral graph theory. Appl. Comput. Harmon. Anal., 30(2):129--150,
%     Mar. 2011.
%     
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_filter_synthesis.php

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

% Author: Nathanael Perraudin
% Testing: test_filter
% Date: 19 March 2014

% TODO: Perfect
  
% Read input parameters
if nargin < 4
    param = struct;
end

if iscell(G)
    NG = numel(G);
    s = cell(NG,1);
    for ii = 1:NG
        warning('Check what happen here')
       s{ii} = gsp_filter_synthesis(G{ii}, filters{ii}, c{ii}, param);
%         if iscell(s)
%             c{ii} = gsp_filter_analysis(G{ii}, fi{ii}, s{ii}, param);
%         else
%             c{ii} = gsp_filter_analysis(G{ii}, fi{ii}, s, param);
%         end
    end
    return
end


if isnumeric(filters)
    Nf = size(filters,2);
else    
    Nf = numel(filters);
end

if isfield(param, 'exact')
    warning('param.exact is not used anymore. Please use param.method instead');
    if param.exact
        param.method = 'exact';
    else
        param.method = 'cheby';
    end
end

if ~isfield(param,'method')
    if gsp_check_fourier(G)
        param.method = 'exact';
    else
        param.method = 'cheby';
    end
end

if ~isfield(param,'order'); param.order = 30; end
if ~isfield(param,'verbose'); param.verbose = 1; end

if isfield(param, 'cheb_order')
    param.order = param.cheb_order;
    warning('param.cheb_order is not used anymore. Please use param.order instead');
end



switch param.method
    case 'exact' 

        if ~gsp_check_fourier(G)
            if param.verbose
                warning(['GSP_FILTER_SYNTHESIS: The Fourier matrix is not ',...
                    'available. The function will compute it for you. ',...
                    'However, if you apply many time this function, you ',...
                    'should precompute it using the function: ',...
                    'gsp_compute_fourier_basis']);
            end
            G=gsp_compute_fourier_basis(G);
        end
        if isnumeric(filters)
            fie = filters;
        else
            fie = gsp_filter_evaluate(filters,G.e);
        end
%         Nv = size(c,2);
%         s =zeros(G.N,size(c,2));
%         for ii=1:Nf
%             s = s + G.U * ...
%                 (repmat(fie(:,ii),1,Nv) ...
%                 .* (G.U' * c((1:G.N)+G.N * (ii-1),:)));
%         end

        chat = gsp_gft(G,gsp_vec2mat(c,numel(filters)));
        shat = squeeze(sum(bsxfun(@times, fie, chat), 2));
        s = gsp_igft(G, shat);

    case 'cheby'
        if ~isfield(G,'lmax');
            G = gsp_estimate_lmax(G);
            if param.verbose
                warning(['GSP_FILTER_ANALYSIS: The variable lmax is not ',...
                    'available. The function will compute it for you. ',...
                    'However, if you apply many time this function, you ',...
                    'should precompute it using the function: ',...
                    'gsp_estimate_lmax']);
            end
        end


        cheb_coeffs = gsp_cheby_coeff(G, filters,...
                param.order, param.order +1);    

        s=zeros(G.N,size(c,2));

        for ii=1:Nf
            s = s + gsp_cheby_op(G,cheb_coeffs(:,ii),c((1:G.N)+G.N * (ii-1),:));
        end

    
    case 'lanczos'
        s=zeros(G.N,size(c,2));
        if ~iscell(filters)
            filters = {filters};
        end
        for ii=1:Nf
            s = s + gsp_lanczos_op(G, filters{ii}, c((1:G.N)+G.N * (ii-1),:), param);
        end
   
    otherwise
        error('Unknown method: please select exact, cheby or lanczos');
end

end


