classdef SDFT <  LinOp
    %% SDFT : Sliced Discrete Fourier operator
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = SDFT(index)
    % Discrete Fourier transform operator
    %
    % FIXME : Should better be merged with DFT
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    %
    %
    % See also LinOp fftn Sfft iSfft ifftn
    
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    properties (SetAccess = protected,GetAccess = public)
        index % index along wich dimension are computed the FFT
        Notindex% ~index
        N % Number of element
    end
    methods
        function this = SDFT(index,sz)
            if (~isempty(index))
                assert(issize(index),'The index should be a conformable  to sz');
                this.index = index;
            else
                this.index = 1:this.ndms;
            end
            
            if nargin>1
            this.sizeout=sz;
            this.sizein=sz;
            end
            this.name ='SDFT';
            this.iscomplex= true;
            this.isinvertible=true;
            this.issquare = true;
        end
        function y = Apply(this,x)
             if (~isempty(this.index))
                dim = 1:ndims(x);
                Iidx = true(ndims(x),1);
                Iidx(this.index) = 0;
                this.Notindex = dim(Iidx);
            else
                this.index = 1:ndims(x);
                this.Notindex = [];
             end
            this.sizein = size(x);
            this.sizeout = size(x);
            this.N= prod(this.sizein(this.index));
            y =  Sfft(x, this.Notindex);
        end
        function y = Adjoint(this,x)
            y =   this.N * iSfft(x,this.Notindex);
        end
        function y = Inverse(this,x)
            y =  iSfft(x, this.Notindex);
        end
        function y = HtH(this,x)
            y = this.N * x;
        end
        function y = AdjointInverse(this,x)
            y =   1/this.N * Sfft(x, this.Notindex);
        end
    end
end
% 
% function Notindex = buildNotindex(Ndims, Index)
% % Build notindex such index+ notindex span the whole dimension Ndims
% Notindex= 1:Ndims;
% if ~isempty(Index)
%     for n = Index
%         Notindex(n==Index)= 0;
%     end
%     Notindex = Notindex(Notindex~=0);
% end
% end
