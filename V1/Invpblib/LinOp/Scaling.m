classdef Scaling <  LinOp
    %% Scaling : scaling operator
    %  Matlab Linear Operator Library 
    %
    % Example
    % Obj = Scaling(scale)
    %
    % Build the scaling operator that scale the input by the scalar factor
    % scale
    %
    % Please refer to the LinOp superclass for documentation
    % See also LinOp
        
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
        scale % scale factor
    end
    methods
        function this = Scaling(scale,sz)
            this.name ='Scaling';
            this.issquare = true;
            if isscalar(scale)
                this.scale = scale;
                if isreal(scale)
                    this.iscomplex=false;
                else
                    this.iscomplex=true;
                end
                if scale == 0;
                    this.isinvertible= false;
                end
            else
                error('Scale value must be a scalar');
            end
            
            if nargin>1
            this.sizeout=sz;
            this.sizein=sz;
            end
            
        end
        function y = Apply(this,x)
            this.sizeout=size(x);
            this.sizein=size(x);
            y =this.scale .* x;
        end
        function y = Adjoint(this,x)
            this.sizeout=size(x);
            this.sizein=size(x);
            if this.iscomplex
                y =this.scale .*x;
            else
                y =conj(this.scale) .*x;
            end
        end
        function y = Inverse(this,x)
            if ( ~this.isinvertible)
                error('Operator non invertible');
            end
            y =(1./this.scale) .*x;
        end
        function y = AdjointInverse(this,x)
            if (~this.isinvertible)
                error('Operator non invertible');
            end
            if this.iscomplex
                y = (1./this.scale) .*x;
            else
                y =conj(1./this.scale) .*x;
            end
        end
    end
end
