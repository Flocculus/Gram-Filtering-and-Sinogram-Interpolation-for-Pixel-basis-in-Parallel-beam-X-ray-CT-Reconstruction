classdef adjoint < LinOp
    %% adjoint : overload of adjoint function for LinOp
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = adjoint(LinOp)
    % Obj is the adjoint of the LinOp
    %
    %
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
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
        TLinOp     % linop
    end
    
    methods 
        function this = adjoint(TLinOp)
            this.name ='adjoint';
            
            
            assert(isa(TLinOp,'LinOp'),'Input should be a  LinOp');
            this.TLinOp = TLinOp;
            this.iscomplex= this.TLinOp.iscomplex;
            this.isinvertible=this.TLinOp.isinvertible;
            this.issquare = this.TLinOp.issquare;
            this.sizein =  this.TLinOp.sizeout;
            this.sizeout =  this.TLinOp.sizein;
            
          end
        
        function y = Apply(this,x) % Apply the operator
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizein);            
            y =this.TLinOp.Adjoint(x);
        end
        function y = Adjoint(this,x) % Apply the adjoint
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizeout);
            y =this.TLinOp.Apply(x);
        end
        function y = HtH(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizein);
            y =this.TLinOp.HHt(x);
        end
        function y = HHt(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizeout);
            y =this.TLinOp.HtH(x);
        end
        function y = Inverse(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizeout);
            y =this.TLinOp.AdjointInverse(x);
        end
        function y = AdjointInverse(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizein);   
            y =this.TLinOp.AdjointInverse(x); 
        end
    end
end

