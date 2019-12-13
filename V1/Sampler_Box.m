classdef Sampler_Box <  LinOp
    %% Template : Template for Linear  Operator
    %  Matlab Linear Operator Library 
    %
    % TODO Write your do here
    % Example:
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
	   xStart
	   xStep
	   xSize
	   upsampleRate
	   D
	   options
	   
	   inside
    end
    methods
        function this = Sampler_Box(xStart, xStep, xSize, upsampleRate, varargin)
            this.name ='';             %
            this.issquare = false;      % is your operator square?
            this.iscomplex= false;      % is your operator complex?
            this.sizein = xSize;          % what is the size of the right hand side
			if length(this.sizein) == 1
				this.sizein = [this.sizein 1];
			end
			
            this.sizeout = xSize * upsampleRate;         % what is the size of the left hand side
			if length(this.sizeout) == 1
				this.sizeout = [this.sizeout 1];
			end
			
			this.D = length(xStart);
			
			this.xStart = xStart;
			this.xStep = xStep;
			this.xSize = xSize;
			this.upsampleRate = upsampleRate;
			
			p = inputParser;
			addOptional(p, 'radius', []);
			parse(p,varargin{:});
			this.options = p.Results;
			
        end
        % MANDATORY METHODS
        function y = Apply(this,x)   
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
			y = kron(x,ones(this.upsampleRate));
			y = this.applyWindow(y);
			
			
        end
        function y = Adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
			
			x = this.applyWindow(x);
			y = blockproc(x,[this.upsampleRate this.upsampleRate],@(x) mean(mean(x.data)));
		end
		
		function x = applyWindow(this, x)
%  			if ~isempty(this.options.radius)
%  				if isempty(this.inside)
%  					X = makeGrid(this.xStart, this.xStep/this.upsampleRate, this.xSize*this.upsampleRate);
%  					rSquared = sum(X.^2);
%  					this.inside = rSquared <= this.options.radius^2;
%  				end
%  				x(~this.inside) = 0;
%  			end
            x=x;

		end
		
		
%         FACULTATIVE METHODS
%         function y = Inverse(this,x)
%         end
%         function y = AdjointInverse(this,x)
%         end
%         function y = HtH(this,x) %  Apply the HtH matrix
%         end
%         function y = HHt(this,x) %  Apply the HHt matrix
%         end
%         function y = HtWH(this,x,W) %  Apply the HtWH matrix
%         end
    end
end
