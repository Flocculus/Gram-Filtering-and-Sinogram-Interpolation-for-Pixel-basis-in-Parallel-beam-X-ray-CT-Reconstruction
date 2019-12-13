classdef SumLinOp < LinOp
    %% SumLinop : Sum of linear operators
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = SumLinop(ALinOp,alpha)
    % Element wise sum  of LinOps:
    % Sum the all linop contained in vector ALINOP weighted by ALPHA
    % (default 1)
    % Obj  sum_n alpha(n) * ALinOp(n)
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
        ALinOp     % Array of linop
        numLinOp   % number of linop
        alpha      % scalar factor
    end
    
    methods
        function this = SumLinOp(ALinOp,alpha)
            this.name ='SumLinOp';
			
			if nargin == 1
				alpha = 1;
			end
			
			this.numLinOp = numel(ALinOp);
			assert(isnumeric(alpha)&& ( isscalar(alpha) || ( isvector(alpha) && (numel(alpha)== this.numLinOp))),'second input should be a scalar or an array of scalar of the same size as the first input');
			if  isscalar(alpha)
				this.alpha = repmat(alpha, 1, this.numLinOp) ;
			else
				this.alpha = alpha;
			end
			
			allLinOps = all( cellfun(@(x)(isa(x, 'LinOp')), ALinOp) );
			assert(iscell(ALinOp) && allLinOps, 'First input should be a cell array LinOp');
			
			% if any of the inputs is itself a sum, expand it
			eLinOp = []; % expanded LinOp list
			newAlphas = [];
			for i = 1:length(ALinOp)
				if isa(ALinOp{i}, 'SumLinOp')
					eLinOp = [eLinOp ALinOp{i}.ALinOp];
					newAlphas = [newAlphas  ALinOp{i}.alpha];
				else
					eLinOp = [eLinOp ALinOp(i)];
					newAlphas = [newAlphas this.alpha(i)];
				end
			end
			this.alpha = newAlphas;
			ALinOp = eLinOp;
			this.numLinOp = length(ALinOp);
			
			this.ALinOp = ALinOp;
			this.iscomplex= this.ALinOp{1}(1).iscomplex;
            this.issquare = this.ALinOp{1}(1).issquare;
            this.isinvertible=false;
            this.sizein =  this.ALinOp{1}(1).sizein;
            this.sizeout =  this.ALinOp{1}(1).sizeout;
            for n =2:this.numLinOp
                assert(isempty(this.ALinOp{n}(1).sizein)  || isequal(this.sizein,this.ALinOp{n}(1).sizein),'%d-th input does not have the right hand side size ', n) ;
                assert(isempty(this.ALinOp{n}(1).sizeout) ||isequal(this.sizeout,this.ALinOp{n}(1).sizeout),'%d-th input does not have the left hand side size ', n);
                this.iscomplex= this.ALinOp{n}(1).iscomplex || this.iscomplex ;
            end
            
                 
        end
        
        function y = Apply(this,x) % Apply the operator
			LinOp.checkSize(x, this.sizein)
            y = zeros(this.sizeout);
            for n = 1:this.numLinOp
                y = y + this.alpha(n) .* this.ALinOp{n}(1).Apply(x);
            end
        end
        function y = Adjoint(this,x) % Apply the adjoint
			LinOp.checkSize(x, this.sizeout);
            y =  zeros(this.sizein);
            for n = 1:this.numLinOp
                y = y + this.alpha(n) .* this.ALinOp{n}(1).Adjoint(x);
            end
        end
    end
end

