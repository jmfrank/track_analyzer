classdef counter < handle
    
    properties (SetAccess = public )
        
        n
        
        f
        
        sum
        
        sum_sq
        
        std
        
        sem
        
        mean
    end
    
    methods
        
        % Constructer
        function obj = counter(n, f)
            
            if nargin < 2
                f = 1:n;
            end
            
            obj.sum = zeros(1,n);
            obj.sum_sq = zeros(1,n);
            obj.n=zeros(1,n);
            obj.f=f;
        end
        
        % Sum and sum sq. 
        function add(obj, D, F )
            
            if length(size(D)) > 2 || ~any(size(D)==1)
                error(' Input dimension is incorrect')
            end
            
            % Reshape input data if necessary. 
            if size(D,1)~=1
                D = D';
            end
            
            % Figure out idx values. 
            if nargin < 3
                F = 1:length(D);
            end
            
            % Find when F == f
            [~, idx_obj, idx_d] = intersect(obj.f,F);
           
            obj.sum(idx_obj) = obj.sum(idx_obj) + D(idx_d);
            
            obj.sum_sq(idx_obj) = obj.sum_sq(idx_obj) + D(idx_d).^2;
            
            % Count 
            obj.n(idx_obj) = obj.n(idx_obj) + 1;
        end
    
        % Std and sem. 
        function stats(obj)
           
            obj.std = sqrt( ( obj.n.*obj.sum_sq - obj.sum.^2 ) ./ (obj.n.*(obj.n-1)) );
            
            obj.sem = obj.std ./ sqrt( obj.n );
            
            obj.mean = obj.sum ./ obj.n;
            
        end
    end
    
end
