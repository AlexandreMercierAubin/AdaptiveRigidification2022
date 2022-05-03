classdef TriangleMaterial < handle
    %Material that can independently set for each element, and/or reused for
    %multiple of them.
    properties
        rho     % Density
        mu      % Lamé parameter
        lambda  % Lamé parameter
        alpha0  % Rayleigh parameter
        alpha1  % Rayleigh parameter
        strainUpperBound % 5% would be 1.05
        strainLowerBound % 5% would be 0.95
        color            % color of the material
        cacheName = 'defaultCacheName';
    end
    
    methods
        function obj = TriangleMaterial(rho, mu, lambda, alpha0, alpha1, color, strainUpperBound, strainLowerBound)
            if nargin > 0
                obj.rho = rho;
            else
                obj.rho = 100;
            end
            
            if nargin > 1
                obj.mu = mu;
            else
                k = 4e4;
                nu = 0.40;
                mu = k / 2 / (1 - nu); 
                obj.mu = mu;
            end
            
            if nargin > 2
                obj.lambda = lambda;
            else
                k = 4e4;
                nu = 0.40;
                lambda = k * nu / (1 + nu) / (1 - 2 * nu);
                obj.lambda = lambda;
            end
            
            if nargin > 3
                obj.alpha0 = alpha0;
            else
                obj.alpha0 = 0.03;
            end
            
            if nargin > 4
                obj.alpha1 = alpha1;
            else
                obj.alpha1 = 0.02;
            end
            
            if nargin > 5
                obj.color = color;
            else
                obj.color = [0.3,0.46,0.8];
            end
            
            if nargin > 6
                obj.strainUpperBound = strainUpperBound;
            else
                obj.strainUpperBound = Inf;
            end
            
            if nargin > 7
                obj.strainLowerBound = strainLowerBound;
            else
                obj.strainLowerBound = Inf;
            end
        end
    end
end
