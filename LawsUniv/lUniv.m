classdef lUniv
    %LUNIV Abstract class from univariate law

    

    
      methods (Abstract=true)
        function x=rv(this)
         %rv() Generate a random variable of Law L
        end
        
        function y=pdf(this,x)
        %L.pdf(x) pdf at point x
        end
        
        function y=cumulative(this,x)
        %L.cumulatice(x) cdf at x
        end
        
        function y=moments(this,q)
        % L.moments(q,k) Moments at order q
        end
        
        function ext=extrema(this)
        % L.ext() extrema points of the pdf functions
        end
        
        function partialMoments(this, q,x)
        % L.partialmoments(q,x) partial moments : E (X^q 1_(X<x) ) 
        end
        
    end
    
end

