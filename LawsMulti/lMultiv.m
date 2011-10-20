classdef lMultiv
    %lMultiv Abstract class for multivariate law
    
      methods (Abstract=true)        
        function x=rv(this)
         %rv() Generate a random variable of Law L
        end
        
        function y=pdf(this,x,k)
        %L.pdf(x,pos) pdf of the k-marginal at point x
        end
        
        function y=cumulative(this,x,k)
        %L.cumulative(x) cdf of the k-marginal at x
        end
        
        function y=moments(this,q,k)
        % L.moments(q,k) Moments of the k-margina at order q
        end
        
        function y=mvPdf(this,v)
        %L.mvpdf(v) multivariate pdf at vector v
        end
        
        function y=mvCumulative(this,v)
        %Rectangular cumulative at vector v
        end
        
        function y=mvMoments(this, points, qs)
        %L.mvMoments(points, qs) Multipoints points moments E prod_k [ X_points(k)^qs(k) ]
        end
        
    end
    
end

