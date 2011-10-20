classdef lVectorized
    %LVECTORIZED Create a multivariate law by repeating independant
    %realisation of a univariate law
    %   Detailed explanation goes here
    
    properties
        Luni
        n
    end
    
    methods
        function this=lVectorized(L,n)
        %lVectorized(L,n) create a vectorized law of size n with univariate
        %law L
        this.Luni=L;
        this.n=n;
        end
        
        
        function v=rv(this)
        %L.rv() Generate a vector of law L;
            v=vRandL(this.Luni,this.n);
        end
        
        function y=pdf(this,x,k)
        %L.pdf(x,k) pdf of the k-marginal at point x
            y=this.Luni.pdf(x);
        end
        
        function y=cumulative(this,x,k)
        %L.cumulatice(x,k) cdf of the k-marginal at x
            y=this.Luni.cumulative(x);
        end
        
        function y=mvPdf(this,v)
        %L.mvpdf(v) multivariate pdf at vector v
            y=  prod(this.Luni.pdf(v));
        end
        
        function y=moments(this,q,k)
        % L.moments(q,k) Moments of the k-marginal at order q
        y=this.Luni.moments(q);
        end
        
        function y=mvMoments(this, points, qs)
        %L.mvMoments(points, qs) Multipoints points moments E prod_k [ X_points(k)^qs(k) ]
            y=1;
            for q=qs
                y=y* this.Luni.moments(q);
            end
        end
        
        
    end
    
end

