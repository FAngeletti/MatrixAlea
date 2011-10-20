classdef lMarginal
    %lUniv Univariate law created  from a component of a multivariate law
    %   ?
    
    properties
        Lm;
        pos;
    end
    
    methods
        function this=lMarginal(Lm,pos)
        %lUniv(Lm,k) Create an univariate law from the k component of a
        %multivariate law Lm
            this.Lm=Lm;
            this.pos=pos;
        end
        
        function x=rv(this)
         %rv() Generate a random variable of Law L
            v=this.Lm.rv();
            x=v(this.pos);
        end
        
        function y=pdf(this,x)
        %L.pdf(x) pdf at point x
            y=this.Lm.pdf(x,this.pos);
        end
        
        function y=cumulative(this,x)
        %L.cumulatice(x) cdf at x
            y=this.Lm.cumulative(x,this.pos);
        end
        
        function y=moments(this,q)
        % L.moments(q,k) Moments at order q
        y=this.Lm.moments(q,this.pos);
        end
        
    end
    
end

