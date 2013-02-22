classdef lGamma
	properties
		alpha
		beta
	end
	methods
		function this=lGamma(alpha,beta)
		this.alpha=alpha;
		this.beta=beta;
		end

		function x=rv(this)
			x=gamrnd(this.alpha,this.beta);
		end

		function y=pdf(this,x)
			y=gampdf(x,this.alpha,this.beta) ;
        	end     
      
		function x=moments(this,q)
		  x=this.beta^q;
                  x=x*gamma(this.alpha+q)/gamma(this.alpha);
        end
        
        function y=cumulative(this,x)
            if(x>=0)
                y= gammainc(x/this.beta,this.alpha);
            else
                y=0;
            end
        end
        
        function y=partialMoments(this,q,x)
            y= ( gamma(this.alpha+q)/gamma(this.alpha) )*gammainc(x/this.beta,this.alpha+q)*(this.beta)^q;
        end
        
        function mds=extrema(this)
            mds=(this.alpha-1)*this.beta;
        end
        
    end
    
		
end
