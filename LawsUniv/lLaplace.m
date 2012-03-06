classdef lLaplace
	properties
		sigma
	end
	methods
		function this=lLaplace(s)
		this.sigma=s;
		end

		function x=rv(this)
			x=Laplacernd(this.sigma);
		end

		function x=pdf(this,x)
			s=this.sigma ;
			x=exp(-abs(x)/s)./ (2*s);
		end

% 		function x=cmoments(this,q)
% 		  x=0;
% 		  if mod(q,2)==0
% 		  x=this.sigma^q;
% 		  for i=1:2:q
% 		    x=x*i;
% 		  end
% 		  end
% 		  end
      
		function x=moments(this,q)
          s=this.sigma ;
          if mod(q,2) == 1 
            x= 0 ; 
          else
            x = factorial(q)*s^q ; 
          end

        end
        
        function x=cumulative(this,x)
            s=this.sigma ;
            if x<0 
                x = exp(x/s)/2 ; 
            else
                x = 1-exp(-x/s)/2 ; 
            end
        end

        function ex=extrema(this)
            ex=[0] ;
        end

	end
		
end
