classdef lNormal
	properties
		mu
		sigma
	end
	methods
		function this=lNormal(m,s)
		this.mu=m;
		this.sigma=s;
		end

		function x=rv(this)
			x=normrnd(this.mu,this.sigma);
		end

		function x=pdf(this,x)
			s2=this.sigma.^2;
			x=exp(-(x-this.mu).^2/(2*s2) ) ./ sqrt(2*pi*s2);
		end

		function x=cmoments(this,q)
		  x=0;
		  if mod(q,2)==0
		  x=this.sigma^q;
		  for i=1:2:q
		    x=x*i;
		  end
		  end
		  end
      
		function x=moments(this,q)
		  x=0;
		  cnp=1;
		  for i=0:q
		    x=x+cnp*this.mu^(q-i)*this.cmoments(i);
		    cnp=cnp*(q-i)/(i+1);
		  end
        end
        
        function x=cumulative(this,x)
            x=normcdf(x,this.mu,this.sigma);
        end

        function ex=extrema(this)
            ex=[ this.mu ];
        end

	end
		
end
