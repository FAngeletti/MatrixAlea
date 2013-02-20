classdef lstudent
	properties
		mu
		sigma
        nu
	end
	methods
		function this=lstudent(m,s,n)
		this.mu=m;
		this.sigma=s;
        this.nu=n;
		end

		function x=rv(this)
            x = random('T',this.nu); x = x*this.sigma + this.mu ; 
			% x=normrnd(this.mu,this.sigma);
		end

		function y=pdf(this,x)
            x = (x-this.mu)/this.sigma ;
			y = pdf('T',x,this.nu)/this.sigma ; 
		end

	function x=cmoments(this,q)
		  x=0;
          if mod(q,2) == 0
 %         tmp = 1 ; 
 %         for i = 1:2:(q-1)
 %             tmp = tmp*i/(this.nu-i-1) ; 
 %         end
 %         tmp = tmp*this.nu^(q/2)*this.sigma^q ; 
 %         x = tmp ;
           x=(this.sigma)^q * gamma((q+1)/2.)*gamma((this.nu-q)/2.)*(this.nu)^(q/2);
           x= x/(sqrt(pi)* gamma(this.nu/2.));
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
        
        function y=cumulative(this,x)
  %          x=normcdf(x,this.mu,this.sigma);
            x = (x-this.mu)/this.sigma ;
			y = cdf('T',x,this.nu) ; 
        end

        function ex=extrema(this)
            ex=[ this.mu ];
        end

	end
		
end
