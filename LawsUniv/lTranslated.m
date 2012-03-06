classdef lTranslated
    %Create a translated law from an original law
   
    properties
        L
        mu
    end

    
      methods 
          function this=lTranslated(L,mu)
              this.mu=mu; this.L=L;
          end
          
          
        function x=rv(this)
         %rv() Generate a random variable of Law L
            x=this.L.rv()+this.mu;
        end
        
        function y=pdf(this,x)
        %L.pdf(x) pdf at point x
            y=this.L.pdf(x-this.mu);
        end
        
        function y=cumulative(this,x)
        %L.cumulatice(x) cdf at x
            y=L.cumulative(x-this.mu);
        end
        
        function y=moments(this,q)
        % L.moments(q,k) Moments at order q
            y=this.mu^q;
          
            cnp=q;
            for i=1:q
                y= y+ this.mu^(q-i) *this.L.moments(i) * cnp ;
                cnp=cnp*(q-i)/(i+1);
            end
        end
        
        function ext=extrema(this)
        % L.ext() extrema points of the pdf functions
            ext=this.L.extrema()-this.mu;
        end
        
        function partialMoments(this, q,x)
        % L.partialmoments(q,x) partial moments : E (X^q 1_(X<x) ) 
            y=this.mu^q;
            cnp=q;
            for i=1:q
                y= y+ this.mu^(q-i) *this.L.partialMoments(i,x) * cnp ;
                cnp=cnp*(q-i)/(i+1);
            end
        end
        
    end
    
end



