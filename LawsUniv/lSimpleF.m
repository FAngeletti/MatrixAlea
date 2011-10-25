classdef lSimpleF
    %lSimpleF Lois à densité de probabilités constantes par morceau
   
    
    properties
        xs
        ws
        dps
    end
    
    methods
        function this=lSimpleF(points, dps)
            this.xs=sort(points);
            larg=diff(points);
            this.ws=(dps.*larg);
            norm=sum(this.ws);
            this.ws=this.ws./norm;
            this.dps=dps./norm ;
        end
        
        function [l,r]=find(this,x)
            l=1; r= length(this.xs);
            if( x < this.xs(l) ) 
                l=-Inf; r=1;
            elseif(x>this.xs(r) )
                l=r; r=Inf;
            else
        
             while (r-l>1)
             mid=floor((l+r)/2);
               if(x>this.xs(mid) ) 
                   l=mid;
               else
                  r=mid;
               end
                
              end
            
            end
            
        end
        
        
        
        function y=pdf(this,x)
            y=0.*x;
            for i=1:length(x)
            if(x(i)<this.xs(1) || x(i)>this.xs(end))
                y(i)=0;
            else
                [l,r]=this.find(x(i));
                 y(i)=this.dps(l);
            end
            end
        end
        
        function y=cumulative(this,x)
             if (x<this.xs(1) )
                y=0; 
             elseif ( x>this.xs(end))
                y=1;
             else
                [l,r]=this.find(x);
                diff=x-this.xs(l);
                y=sum(this.ws(1:(l-1)))+this.dps(l)*diff;                 
             end
        end
        
        function r=rv(this)
            p=rvFinite(this.ws);
            u=rand(1,1);
            r= (1-u)*this.xs(p)+u*this.xs(p+1);
        end
        
        function y=moments(this,q)
            tmp=this.xs;
            tmp= tmp.^(q+1)./(q+1);       
            y= sum ( diff(tmp) .* this.dps); 
        end
        
        function y=partialMoments(this,q,x)
             if (x<this.xs(1) )
                y=0; 
             elseif ( x>this.xs(end))
                y=this.moments(q);
             else
                [l,r]=this.find(x);
                tmp=([this.xs(1:l) x].^(q+1))./(q+1);  
                y=sum(this.ws(1:r).*diff(tmp) );                 
             end
        end
        
        
        
    end
    
end

