classdef lRemainder
    %lRemainder  Loi reliquat dérivé d'une distribution de points à partir
    %d'une loi cible
    %   Detailed explanation goes here
    
    properties
        L
        points
        cross;
        dps
        ws
    end
    
    methods
        
        function this=lRemainder(Lt,pts)
            
            this.L=Lt;
            npts=length(pts);
            
            tmp=Lt.extrema();
            [this.points,ind]=sort([pts tmp] );
     
         
            this.cross=find(ind>npts); 
    
            [this.points,fromA,toA]=unique(this.points);
            this.cross=toA(this.cross);

            n=length(this.points);
            
            this.dps=zeros(1,length(this.points)-1);
            incr=1;
            pos=1;
            ndom=length(this.cross);
            extcross=[1 this.cross n ];
            for d=1:(ndom+1)
                for i=extcross(d):(extcross(d+1)-1)
                    if (incr==1)
                        this.dps(pos)=Lt.pdf(this.points(pos) );
                    else    
                        this.dps(pos)=Lt.pdf(this.points(pos+1) );
                    end
                    pos=pos+1;
                end
                incr=1-incr;
            end
            
            this.ws=zeros(1,n+1);
            this.ws(1)=Lt.cumulative(this.points(1));
            this.ws(n+1)=1-Lt.cumulative(this.points(n));
            for k=1:(n-1)
                l=this.points(k);
                r=this.points(k+1);
                this.ws(k+1)= this.L.cumulative(r)-this.L.cumulative(l);
                this.ws(k+1)=this.ws(k+1) - this.dps(k) * (r-l );
            end
          
        end
        
        function x=rv(this)
            p=rvFinite(this.ws);
            u=rand(1);
            n=length(this.points);
            Lt=this.L;
            if (p==1)
                r=this.points(1);
                diff=Lt.cumulative(r);
               
                x=fzero(@(t) (Lt.cumulative(t))/diff-u,r);
                
            elseif (p==n+1)
                l=this.points(n);
                off=Lt.cumulative(l);
                diff=1-off;
                x=fzero(@(t) (Lt.cumulative(t)-off)/diff-u, l);
            else
                l=this.points(p-1);
                r=this.points(p);
                off=Lt.cumulative(l);
                w=this.dps(p-1);
                diff=Lt.cumulative(r ) -off -w*(r-l);
                x=fzero(@(t) (Lt.cumulative(t)-off-w*(t-l))/diff-u, [l r]);
            end
                
        end
        
        function l=find(this,x)
            l=1; r= length(this.points);
            if( x < this.points(l) ) 
                l=-Inf;
            elseif(x>this.points(r) )
                l=r;
            else
        
             while (r-l>1)
             mid=floor((l+r)/2);
               if(x>this.points(mid) ) 
                   l=mid;
               else
                  r=mid;
               end
                
              end
            
            end
            
        end
        
        function y=pdf(this,x)
            np=length(x);
            n=length(this.points);
            Lt=this.L;
            pts=this.points;
            
            w=sum(this.ws);
            y=zeros(np,1);
            for k=1:np
                if(x(k)<pts(1) || x(k) > pts(n) )
                    y(k)=Lt.pdf(x(k))/w;
                else
                 p=this.find(x(k));
                 y(k)=(Lt.pdf(x(k)) - this.dps(p))/w;
                end
            end
        end
        
        function y=cumulative(this,x)
            n=length(this.points);
            Lt=this.L;
            pts=this.points;
            wts=this.ws;
            w=sum(wts);
            
            if(x<pts(1)  )
                y=Lt.cumulative(x)/w;
            elseif (x>pts(n) )
                y=( Lt.cumulative(x)-sum(this.dps.*diff(pts) ))/w;
            else
                p=this.find(x);
                y=sum(wts(1:p) ); 
                xp=pts(p);
                y=y+Lt.cumulative(x)-Lt.cumulative(xp) - (x-xp)*this.dps(p);
                y=y/w;
            end
        end
        
        function y=moments(this,q)
            y=this.L.moments(q);
            tmp=this.points.^(q+1)./(q+1);
            y= y - sum(this.dps.*diff(tmp));
            w=sum(this.ws);
            y=y/w;
        end
        
        function w=weight(this)
            w=sum(this.ws);
        end
        
         function ws=weights(this)
            ws=this.ws;
         end
        
         
    end
    
end

