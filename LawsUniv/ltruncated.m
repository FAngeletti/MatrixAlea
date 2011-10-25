classdef ltruncated
    %ltruncated Truncated law generated from a Law L on an intervall I.
    %  e
    
    properties
        I
        L
        w
    end

    
    methods
        function this=ltruncated(L,I)
            this.I=sort(I);
            this.L=L;
            this.w= L.cumulative(this.I(2))- L.cumulative(this.I(1) );
        end
        
        function x=rv(this)
            Lt=this.L;
            It=this.I;
            x=Lt.rv();
            while ( x>It(2) || x<It(1) )
                x=Lt.rv();
            end
        end
        
        function y=pdf(this,x)
            It=this.I;
            if( x<It(2) && x>It(1) )
                y=this.L.pdf(x)/this.w;
            else
                y=0;
            end
        end
        
        function y=moments(this,q)
            y=this.L.partialMoments(q,this.I(2))-this.L.partialMoments(q,this.I(1));
            y=y/this.w;
        end
        
        function y=partialMoments(this,q,x)
             It=this.I;
            if( x<It(1) )
                y=0;
            elseif ( x<It(2) )
                y=this.L.partialMoments(q,x)-this.L.partialMoments(q,this.I(1));
                y=y/this.w;
            else
                this.moments(q);
            end
        end
        
    end
    
end

