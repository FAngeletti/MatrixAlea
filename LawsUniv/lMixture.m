classdef lMixture
    %LMIXTURE mixture of laws
    %   Detailed explanation goes here
    
    properties
        Laws
        ws
    end
    
    methods
        function this=lMixture(Lts,wts)
            this.Laws=Lts;
            this.ws=wts./sum(wts);
        end
        
        function x=rv(this)
            p=rvFinite(this.ws);
            x=this.Laws{p}.rv();
        end
        
        function y=pdf(this,x)
            n=length(x);
            wts=this.ws;
            Lts=this.Laws;
            y=zeros(1,n);
            for i=1:n
                for k=1:length(wts)
                y(i)=y(i)+  this.ws(k).*Lts{k}.pdf(x(i));
                end
            end
        end
        
        function y=cumulative(this,x)
            wts=this.ws;
            Lts=this.Laws;
            y=0;
            for k=1:length(wts)
                y=y+  wts(k).*Lts{k}.cumulative(x);
            end
        end
        
        function y=moments(this,q)
            y=0;
                        wts=this.ws;
            Lts=this.Laws;
            for k=1:length(wts)
                y=y+  wts(k).*Lts{k}.moments(q);
            end
        end
        
         function y=partialMoments(this,q,x)
            wts=this.ws;
            Lts=this.Laws;
            y=0;
            for k=1:length(wts)
                y=y+  wts(k).*Lts{k}.partialMoments(q,x);
            end
        end
        
    end
    
end

