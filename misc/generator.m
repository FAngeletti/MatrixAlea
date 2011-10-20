classdef generator
    properties
        A
	E
        Laws
	d
    end
    methods
        
        function obj=generator(A,E,Laws)
	   [d,~]=size(E);
	   obj.d=d;
           obj.E=E;
	   obj.A=A;
           obj.Laws=Laws;
        end


	function [source,sink]=initie(obj,n)      
		v=reshape(obj.A.*(obj.E^n),obj.d*obj.d,1);
		n=rvFinite(v);
		source=1+floor((n-1)/obj.d);
		sink=1+mod(n,obj.d);

	end    
        
        function x=rv(obj,n)
            x=zeros(n,1);
            [state,sink]=obj.initie(n);
            CurrentT=obj.E^n;
            for i=1:(n-1)
               CurrentT=obj.E^(n-i);
                newState=rvFinite( obj.E(state,:).* CurrentT(:,sink)' );
                Law=obj.Laws{state,newState};
                x(i)=Law.rv();
                state=newState;
            end
	end
	
  

    end
        
end
