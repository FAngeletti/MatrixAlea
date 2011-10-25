% Crée une loi vectorielle de proabibilité jointe 
% P(X) = tr( A (E*P)(X_1)...(E*P)(X_N) tr( E^N)
classdef matrixLaw
    properties
        A
	E
        Laws
	d
	n
	norm
	stoMatMq=cell(1,1);

    end
    methods
        
        function obj=matrixLaw(A,E,Laws,n)
	   [dtmp,nihil]=size(E);
	   obj.d=dtmp;
           obj.E=E;
	   obj.A=A;
           obj.Laws=Laws;
	   obj.n=n;
	   obj.norm=obj.projection(E^n);
        end

	function m=matMq(this,q)

		if(fix(q)==q && q>0 )
			
		if( q>length(this.stoMatMq) || isempty(this.stoMatMq{q}) )
			this.stoMatMq{q}=this.matInterm( @(L) L.moments(q) );
		end 
		m=this.stoMatMq{q};
		else
		warning 'Cache miss : float moments';
		m=this.matInterm(@(L) L.moments(q) );
		end
	end
	

	function set_n(this,n)
		this.n=n;
		this.norm=obj.projection(this.E^n);
	end
	

	function x=projection(this,M)
		x=sum(sum(this.A.*M));
	end

	function [source,sink]=initie(obj)
		nt=obj.n;      
		v=reshape(obj.A.*(obj.E^nt),obj.d*obj.d,1);
		nt=rvFinite(v);
		source=1+floor((nt-1)/obj.d);
		sink=1+mod(nt-1,obj.d);

	end  

  
        
        function x=rv(obj)
	    nt=obj.n;
            x=zeros(nt,1);
            [state,sink]=obj.initie();
          %  CurrentT=obj.E^nt;
            CLaws=obj.Laws;
            for i=1:(nt-1)
               CurrentT=obj.E^(nt-i);
                newState=rvFinite( obj.E(state,:).* CurrentT(:,sink)' );
                Law=CLaws{state,newState};
                x(i)=Law.rv();
                state=newState;
            end
	end


	function x=matproj(this,ks,Ms)
		nt=this.n;
		Et=this.E;
		m=Et^(ks(1)-1)*Ms{1};
		for i=2:length(ks)
			Edk=Et^( ks(i)-ks(i-1)-1 );
			m= m * Edk * Ms{i};
		end
		m= m * Et^(nt-ks(length(ks) ) );
		x=this.projection(m)/this.norm;
	end

	function m=matInterm(this, f)
			m=zeros(this.d);
			tmp=this.Laws;
			Et=this.E;
			for i=1:this.d
			for j=1:this.d
				L=tmp{i,j};
				if (Et(i,j)== 0)
					m(i,j)=0;
				else
					m(i,j)=Et(i,j).*f(L);
				end	
			end
			end
	end

	
	function y=pdf(this,x,pos)
		if (nargin<3)
			pos=1;
		end
	
		len=length(x);
		y=zeros(len);
	
		for k=1:len
		m=this.matInterm( @(L) L.pdf(x(k)) ); 	
		y(k)= this.matproj(pos,{m});
		end

	end


	function x=moments(this,q, pos)
	  if nargin<3
	    pos=1;
	  end
	    m=this.matMq(q);
	    x=this.matproj(pos,{m} );
	end

      
	function x=average(this, pos )
		if nargin<2
			pos=1;
		end
		x=this.matproj(pos,{this.matMq(1)});
	end

	function x=variance(this, pos )
		if nargin<2
			pos=1;
		end
		x=this.matproj(pos,{this.matMq(2)});
		x=x-this.average(pos)^2;
	end


	function r=corr(this,k,pos)
		
	
		if nargin<3
			pos=1;
		end

		if k<pos
			tmp=pos;
			pos=k;
			k=tmp;
		end


		if k==pos
			m=this.matMq(2);
			r=this.matproj(pos,{m});  
			r=r-this.moments(1,pos).^2;
			r=r/(this.variance(pos));
		else
			m=this.matMq(1);
			r=this.matproj( [pos k], {m m} );
			r=r-this.moments(1,pos)*this.moments(1,k);
			r=r/sqrt(this.variance(pos)*this.variance(k));
		end
	end


	function v=corrVect(this,pos)
		
		v=zeros(this.n,1);
		if nargin<2
			pos=1;
		end

		matA=this.matMq(1);
		mat2=this.matMq(2);
		moy=this.moments(1,pos);
		var=this.variance(pos);
		for k=1:this.n


			if k==pos
				r=this.matproj(pos,{mat2});  
				r=r-moy.^2;
				v(k)=r/var;
			else
				r=this.matproj( sort([pos k]), {matA matA} );
				r=r-moy*this.average(k);
				v(k)=r/sqrt(var*this.variance(k));
			end
		end
	end


% Calcul les moments à n points à un ordre général
% orders(k) : l'ordre des moments au point k
% pos (k) : position du point k
      function r=genMoments(this, orders,pos )

	  np=length(orders);
	  ms=cell(np,1);
	  for i=1:np
	      ms{i}=this.matMq(orders(i));
	  end
	r=this.matproj(pos, ms);
	      
      end
  

    end
        
end
