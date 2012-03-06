% Crée une loi vectorielle de proabibilité jointe 
% P(X) = tr( A (E*P)(X_1)...(E*P)(X_N) tr( E^N)
classdef matrixLawNS 
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
        
        function obj=matrixLawNS(A,E,Laws,n)
            %matrixLawNS(A,E,Laws,n) Crée une loi matricielle non-stationnaire
            %  A : vecteur dual de L : L(M) =tr( A' M);
            %  E : matrice de structure : définit la forme du graphe
            %  markovien
            %  Laws : fonction telle que MP_k= Laws(k) où MP_k représente la matrice de loi au temps k
            %  n  :taille du signal. 
	   [dtmp,nihil]=size(E);
           obj.d=dtmp;
           obj.E=E;
           obj.A=A;
           obj.Laws=Laws;
           obj.n=n;
           obj.norm=obj.projection(E^n);
        end

        
        function x=rv(obj)
        %rv() Generate a random variable of Law L
    %    disp('Begining generation');
	    nt=obj.n;
            x=zeros(nt,1);
            [state,sink]=obj.initie();
          %  CurrentT=obj.E^nt;
          fMLaw=obj.Laws;
            for i=1:(nt)
               CurrentT=obj.E^(nt-i);
                newState=rvFinite( obj.E(state,:).* CurrentT(:,sink)' );
     %           fprintf(1,'Time %d : %d -> %d \n', i, state,newState); 
                MLaw=fMLaw(i);
                Law=MLaw{state,newState};
                x(i)=Law.rv();
                state=newState;
            end
        end
    
        

	
    function y=mvPdf(this,v)
    %L.mvpdf(v) multivariate pdf at vector v 
        m=cell(this.n,1);
        for i=1:this.n
            m{i}=this.matInterm(@(L) L.pdf(v(i)),i );
        end
        y=this.matproj(1:this.n,m);
    end
    
    function y=mvPartialPdf(this,pos,v)
    %L.mvpdf(v) partial multivariate pdf at vector v at point pos 
        nv=length(pos);
        m=cell(nv,1);
        for i=1:nv
            m{i}=this.matInterm(@(L) L.pdf(v(i)),pos(i) );
        end
        y=this.matproj(pos,m);
    end
    
    
	function y=pdf(this,x,pos)
    %L.pdf(x,pos) pdf of the k-marginal at point x
		if (nargin<3)
			pos=1;
		end
	
		len=length(x);
		y=zeros(len);
	
		for k=1:len
		m=this.matInterm( @(L) L.pdf(x(k)), pos ); 	
		y(k)= this.matproj(pos,{m});
		end

    end
    
    function y=cumulative(this,x,pos)
    %L.cumulative(x,pos) cdf of the k-marginal at point x
		if (nargin<3)
			pos=1;
		end
	
		len=length(x);
		y=zeros(len);
	
		for k=1:len
		m=this.matInterm( @(L) L.cumulative(x(k)), pos ); 	
		y(k)= this.matproj(pos,{m});
		end

    end

    
    function m=matMq(this,q,pos)
    % Moments matrix of order q at position pos    
 %       if(nargs<3)
 %           pos=1;
 %       end
		m=this.matInterm(@(L) L.moments(q),pos );
		
	end

	function x=moments(this,q, pos)
    % L.moments(q,k) Moments of the k-margina at order q
	  if nargin<3
	    pos=1;
	  end
	    m=this.matMq(q,pos);
	    x=this.matproj(pos,{m} );
	end

      


      function r=mvMoments(this, points, qs )
      % L.mvMoments(points,qs) Calcul les moments à n points à un ordre général
      % points (k) : position du point k
      % qs(k) : l'ordre des moments au point k
      
	  np=length(qs);
      ms=cellgen( @(i) this.matMq(qs(i),points(i)), np);
      r=this.matproj(points, ms);
	      
      end
  
    end
    
    
    methods
    
	function x=average(this, pos )
		if nargin<2
			pos=1;
		end
		x=this.matproj(pos,{this.matMq(1,pos)});
	end

	function x=variance(this, pos )
		if nargin<2
			pos=1;
		end
		x=this.matproj(pos,{this.matMq(2,pos)});
		x=x-this.average(pos)^2;
	end


	function r=corr(this,k,pos)
		
	
		if nargin<3
			pos=1;
		end


		if k==pos
			m=this.matMq(2,pos);
			r=this.matproj(pos,{m});  
			r=r-this.moments(1,pos).^2;
			r=r/(this.variance(pos));
		else
			mk=this.matMq(1,k);
            mpos=this.matMq(1,pos);
			r=this.matproj( sort([pos k]), {mk mpos} );
			r=r-this.moments(1,pos)*this.moments(1,k);
			r=r/sqrt(this.variance(pos)*this.variance(k));
		end
	end


	function v=corrVect(this,pos)
		
		v=zeros(this.n,1);
        if nargin<2
			pos=1;
        end


        
        matAp=this.matMq(1,pos);
		mat2p=this.matMq(2,pos);
        
		moy=this.moments(1,pos);
		var=this.variance(pos);
		for k=1:this.n

        
       

			if k==pos
				r=this.matproj(pos,{mat2p});  
				r=r-moy.^2;
				v(k)=r/var;
            else
                matAk=this.matMq(1,k);
                if(k>pos)
                    r=this.matproj( [pos k], {matAp matAk} );
                else
                    r=this.matproj( [k pos], {matAk matAp} );
                end
                
				r=r-moy*this.average(k);
				v(k)=r/sqrt(var*this.variance(k));
			end
		end
    end
    end

    
    methods (Hidden=true)
	

	function x=projection(this,M)
    %L.projection(M) compute tr(A' M) 
		x=sum(sum(this.A.*M));
    end

    
	function [source,sink]=initie(obj)
    % L.initie() choose an initial state    
		nt=obj.n;      
        dt=obj.d;
		v=reshape(obj.A.*(obj.E^nt),dt^2,1);
      %  disp('Flat init vect'); disp(v); disp('n');
		nt=rvFinite(v);
     %   fprintf(1, 'Chosen flat state %d\n',nt); 
		source=1+fix((nt-1)/dt);
		sink=1+mod(nt-1,dt);

    end
    
    
    
	function x=matproj(this,ks,Ms)
    %L.matptoj(ks,Ms) compute tr ( A' E .. Ms(1) .. Ms(2) ... E ) where Ms(t) appears at position ks(t)    
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

	function m=matInterm(this, f, pos)
    %M=matInterm(this, f, pos) For a function f,  compute the matrix M(i,j)= E(i,j) * f(MP(pos)(i,j))    
			m=zeros(this.d);
			tmp=this.Laws(pos);
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
      
      
    end
        
end
