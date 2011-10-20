%Dimension de la matrice
d=5;

p=0.9;
q=1-p;

J=zeros(d);
for(i=1:(d-1))
J(i,i+1 )=1;
end 
J(d,1)=1;

Id=zeros(d);
for(i=1:d)
Id(i,i)=1;
end 


E= p.*Id+q.*J ;

%Définition de l'opérateur de projection L(M) = <A,M> 
A=ones(d);

% lgen : générateur de génrateur aléatoire 
lgen=@(m) (@() normrnd(m,1));

% Matrice des lois (P)
Laws=cell(d,d);

for i=1:d
    for j=1:d
	s=floor(i/2);
	if(i==j)
        	Laws{i,j}=lgen((-1)^i+(-1)^s);
	else 
		Laws{i,j}=lgen(0);
	end
    end
end


% Définition du générateur
g=generator(A,E,Laws);

%Taille du signal
n=1000;

% Synthèse d'un signal de longueur n
x=g.rv(n);

plot(x)
