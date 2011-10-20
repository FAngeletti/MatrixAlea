addpath /Users/patriceabry/MATLAB/UTILS_STAT/

%Dimension de la matrice
clear all
close all
fontsize = 18 ; 
linewidth = 1 ; 

d=5 ;

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
        (-1)^i+(-1)^s
        	Laws{i,j}=lgen((-1)^i+(-1)^s);
	else 
		Laws{i,j}=lgen(0);
	end
    end
end

% Définition du générateur
g=generator(A,E,Laws);

%Taille du signal
n=20000 ;

% Synthèse d'un signal de longueur n
x=g.rv(n);

figure(1) ; clf ; 
plot(x)

figure(2) ; clf 
  [hh,bh,gh]=hist1d(x,100);
  plot(bh,hh); hold; plot(bh,gh,'--')

window = 512 ; Fs = 1 ; 
% ydata = spectrum(signal,ny,ny/2) ; freq = [0:1:ny/2]/ny/Te ; ydata = ydata(:,1) ; 
ydata = pwelch(x,window,window/2) ; 
freq = Fs*(0:1/window:1/2); ydata = ydata(:,1)*4 ; 
size(ydata)
H = figure(6) ; clf ; 
set(gca,'FontSize',fontsize) ; 
plot((freq),(ydata),'k-','LineWidth',linewidth) ; grid on ; % hold on ; 

H = figure(7) ; clf ; 
set(gca,'FontSize',fontsize) ; 
plot(log10(freq),log10(ydata),'k-','LineWidth',linewidth) ; grid on ; % hold on ; 
% axis([0 200 0 20])

y = xcov(x,'coeff') ; 
figure(8) ; clf 
plot(y) ; grid on
