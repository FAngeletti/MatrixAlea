addpath /Users/patriceabry/MATLAB/UTILS_STAT/

%Dimension de la matrice
clear all
% close all
fontsize = 18 ; 
linewidth = 1 ; 

d=2 ; 

% p=0.98;
p=0.1;
% p = 0.3
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
% lgen=@(m) (@() normrnd(m,1));
lgen=@(m,s) (@() normrnd(m,s));

% Matrice des lois (P)
Laws=cell(d,d);

% for i=1:d
%     for j=1:d
% 	s=floor(i/2);
% 	if(i==j)
%         (-1)^i+(-1)^s
%         	Laws{i,j}=lgen((-1)^i+(-1)^s);
% 	else 
% 		Laws{i,j}=lgen(0);
% 	end
%     end
% end

% m = [0 0 0 0 ] ; sigm = [1, 1,1, 1] ; 
% m = [0 1 1 2] ; sigm = [1, 0.75, 0.75, 0.6] ; 
% m = [0 2 2 6] ; sigm = [2, 1, 1, 0.5] ; 
% m = [0 2 2 0] ; sigm = [3, 1, 1, 0.5] ; 
m = [0 -2 2 0] ; sigm = [3, 1, 1, 0.5] ; 


%Laws{1,1}=lgen(m(1),sigm(1));
%Laws{1,2}=lgen(m(2),sigm(2));
%Laws{2,1}=lgen(m(3),sigm(3));
%Laws{2,2}=lgen(m(4),sigm(4));

Laws{1,1}=lnormal(m(1),sigm(1));
Laws{1,2}=lnormal(m(2),sigm(2));
Laws{2,1}=lnormal(m(3),sigm(3));
Laws{2,2}=lnormal(m(4),sigm(4));



%Taille du signal
n=4096 ;

% Définition du générateur
g=matrixLaw(A,E,Laws,n);


% Synthèse d'un signal de longueur n
x=g.rv();

figure(1) ; clf ; 
plot(x)

figure(2) ; clf 
  [hh,bh,gh]=hist1d(x,100);
  plot(bh,hh); hold on ; plot(bh,gh,'--') ; grid on ; 
ghtheo=g.pdf(bh);
 % ghtheo = p*normpdf(bh,m(1),sigm(1)) + q*normpdf(bh,m(2),sigm(2)) +  q*normpdf(bh,m(3),sigm(3)) +p*normpdf(bh,m(4),sigm(4)); 
%ghtheo = ghtheo/d ; 
  plot(bh,ghtheo,'r--') ; 


y = xcov(x,'coeff') ; 
figure(3) ; clf 
plot(y(n:end) ) ;  grid on; hold on;


ytheo= g.corrVect();

plot(ytheo,'--r');

figure(4) ; clf 
debut=n;
stride=20;
fin=debut+20;
plot(y(debut:fin)) ; grid on; hold on
plot(ytheo(1:stride), 'r--')

stheo=2*real(fft(ytheo));
stheo=stheo(1:(n/2));

sfreq=1./n:1./n:0.5;
window = 256  ; Fs = 1 ; 
% ydata = spectrum(signal,ny,ny/2) ; freq = [0:1:ny/2]/ny/Te ; ydata = ydata(:,1) ; 
ydata = pwelch(x,window,window/2) ; 
freq = Fs*(0:1/window:1/2); ydata = ydata(:,1) ; 
size(ydata)

H = figure(6) ; clf ; 
plot((freq),(ydata),'k-') ; grid on ;  hold on ; 
plot( sfreq, stheo, 'r--');

H = figure(7) ; clf ; 
plot(log10(freq),log10(ydata),'k-') ; grid on ; % hold on ; 
% axis([0 200 0 20])


