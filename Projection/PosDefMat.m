


d = 6 ; 
rho = 0.4 ; % rho < 1/2

P = eye(d) ; 

P(1,d) = rho ; 
P(d,1) = rho ; 
for k =1:1:d-1
    P(k,k+1) = rho ; 
    P(k+1,k) = rho ; 
end

[V, l] = eig(P) 
inv(P)
QQ = inv(P)  ; 
tmp = sqrt(diag(QQ)); 
tmp=tmp*tmp'
QQ./tmp

%% 
rho1 = 0.6 ; 
rho2 = 0.2 ; 
Q = eye(d) ; 
Q(1,2) = rho1 ; 
Q(1,3) = rho1 ; 
Q(2,4) = rho1 ; 
Q(3,4) = rho1 ; 
Q(3,5) = rho1 ; 
Q(4,6) = rho1 ; 
Q(5,6) = rho1 ; 
Q =(Q+Q')/ 2 

[V, l] = eig(Q) 
QQ = inv(Q)  ; 
tmp = sqrt(diag(QQ)); 
tmp=tmp*tmp'
QQ./tmp


%% 
rho1 = 0.8 ; 
rho2 = 0.6 ; 
Q = eye(d) ; 
Q(1,2) = rho1 ; 
Q(1,3) = rho1 ; 
Q(2,4) = rho1 ; 
Q(3,4) = (rho1+rho2)/2 ; 
Q(3,5) = rho2 ; 
Q(4,6) = rho2 ; 
Q(5,6) = rho2 ; 
Q =(Q+Q')/ 2 

[V, l] = eig(Q) 
inv(Q)
QQ = inv(Q)  ; 
tmp = sqrt(diag(QQ)); 
tmp=tmp*tmp'
QQ./tmp

