clear
clc

dim=7;

prop=[0.054,0.407,0.704,0.813,0.704,0.407,0.054]
prop=diag(prop)
//prop=eye(dim,dim);

t3f=0.067*ones(7,7)
t3b=t3f


r3f=0.198*ones(7,7)
r3b=r3f

t1f=diag([15.99,2.1,1.21,1.02,1.1,1.31,2.84])
t1b=t1f
t2b=t1f
t2f=t1f

r1f=diag([1.44,0.189,0.112,0.114,0.207,0.952,15.28])
r1b=r1f
r2b=r1f
r2f=r1f


//Paso inicial: inicializamos para una pura capa

Tf=eye(dim,dim)
Tb=eye(dim,dim)
Rf=eye(dim,dim)
Rb=eye(dim,dim)


//Existe una segunda capa, así que iteramos:

tf=eye(dim,dim)
tb=eye(dim,dim)
rf=eye(dim,dim)
rb=eye(dim,dim)

aux1=inv(eye(dim,dim)-prop*Rb*prop*rf)
aux2=inv(eye(dim,dim)-prop*rf*prop*Rb)




Rf=Rf+Tb*aux2*prop*rf*prop*Tf
Tf=tf*aux1*prop*Tf

Tb=Tb*aux2*prop*tb
Rb=rb+tf*aux1*prop*Rb*prop*tb


// SIGUIENTE CAPA

tf=eye(dim,dim)
tb=eye(dim,dim)
rf=eye(dim,dim)
rb=eye(dim,dim)

aux1=inv(eye(dim,dim)-prop*Rb*prop*rf)
aux2=inv(eye(dim,dim)-prop*rf*prop*Rb)


Rf=Rf+Tb*aux2*prop*rf*prop*Tf
Tf=tf*aux1*prop*Tf

Tb=Tb*aux2*prop*tb
Rb=rb+tf*aux1*prop*Rb*prop*tb

disp("Tf")
disp(Tf)

disp("Tb")
disp(Tb)

disp("Rf")
disp(Rf)

disp("Rb")
disp(Rb)


