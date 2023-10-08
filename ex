Cálculo del vector PageRank de una red creada de forma aleatoria
function PageRank
0) Re-inicializacion de la sucesión de números aleatorios
%  Hacer, si se desea, una vez por sesion de MATLAB,
%  para cambiar la "semilla" de la sucesion de numeros aleatorios
rng('shuffle');
a) Datos
%  N     = numero de paginas de la red
%  nme   = numero medio de enlaces por pagina
%  alpha = factor de amortiguacion
%  tol   = tolerancia para el test de parada
%
N     = 5000;
% N     = 10;  % para hacer pruebas con matriz pequeña
% nme   = 2;   % para hacer pruebas con matriz pequeña
nme   = 7;  
alpha = 0.98;
tol   = 1.e-15;
b) Construcción de las matrices
%  Creacion de la matriz de adyacencia.
%  Se genera una matriz de numeros aleatorios con probabilidad uniforme.
%  Se crea una matriz logica detectando los elementos menores que nme/N.
%
P = rand(N,N);
L = P < nme/N;
%  Conversion de L a tipo numerico.
%  Los valores son mostrados por MATLAB como 1 y 0, pero NO SON NUMEROS.
%  Cualquier operación aritmetica con datos logicos produce un resultado numerico.
%  Para conseguir una matriz P con 1 donde L tiene true y 0 donde L tiene
%  false, basta con sumarle 0.
P = +L;
clear L
%  La matriz de adyacencia debe tener diagonal = 0, ya que una
%  pagina no se enlaza a si misma. 
for k = 1:N
    if (P(k,k) ~= 0) 
        P(k,k) = 0;
    end
end
%  Construcción de la matriz de probabilidades de transicion.
%  Se calcula la suma de cada columna de P.
%  Se divide cada elemento de la matriz por la suma de su columna.
%  Para evitar la division por cero, se sustituyen los ceros de z por unos.
%  
z = sum(P);  
z = z + (z==0);
P = P * diag(1./z);
c) Transformación de P en una matriz "sparse"
PS = sparse(P);
%  Comprobacion de la cantidad de memoria ocupada por la matriz llena 
%  y por la matriz hueca.
%  Descomentar las ordenes siguientes
%  Bytes de memoria ocupados por la matriz llena
%  whos P
%  Bytes de memoria ocupados por la matriz "sparse"
%  whos PS
%  Visualizacion de la estructura de la matriz "sparse"
%  spy(PS)
%  shg
clear P
d) Preparación de las iteraciones
v = 1/N;
%  Diversas opciones de inicializacion
%
%  1) Inicializacion con un vector de probabilidad uniforme
x0 = v*ones(N,1);
%  2) Inicialización con un vector de probabilidad aleatorio
% x0  = rand(N,1);
% x0  = x0/norm(x0,1);
%  3) Inicializacion con un vector de la base canonica de R^n
% x0  = zeros(N,1);
% x0(1) = 1;
e) Iteraciones
for k = 1:50
    x1 = alpha*(PS*x0);
    gamma = 1 - norm(x1,1);
    x1 = x1 + gamma*v;
    err = norm(x1-x0);
    fprintf(' k=%5i,  error = %20.16f \n', k, err)
    if err < tol
        fprintf('\n')
        fprintf(' El algoritmo converge en %3i iteraciones \n', k)
%       h = histogram(x1);
%       h = histogram(x1,length(x1))
%       [m,Ind] = max(x1)
        return
    end
    x0 = x1;
end
end


function DosParticulas(mA, mB, g, G, y0, yp0, T)

%------------------------------------------------------------------
%   Calculo de las trayectorias de dos particulas
%   que interaccionan
%   y que estan sometidas a la accion de la gravedad
%   La particula A tiene masa m_A y se encuentra inicialmente
%   quieta en el punto (0,0,1000); su posicion es x(t)
%   La particula B tiene masa m_B, su posicion inicial es y0,
%   y su velocidad inicial es yp0; su posicion viene dada por y(t)
%
%  | x'' = - \frac{G m_B}{|x-y|^3} (x-y) - g e_3
%  | y'' =   \frac{G m_A}{|x-y|^3} (x-y) - g e_3
%  | x(0)  = (0, 0, 1000)
%  | x'(0) = (0, 0, 0)
%  | y(0)  = y0
%  | y'(0) = yp0
%------------------------------------------------------------------
%   Argumentos de entrada:
%   mA, mB  : masas respectivas de las particulas A y B
%   g       : aceleracion de la gravedad terrestre
%   G       : constante de proporcionalidad (cte. gravitacion universal)
%   y0, yp0 : posicion y velocidad iniciales de la particula B
%   T       : tiempo final para la resolucion numerica.
%------------------------------------------------------------------
%   Test:
%  Enunciado:
%   DosParticulas(1, 1, 9.8, 3.e4, [100; 0; 1000], [0;  1; 0], 40)
%
%   Distintas velocidades iniciales:
%    DosParticulas(1, 1, 9.8, 1.e4, [100; 0; 1000], [ 0; 1; 0], 40)
%    DosParticulas(1, 1, 9.8, 1.e4, [100; 0; 1000], [-1; 1; 0], 40)
%
%   Menor constante de atraccion:
%    DosParticulas(1, 1, 9.8, 1.e3, [100; 0; 1000], [0;  1; 0], 40)
%
%   Mayor constante de atraccion:
%    DosParticulas(1, 1, 9.8, 3.e4, [100; 0; 1000], [0;  1; 0], 40)
%
%   Distintas masas:
%    DosParticulas(1, 5, 9.8, 1.e4, [100; 0; 1000], [ 0; 1; 0], 40)
%
%   Distinta velocidad inicial de B:
%    DosParticulas(1, 1, 9.8, 1.e4, [100; 0; 1000], [ 0; 10; 10], 50)
%------------------------------------------------------------------

Preparacion

param.g = g;
param.G = G;
param.mA = mA;
param.mB = mB;

x0  = [0; 0; 1000];
xp0 = zeros(3,1);
z0 = [ x0; xp0; y0; yp0 ];

Resolucion numerica

fun = @(t,y) dzdt(t,y,param);
[t,z] = ode45(fun, [0,T], z0 );

% Seleccionamos el indice maximo para el que
% z_3(t) (=x_3(t)) >=0
% z_9(t) (=y_3(t)) >=0
%
finx = find(z(:,3)>=0, 1, 'last');
finy = find(z(:,9)>=0, 1, 'last');
fin  = min(finx,finy);
pos  = 1:fin;

Representaciones graficas

close all

%  Orbitas  en 3D
%
figure(1)
plot3(x0(1), x0(2), x0(3), 'b.', 'MarkerSize',20)
hold on
plot3(y0(1), y0(2), y0(3), 'r.', 'MarkerSize',20)
plot3(z(pos,1), z(pos,2), z(pos,3), 'b');
plot3(z(pos,7), z(pos,8), z(pos,9), 'r');
legend('Particula A', 'Particula B')
xlabel('Eje X')
ylabel('Eje Y')
zlabel('Eje Z')
title(' Orbitas de las dos particulas')
grid on
view(140, 30)
hold off

% Proyeccion de orbitas sobre plano OXY
%  Lo siguiente se podria hacer, simplemente,
%  con las ordenes
%  figure(2)
%  view(2)
%
figure(2)
plot(x0(1), x0(2), 'b.','MarkerSize',20)
hold on
plot(y0(1), y0(2), 'r.','MarkerSize',20)
plot(z(pos,1), z(pos,2),'b');
plot(z(pos,7), z(pos,8),'r');
legend('Trayect. A', 'Trayect. B')
xlabel('X')
ylabel('Y')
hold off
shg


% Trayectorias de la solucion
%  No se ve bien por la diferencia de magnitud entre
%  las dos primeras componentes y la tercera.
%  Hacemos una grafica con dos ejes OY a distintas escalas
%
figure(3)
subplot(2,1,1)
yyaxis left
plot(t(pos), z(pos,1:2))
yyaxis right
plot(t(pos), z(pos,3))
title(' Trayectorias de la particula A')
legend('x_1','x_2','x_3')
%
subplot(2,1,2)
yyaxis left
plot(t(pos), z(pos,7:8))
yyaxis right
plot(t(pos), z(pos,9))
title(' Trayectorias de la particula B')
legend('y_1','y_2','y_3')

end

Funcion auxiliar

function [zp] = dzdt(~,z,param)
%
%  Segundo miembro del SDO (12 ecuaciones)
%

zdif   = z(1:3) - z(7:9);
Gdif3 = param.G/norm(zdif)^3;

zp = zeros(size(z));

zp(1:3) = z(4:6);
zp(4:6) = - Gdif3*param.mB*zdif;
zp(7:9) = z(10:12);
zp(10:12) = Gdif3*param.mA*zdif;

zp(6)  = zp(6)  - param.g;
zp(12) = zp(12) - param.g;

end





function Autofun(n)

%------------------------------------------------------------------
%  Calculo de un autovalor próximo a n^2 y una autofuncion de
%
%  | - y'' + 2*cos(2x)) y = \lambda y
%  | y(0)  = 0
%  | y(pi) = 0
%  Se impone la condicion adicional y'(0) = 1
%  para normalizar la autofuncion
%------------------------------------------------------------------
%  Los autovalores de estos problemas se comportan asintoticamente
%  como n^2 y los ceros de las autofunciones separan los ceros de
%  las autofunciones del siguiente autovalor.
%------------------------------------------------------------------
%  Argumento de entrada:
%     n    : valor para proporcionar la aproximación inicial
%            Se toma n^2 como aprox. inicial del autovalor y
%            sin(n*x) como aproximación de la
%            autofunción.
%            Para elegir una entre todas las autofunciones se
%            impone la condición adicional y'(0) = 1.
%------------------------------------------------------------------

Preparacion

dydx = @(x,y,lambda) [ y(2); (2*cos(2*x)-lambda)*y(1)];

condcon = @(ya, yb, ~) [ya(1); yb(1); ya(2)-1 ];

tguess    = linspace(0, pi, 10);
yfunguess = @(x) [ sin(n*x); n*cos(n*x) ];
solinit = bvpinit(tguess, yfunguess, n^2);

options = bvpset('stats','on', 'RelTol', 1.e-4);

Resolucion numerica

sol = bvp4c(dydx, condcon, solinit, options);

The solution was obtained on a mesh of 83 points.
The maximum residual is  9.952e-05. 
There were 3011 calls to the ODE function. 
There were 66 calls to the BC function. 

Representacion grafica

xx     = linspace(0,pi,150);
yy     = deval(sol, xx, 1);
yyg    = yfunguess(xx);
lambda = sol.parameters;

close all
axis([0,pi,-1.5, 1.5])
hold on
leg1 = ['Autofuncion asoc. a autovalor = ',num2str(lambda)];
leg2 = 'Inicializacion de la solucion';
plot(xx, yyg(1,:),'g', 'LineWidth',1)
plot(xx, yy,      'b', 'LineWidth',2)
lgd = legend(leg1, leg2);
set(lgd, 'FontSize', 14)

hold off
shg

end




function Autofun2(n, funcoef)

%------------------------------------------------------------------
%  Calculo de un autovalor y una autofuncion de
%
%  | -y'' + f(x) y = \lambda y
%  | y(0)  = 0
%  | y(pi) = 0
%  Se impone la condicion adicional y'(0) = 1
%  para normalizar la autofuncion
%------------------------------------------------------------------
%  Los autovalores de estos problemas se comportan asintoticamente
%  como n^2 y los ceros de las autofunciones separan los ceros de
%  las autofunciones del siguiente autovalor.
%------------------------------------------------------------------
%  Argumento de entrada:
%     n       : valor para proporcionar la aproximación inicial
%               Se toma n^2 como aprox. inicial del autovalor y
%               sin( n*x ) como aproximación de la autofunción.
%               Para elegir una entre todas las autofunciones se
%               impone la condición adicional y'(0) = 1.
%     funcoef : un manipulador de la funcion f(x) que multiplica a y
%------------------------------------------------------------------
%

Preparacion

dydx = @(x,y,lambda) [ y(2); (funcoef(x)-lambda)*y(1)];

condcon = @(ya, yb, ~) [ya(1); yb(1); ya(2)-1 ];

tguess    = linspace(0, pi, 10);
yfunguess = @(x) [ sin(n*x); n*cos(n*x) ];
solinit = bvpinit(tguess, yfunguess, n^2);

options = bvpset('stats','on', 'RelTol', 1.e-4);

Resolucion

sol = bvp4c(dydx, condcon, solinit, options);

%  Representacion grafica
%
xx = linspace(0,pi,150);
yy = deval(sol, xx, 1);
yyg = yfunguess(xx);
lambda = sol.parameters;

close all
axis([0,pi,-1.5, 1.5])
hold on

leg1 = [' Autofuncion asoc. a autovalor= ',num2str(lambda)];
leg2 =  ' Inicializacion de la solucion ';
plot(xx, yy,      'b', 'LineWidth',2)
plot(xx, yyg(1,:),'g', 'LineWidth',1)

lgd = legend(leg1, leg2);
set(lgd, 'FontSize', 14)

hold off
shg

end






