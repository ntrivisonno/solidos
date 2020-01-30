if 1  % fijarse q no es def positiva
T=rand(3);
% lo generamos simetrico al tensor
T=T+T'

[v,d] = eig(T)

% lo generamos random y llegamos a uno def positivo
T = [1.9923 0.1848 1.2176; 0.1848 1.9238 0.8219; 1.2176 0.8219 1.7374]
if 0
A=[0.360819   0.221070   0.519906 ; 0.999279   0.623326   0.077459 ;0.121149   0.693774   0.276793] % otro def pos
endif

% descomposicion espectral
desp_espec = d(1,1)*v(:,1)*v(:,1)'+d(2,2)*v(:,2)*v(:,2)'+d(3,3)*v(:,3)*v(:,3)'

%diadic(v) = v(:,1)*v(:,1) % este es el vector diadico si v est'a def como vector columna

U = sqrt(d(1,1))*v(:,1)*v(:,1)'+sqrt(d(2,2))*v(:,2)*v(:,2)'+sqrt(d(3,3))*v(:,3)*v(:,3)'
% sqrtm % funcion standar para calcular UU=T odnde T es def positivo

T-U*U

% c=a(:)*b(:)' producto dyadic by victor
endif

if 0

A = rand(3);det(A);
A =[0.360819   0.221070   0.519906;
   0.999279   0.623326   0.077459;
   0.121149   0.693774   0.276793]%det(A)=0.30497
disp('det(A)')
det(A)

U2 = A' * A; % es U^2
 

[v,d] = eig(U2)
disp('descomposicion espectral')
U = sqrt(d(1,1))*v(:,1)*v(:,1)'+sqrt(d(2,2))*v(:,2)*v(:,2)'+sqrt(d(3,3))*v(:,3)*v(:,3)'
disp('U = sqrtm(U2)')
U = sqrtm(U2) %nos evitamos la descomposicion espectral
V = sqrtm(A*A')

R = A*inv(U);
V = A*R';
endif
