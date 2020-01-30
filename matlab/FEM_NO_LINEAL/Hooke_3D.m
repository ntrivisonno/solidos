function C = Hooke_3D(E,nu,E2,E3,nu23,nu13,G12,G23,G13)
  % C = Hooke_3D(E,nu)
  % C = Hooke_3D(E1,nu12,E2,E3,nu23,nu13,G12,G23,G13)
  % Funcion que calcula la matriz C de 6x6, representacion del tensor de modulos
  % elasticos en notacion de VOigt
  % ENTRADAS:
  % Para solidos isotropos:
  %   E: modulo de Young
  %   nu: coeficiente de Poisson
  % Para solidos ortotropos:
  %   Ei: modulos de Young
  %   nuij: coeficientes de Poisson
  %   Gij: modulos de corte

  isotropic = nargin==2 || isempty(E2);
  
  if isotropic
    c = E * (1 - nu) / (1 + nu) / (1 - 2 * nu);
    c1 = nu / (1 - nu) * c;
    c2 = E / (nu + 1.0) / 2;
    
    C(6,6) = c2;
    C(1,1) = c;
    C(2,2) = c;
    C(3,3) = c;
    C(4,4) = c2;
    C(5,5) = c2;
    
    C(2,1) = c1;
    C(3,1) = c1;
    C(3,2) = c1;
    
    C(1,2) = c1;
    C(1,3) = c1;
    C(2,3) = c1;
  else
    E1 = E;
    nu12 = nu;
    nu32 = nu23*E3/E2;
    nu31 = nu13*E3/E1;
    nu21 = nu12*E2/E1;
    d12 = nu12*nu21;
    d13 = nu13*nu31;
    d23 = nu23*nu32;
    H11 = 1-d23;
    H22 = 1-d13;
    H33 = 1-d12;
    
    if H11*H22*H33 < 0
      error('(1-nu23*nu32)*(1-nu13*nu31)*(1-nu12*nu21) < 0')
    end
    Delta = (H33-d13-d23-2*nu21*nu32*nu13)/(E1*E2*E3);
    if Delta < 0
      error('Delta < 0')
    end
    C(6,6) = G13;
    C(5,5) = G23;
    C(4,4) = G12;
    C(1,1) = H11/Delta/E2/E3;
    C(1,2) = (nu12+nu32*nu13)/Delta/E1/E3;
    C(1,3) = (nu13+nu12*nu23)/Delta/E1/E2;
    C(2,1) = C(1,2);
    C(2,2) = H22/Delta/E1/E3;
    C(2,3) = (nu23+nu21*nu13)/Delta/E1/E2;
    C(3,1) = C(1,3);
    C(3,2) = C(2,3);
    C(3,3) = H33/Delta/E1/E2;
  end
end 
