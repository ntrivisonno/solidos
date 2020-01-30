function V = volumen_hexa8(x,y,z,xnod)
 
  v = 0;
  npg = 8;
  for i = 1:npg
      [w, psi, eta, ceda] = WeightAndCoordHexa(npg,i);
      J = jacobiano (x,y,z,psi,eta,ceda,xnod);
      f = det(J);
      v += f*w;
  endfor
  
endfunction
