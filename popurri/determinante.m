% TEOREMA -> si una matriz cualquiera (va creo q positiva), se la divide por su determinante elevado a (1/3), y a esa matriz le sacamos el determinante, este ser'a uno
A=rand(3)

A1 = A/det(A)^(1/3)

det(A1)
