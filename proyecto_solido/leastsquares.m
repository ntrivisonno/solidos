function [a,fminres]=leastsquares(x,y)
format compact
format long
%Function to calculate the sum of residuals for a given p1 and p2
fun = @(a) sum((y -a(1)*(x+a(2)).^a(3)).^2);
aguess=[40,0.2,0.1];
[a,fminres]=fminsearch(fun,aguess);