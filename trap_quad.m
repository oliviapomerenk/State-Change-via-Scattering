function [x,w] = trap_quad(N,a,b)
x = linspace(a,b,N);
dx = diff(x);
w = ([0,dx] + [dx,0])/2;
x = x';
w = w';
end

