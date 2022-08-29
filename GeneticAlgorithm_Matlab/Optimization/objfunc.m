function f = objfunc(x)

%% reading the objective function from a file
eqn = fileread('func.txt');
fh = str2func(eqn);
f1 = fh(x(1),x(2));
[x,y]= meshgrid(-10:1:10, -10:1:10);
f2= fh(x, y);
mesh(x, y, f2);
drawnow;
f = 1/(1 + f1) ;                              %% to change it to maximization