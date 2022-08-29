function y= Mutate(x, mu)

%changing a value inside our chromosme with a probability of mu
flag = (rand(size(x)) < mu);
y=x;
y(flag) = 1-x(flag);


end