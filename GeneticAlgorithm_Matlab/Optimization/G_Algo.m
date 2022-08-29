function [max_x, max_fit, min_fit, avg_fit, timeElapsed]=G_Algo(Pc, Pm, N, gnMax, eqn)
tic
n = 2;                       %% number of variables
nb = 20;                     %% number of bits to represent each variable
Xmax = 0.5;                  %% geometric constraints
Xmin = 0;
Gn_max = gnMax;                 %% maximum number of generations
N= N;
Pm= Pm;
Pc= Pc;

%generating the intial random population
s = zeros(N,n*nb);           % matrix containing the solutions

for i=1:N                    
    for j=1:(n*nb)
        if( rand()<0.5 )
           s(i,j) = 0;
        else
           s(i,j) = 1;
        end
    end
end

Gn = 0;
x = zeros(n,1);
x_temp = zeros(1,n);
avg_fit = zeros(1,Gn_max);
max_fit = zeros(1,Gn_max);
min_fit = zeros(1,Gn_max);
min_x = zeros(Gn_max,n);
max_x = zeros(Gn_max,n);

while Gn < Gn_max            

%% fitness evalauation
fit = zeros(N,1);            %% array containing the fitness value of all solutions

for i=1:N                    
    dv = zeros(n,1);
    for k=1:n
        for j=1:nb
            dv(k,1) = dv(k,1) + power(2,nb-j)*s(i,j+((k-1)*nb) );
        end
        x_temp(1,k) = Xmin + ( (Xmax - Xmin)*dv(k,1) )/( power(2,nb)-1 );
    end

    fit(i) = objfunc(x_temp);
end

%% Reproduction
sum=0;
for i=1:N
    sum = sum + fit(i);
end
K = 1/sum ;

P = zeros(N,1);            %% array containing the probability of selecting the solution
for i=1:N
    P(i) = K*fit(i);
end

RW = zeros(N,1);           %% array contains the Roulette wheel with angle of each element 
for i=1:N
   if(i==1)
       RW(i) = P(i)*360;
   else
       RW(i) = RW(i-1) + P(i)*360;
   end
end

pool = zeros(N,n*nb);      %% matrix containing the mating pool solutions

for i=1:N                  %% loop for mating pool selection
    r = randi(360);
    for j=1:N
        if( r < RW(j) )
            index = j;
            break;
        end
    end
    pool(i,:) = s(j,:); 
end

%% Two point Crossover
par = zeros((N/2),3);    %% N/2 parent pairs with corresponding random numbers for croos over   
temp = randperm(N,N);    %% Generate N random numbers between 1 and N 

for i=1:N/2              %% Loop to find the pair of parents and their random number
   par(i,1) = temp(i);
   par(i,2) = temp(i + (N/2));
   par(i,3) = rand();
end

rr = zeros(1,2);
for i=1:N/2                %% loop for checking crossover probability
    if( par(i,3) <= Pc )
        rr(1,:) = randperm((n*nb-1),2);
        
        if( rr(1,1) < rr(1,2) )
            r1 = rr(1,1);
            r2 = rr(1,2);
        else
            r1 = rr(1,2);
            r2 = rr(1,1);
        end

        temp = zeros(1,(r2-r1));
        temp(1,:) = pool( par(i,1) , r1+1 : r2 );          
        pool( par(i,1) , r1+1 : r2 ) = pool( par(i,2) , r1+1 : r2 );    %%cross over
        pool( par(i,2) , r1+1 : r2 ) = temp(1,:);
    end
end

%% Bitwise mutation
for i=1:N                  %% loop to carry out bitwise mutation
    for j=1:n*nb
       if( rand() <= Pm )
           if( pool(i,j) == 1)
               pool(i,j) = 0;
           else
               pool(i,j) = 1;
           end
       end
    end
end

s = pool;                   %% carrying the mating pool to initial solution of next iteration
Gn = Gn +1;

%% Calculation of average fitness value

avg = 0;                %% contains the average fitness avlue
max = -1000;
min = 1000;

for i=1:N                    
    dv = zeros(n,1);
    for k=1:n
        for j=1:nb
            dv(k) = dv(k) + power(2,nb-j)*s(i,j+((k-1)*nb) );
        end
        x_temp(1,k) = Xmin + ((Xmax - Xmin)*dv(k))/( power(2,nb)-1 );
    end
    
    fitness = objfunc(x_temp);
    
    if( fitness < min )
        min = fitness;
        min_x(Gn,:) = x_temp;
    end
    
    if( fitness > max )
        max = fitness;
        max_x(Gn,:) = x_temp;%scatter3(X,Y,Z,'.')
    end
   %scatter3(0,0,max_x(Gn,1),'MarkerEdgeColor',[0 .5 .5],...
            % 'MarkerFaceColor',[0 .7 .7],...
            % 'LineWidth',1.5);
    %drawnow;
    %disp(x_temp);
    avg = avg + fitness;
end

avg = avg /N ;

avg_fit(Gn) = avg;
max_fit(Gn) = max;
min_fit(Gn) = min;

end                         %% End of while loop

sol = zeros(N,n);           %% matrix haing the final solutions

for i=1:N                   %% loop to calculate final solution    
    dv = zeros(n,1);
    for k=1:n
        for j=1:nb
            dv(k) = dv(k) + power(2,nb-j)*s(i,j+((k-1)*nb) );
        end
        x_temp(1,k) = Xmin + ((Xmax - Xmin)*dv(k))/( power(2,nb)-1 );
        sol(i,k) = x_temp(1,k);
    end
end




%fi = fopen('OUTPUT','w'); 
%fprintf(fi,'\nThe value of ( x1, x2 ) having maximum fitness after each generation  %d\n');

   % for i=1:Gn_max
     %   for j=1:n
       %     fprintf(fi,'%d\t',max_x(i,j));
       % end
       % fprintf(fi,'\n');
   % end
%colormap(customMap);
%axes('position',[0 0 1 1]);
%plot1 = scatter(lon(1),lat(1),30,concentration(1),'o');
%xlim([-10 10]);
%ylim([-10 10]);
%set(gca,'Color','none');
%set(gca,'CLim',[0, 1E-4]);
%for k = 2:gnMax 
     %plot1.XData = max_x(k,1); 
     %plot1.YData = max_x(k,2); 
     %plot1.CData = concentration(k); 
    
     % pause 2/10 second: 
     %pause(0.2)
%end
timeElapsed = toc;
disp(timeElapsed);
for k = 1:gnMax
plot(max_x(k,1),max_x(k,2),'MarkerFaceColor','auto','Marker','diamond','Color',[0.190 .80 0.120]);
    hold on
axis([-10 10 -10 10])
    grid on
    xlabel('x1')
    ylabel('x2')
    pause(0.5)
    
    if k ~= gnMax
        clf
    end
end
    
    