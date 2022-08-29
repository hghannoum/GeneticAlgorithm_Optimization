function out= runGA(problem, params, choice)

%problem
CostFunction = problem.CostFunction;
nVar = problem.nVar;



%params extraction
MaxIt = params.MaxIt;
nPop = params.nPop;
beta=params.beta
pC = params.pC; %NUMBER OF CHILDREN
nC= round(pC*nPop/2)*2;
mu= params.mu;  %modify one gene out of 10

%template for empty individuals


empty_individual.Position = [];
empty_individual.Cost = [];

%best solution ever found
bestsol.Cost = inf; 

%Initialization
pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    %generate random solution
    pop(i).Position = randi([0,1],1,nVar);

    %evaluate solution
    pop(i).Cost = CostFunction(pop(i).Position);
    %compare solution to best solution 
 
    if pop(i).Cost < bestsol.Cost

        bestsol= pop(i);
    end

end

%best cost of iterations
bestcost = nan(MaxIt, 1);

%main loop
for it= 1:MaxIt

    %selection probabilities
    c=[pop.Cost];
    avgc= mean(c);
    if avgc~=0

    c=c/avgc;
    end

    probs= exp(-beta*c);
%initialization
popc= repmat(empty_individual, nC/2, 2);
%crossover
for k=1:nC/2

    %select parents
   %%randomly    p1 = pop(q1)  p2=pop(q2);
    p1 = pop(RouletteWheelSelection(probs));
    p2 =pop(RouletteWheelSelection(probs));
%perform crossover
    switch(choice)
    case 'SinglePoint'
    [popc(k,1). Position , popc(k,2).Position]= ...
       SinglePointCrossover(p1.Position, p2.Position);
    case 'DoublePoint'
    [popc(k,1). Position , popc(k,2).Position]= ...
       DoublePointCrossover(p1.Position, p2.Position);
    case 'Uniform'
    [popc(k,1). Position , popc(k,2).Position]= ...
       UniformCrossover(p1.Position, p2.Position);
    otherwise
    [popc(k,1). Position , popc(k,2).Position]= ...
       MyCrossover(p1.Position, p2.Position);
    

    end

    


end

%convert popc to single column matrix

popc = popc(:);

%mutation
for l = 1:nC

    %perform mutation
popc(l).Position = Mutate(popc(l).Position, mu);

%evaluate the solution
popc(l).Cost =CostFunction(popc(l).Position);

  %compare solution to best solution 
 
    if popc(l).Cost < bestsol.Cost

        bestsol= popc(l);
    end

end


       % merge and sort population
    pop= SortPopulation([pop; popc]); %so means sorted by cost

    %remove extra members from the population
    pop =pop(1:nPop);
%update best cost of iteration
bestcost(it) = bestsol.Cost;
    %display iteration iteration
    %disp(['iteration' num2str(it) ':Best Cost =' num2str(bestcost(it))]);
end


% results
out.pop = pop;
out.bestsol= bestsol;
out.bestcost= bestcost;
end