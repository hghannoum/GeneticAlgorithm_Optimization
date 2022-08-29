function plt=application(choice, nPop, MaxIt, mu)


%%problem definition

problem.CostFunction = @(x) MinOne(x);
problem.nVar = 100;



%%GA parameters(number of iteration, population size)

params.MaxIt=MaxIt;
params.nPop = nPop;
params.pC = 1;
params.mu= mu;
params.beta=1;

%%Run GA 

out = runGA(problem, params, choice);
plt= out.bestcost;
%%results section
end