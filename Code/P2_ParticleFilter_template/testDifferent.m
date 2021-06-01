simConst = SimulationConst;
estConst = EstimatorConst;

rng(42, 'twister');

iteration = 0;
seed = [];
epsilon = [];
trackError = [];

while true
    disp("Iteration: " + num2str(iteration));
    seed = [seed, randi(10e7, 1)];
    
    epsilon = [epsilon, rand*0.02];
    simConst.epsilon = epsilon(end);
    estConst.epsilon = epsilon(end);
    
    trackError = [trackError, run(simConst,estConst,false,seed(end))];
    iteration = iteration+1;
end

