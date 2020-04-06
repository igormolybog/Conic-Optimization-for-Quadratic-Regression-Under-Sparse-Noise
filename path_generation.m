
function process = path_generation(start,len)

dims = size(start);
width = dims(1);
contVariance = 1/len;  % Picked to match fig.4 https://www.ocf.berkeley.edu/~madani/paper/dc_state_estimation.pdf
jumpSize = 1; % jumps are uniformly distributed between -1 and 1
totalJumpsNum = width*5; % on average, 5 jumps per bus

jumpLocationX = randi(width, totalJumpsNum, 1);
jumpLocationY = randi(len, totalJumpsNum, 1);

jumpValues = jumpSize*(rand(totalJumpsNum,1)-1/2)*2;

continuousAdvances = contVariance.*randn(width, len);

for i = 1:totalJumpsNum
    continuousAdvances(jumpLocationX(i), jumpLocationY(i)) = jumpValues(i);
end
contAdvAndJumps = cumsum(continuousAdvances, 2);

process = start + contAdvAndJumps;


weatherDependentBusesNum = 0;% width; % ceil(0.1*width);
weatherDependentBuses = datasample(1:width,weatherDependentBusesNum,'Replace',false);

% Load variabiility example fig.7-11 https://matpower.org/docs/MOST-manual.pdf

P = [0.9928 0.0072 0 0 0 0 0 0 0 0;
0.0140 0.9679 0.0181 0.0000 0 0 0 0 0 0
0 0.0212 0.9560 0.0228 0 0 0 0 0 0
0 0 0.0318 0.9385 0.0296 0.0001 0 0 0 0
0 0 0 0.0374 0.9297 0.0328 0.0001 0 0 0
0 0 0 0 0.0405 0.9187 0.0408 0 0 0
0 0 0 0 0 0.0435 0.9161 0.0403 0.0001 0
0 0 0 0 0 0.0001 0.0459 0.9076 0.0464 0
0 0 0 0 0 0 0.0001 0.0359 0.9376 0.0265
0 0 0 0 0 0 0 0 0.0191 0.9809]; % obtained from Table 7 of https://pdfs.semanticscholar.org/61b0/1e7cd3004cf3c251bac7ce41d17b3d3ec3d4.pdf
weatherTypeNum = size(P,1);
for i = 1:weatherDependentBusesNum
    mc = dtmc(P);
    weather = simulate(mc,len-1);
    process(weatherDependentBuses(i), :) = process(weatherDependentBuses(i),:).*weather.'/weatherTypeNum;
end
