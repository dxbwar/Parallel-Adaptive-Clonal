function [pop ,rng] = initializeCS(num, bounds, options)
% function [pop]=initializega(populationSize, variableBounds,
%                           evalOps,options)
%    initializega creates a matrix of random numbers with 
%    a number of rows equal to the populationSize and a number
%    columns equal to the number of variables
%
% pop            - the initial, random population 
% populatoinSize - the size of the population, i.e. the number to create
% variableBounds - a matrix which contains the bounds of each variable, i.e.
%                  [var1_high var1_low; var2_high var2_low; ....]
% options        - options to the initialize function, ie. 
%                  [type prec] where eps is the epsilon value 
%                  and the second option is 1 for float and 0 for binary, 
%                  prec is the precision of the variables defaults [1e-6 1]

if nargin<3
  options=1;
end

nargin;
%pause

numVars     = size(bounds,1);                %Number of variables
rng         = (bounds(:,2)-bounds(:,1))';    %The variable ranges'
%pause

if options==1 %Float GA
    %pop = zeros(num,numVars); 	%Allocate the new population
    pop=(ones(num,1)*rng).*(rand(num,numVars))+(ones(num,1)*bounds(:,1)');
    %pause
else %Binary GA
    bits=calcbits(bounds,options(1))
    pop = round(rand(num,sum(bits)));
    %pause
end
