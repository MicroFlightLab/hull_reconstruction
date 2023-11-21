function [results] = fetchNext_waitbar(parclean,frms,mov,codepart,varargin)
% used for parallel computing, fetch next to get results. Also shows
% waitbar. to sapress waitbar define 'waitbar' as 0

parser = inputParser;
addParameter(parser,'waitbar',1); 
parse(parser, varargin{:})

h = waitbar(0,sprintf('mov%d : Waiting for FevalFutures to complete...',mov),'Name',codepart);
N = length(parclean);
results = cell(1,N);
for idx = 1:N
    
    [completedIdx,thisResult] = fetchNext(parclean);
    % Store the result
    results{completedIdx} = thisResult;
    % Update waitbar
    waitbar(idx/N,h,sprintf('mov%d : frame: %d',mov,frms(completedIdx)));
end
delete(h);

end

