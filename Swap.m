function [AA,AB] = Swap(AA,AB) %% Taken from Mathworks
error(nargchk(2,2,nargin)) ;
if ischar(AA) && ischar(AB),
    % command syntax: SWAP VAR1 VAR2
    evalstr = sprintf('[%s,%s] = deal(%s,%s) ;',AB,AA,AA,AB) ;
elseif ~isempty(inputname(1)) && ~isempty(inputname(2)),    
    % function syntax: SWAP(VAR1,VAR2)
    evalstr = sprintf('[%s,%s] = deal(%s,%s) ;',inputname(2),inputname(1),inputname(1),inputname(2)) ;
else
    % No valid command syntax, force error 
    evalstr = '[:)' ;
end
evalin('caller',evalstr,'error(''Requires two (string names of) existing variables'')') ;   
end