function [tracklist] = savepos2tracklist(NT,savePos)
% This function converts from coordinates in terms of edge index and
% position along edge to real-space coordianates
% tracklist{pc}(tc,:) contains the real-space position of particle pc on
% step tc

% savePos should be 3d matrix, with dimensions: particles, (info), time index
% savePos(pc,1,tc) = edge particle pc is on at time tc
% savePos(pc,2,tc) = how far along the edge it is on

% network NT should have cumedgelen set

[nPart, ~, nCount] = size(savePos);
tracklist = {};

for np = 1 : nPart                                  % np is particle count
    tracklist{np} = zeros(nCount, NT.dim);
    for nc = 1 : nCount                             % nc is time count
        ec = savePos(np,1,nc); % what edge are we on
        lvc = NT.cumedgelen{ec};      % lvc is array of cumulative length
        
        tracklist{np}(nc, :) = interp1(lvc,NT.edgepath{ec},savePos(np,2,nc));               
    end
end

end

