function [dist,newTxPox] = wDistCalc(R,TxPos,RxPos,varargin)
    %Calculate wraparound distances and return wrapped TxPos
    %R - Cell radius
    %TxPos - Tx position
    %RxPos - Rx position
    %To disable wraparound mode, simply add any input argumments in
    %addition to the previous three.
    if isempty(varargin)
        seventxpos = cell(7,1);
        % 21-cell wrap around
        seventxpos{1} = TxPos;
        seventxpos{2} = TxPos + directedCalc(R,'UR',6,'R',3);
        seventxpos{3} = TxPos + directedCalc(R,'UL',6,'UR',3);
        seventxpos{4} = TxPos + directedCalc(R,'L',6,'UL',3);
        seventxpos{5} = TxPos + directedCalc(R,'DL',6,'L',3);
        seventxpos{6} = TxPos + directedCalc(R,'DR',6,'DL',3);
        seventxpos{7} = TxPos + directedCalc(R,'R',6,'DR',3);

        sevendist = cellfun(@(x) norm(x-RxPos), seventxpos);

        [dist,I] = min(sevendist(:));
        newTxPox = seventxpos{I};
    else
        newTxPox = TxPos;
        dist = norm(TxPos-RxPos);
    end
    
    
end

function dplm = directedCalc(R,varargin)
    %R is cell radius
    nvar = length(varargin);
    dplm = [0,0,0];   %Displacement
    for ii = 1:2:(nvar-1)
        direct = varargin{ii};
        numstep = varargin{ii+1};
        vlen = R*numstep;
        switch direct
            case 'R'    %Right
                dplm = dplm + vlen*[1,0,0];
            case 'L'    %Left
                dplm = dplm + vlen*[-1,0,0];
            case 'UR'   %Up Right
                dplm = dplm + vlen*[cosd(60),sind(60),0];
            case 'DR'   %Down Right
                dplm = dplm + vlen*[cosd(-60),sind(-60),0];
            case 'UL'   %Up Left
                dplm = dplm + vlen*[cosd(120),sind(120),0];
            case 'DL'   %Down Left
                dplm = dplm + vlen*[cosd(-120),sind(-120),0];
        end
    end
    
end
