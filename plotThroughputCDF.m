function plotThroughputCDF(YUO,varargin)
    % Plot CDF of throughputs
    CellIntIdx = YUO.cellOfInterestIdx;
    DLThpt = YUO.Throughput(1,1:(end-2),CellIntIdx,:);
    ULThpt = YUO.Throughput(2,1:(end-2),CellIntIdx,:);
    if isempty(varargin)
        plotCDF(DLThpt(:),'DL');
        plotCDF(ULThpt(:),'UL');
    else
        switch varargin{1}
            case 'DL'
                plotCDF(DLThpt(:),'DL');
            case 'UL'
                plotCDF(ULThpt(:),'UL');
        end
    end
end

function plotCDF(x,msg)
    figure
    ecdf(x)
    title([msg,' Throughputs'])
end