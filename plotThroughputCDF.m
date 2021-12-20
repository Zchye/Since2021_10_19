function plotThroughputCDF(YUO,varargin)
    % Plot CDF of throughputs
    if isempty(varargin)
        DLThpt = YUO.Throughput(1,:,:,:);
        ULThpt = YUO.Throughput(2,:,:,:);
        plotCDF(DLThpt(:),'DL');
        plotCDF(ULThpt(:),'UL');
    else
        switch varargin{1}
            case 'DL'
                DLThpt = YUO.Throughput(1,:,:,:);
                plotCDF(DLThpt(:),'DL');
            case 'UL'
                ULThpt = YUO.Throughput(2,:,:,:);
                plotCDF(ULThpt(:),'UL');
        end
    end
end

function plotCDF(x,msg)
    figure
    ecdf(x)
    title([msg,' Throughputs'])
end