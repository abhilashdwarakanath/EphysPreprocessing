function samples = readDecimate(params,chan,startSamp,endSamp)

filename = ['NSX' num2str(chan) '.NC5'];
    samples = read_NC5(filename,startSamp,endSamp);
    for decs = 1:length(params.decimFacs)
        sprintf('Decimating channel %d, part %d of %d',chan,decs,numel(params.decimFacs))
        samples = decimate(samples,params.decimFacs(decs));
    end
end