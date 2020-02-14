start = 1;
while start ~= 0
    idx = 1:size(data,2);
    idx(1:start-1) = [];
    f = figure;
    f.WindowState = 'maximized';
    imagesc('CData',data(:,start:end),'XData',start:size(data,2),'YData',1:size(data,1));
    xlim([1 size(data,2)]);
    colorbar;
    colormap jet;
    prompt = {'First valid sample:'};
    dlgtitle = 'Input';
    dims = [1 50];
    definput = {'1'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    start = str2double(answer{1});
    close;
end
data = data(:,idx);
POD.idx = idx;
