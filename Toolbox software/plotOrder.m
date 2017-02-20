function plotRows= plotOrder(mapData)
prompt = {'Total rows to plot '};
dlg_title = 'Plot order ';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
plotRows=str2num(answer{1});
