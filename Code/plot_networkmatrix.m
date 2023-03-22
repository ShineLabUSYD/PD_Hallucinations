%% Linking Axes
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\schaef_id.mat');
% 3-Dimensional plots
hexMap = [0 0.4470 0.7410 0.8500 0.3250 0.0980 0.9290 0.6940 0.1250 0.4940 0.1840 0.5560 0.4660 0.6740 0.1880 0.3010 0.7450 0.9330 0.6350 0.0780 0.1840];
myColorMap = zeros(length(hexMap)/3, 3); % Preallocate
a = 1;
for k = 1 : length(myColorMap)
	r = hexMap(a);
	g = hexMap(a+1);
	b = hexMap(a+2);
	myColorMap(k, :) = [r, g, b];
    a = a + 3;
end

[~,order] = sort(schaef_id);
schaef_id_sort = schaef_id(order);
matrix_order = unique_diff_conpd(order,order);


f1=figure('color','w','position',[500 200 500 300]);
ax1 = axes('parent',f1,'position',[0.2 0.1 0.7 0.8],'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[]);
imagesc(matrix_order);
colormap(ax1,bluewhitered()); % Include colormap as input for function
cbar = colorbar;
cbar.Position = [0.91 0.3005 0.021 0.4152];
hold(ax1,'all')
set(gca, 'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])
%axis square

hold on
ax2 = axes('parent',f1,'position',[0.185 0.1 0.015 0.8],'xtick',[],'ytick',[]);
imagesc(schaef_id_sort);
colormap(ax2,myColorMap);
hold(ax2,'all')
set(gca, 'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])

ax3 = axes('parent',f1,'position',[0.2 0.9 0.7 0.02],'xtick',[],'ytick',[]);
imagesc(schaef_id_sort'); 
colormap(ax3,myColorMap);
set(gca, 'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])

