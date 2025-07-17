# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import xarray as xr

var_name = 'density_bulk'
run_name = 'KAN_M_100_layers_0'
site ='KAN_M'
zero_surf=False
ylim=10
filename = "C:/Data_save/Data JoG 2020/DYE-2_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10/snic_bin_1.nc"

print('plotting',var_name, 'from',run_name)
ds = xr.open_dataset(filename).transpose()
ds = ds.resample(time='6H').nearest()

if not zero_surf:
    ds['surface_height'] = -ds.depth.isel(level=-1) + ds.depth.isel(level=-1).isel(time=0)
    ds['depth'] = ds.depth + ds.surface_height

if var_name == "slwc":
    # change unit to mm / m3
    ds[var_name] = ds[var_name] * 1000 / ds.depth

# default plot infos
label = var_name
cmap = 'magma'
vmin = np.percentile(ds[var_name], 5)
vmax = np.percentile(ds[var_name], 95)
        
# updating for pre-set values
plot_info = pd.read_csv('lib/plot_info.csv', skipinitialspace=True)
if var_name in plot_info.variable_name.to_list():
    plot_info = plot_info.set_index('variable_name')
    label = plot_info.loc[var_name, 'label']
    cmap = plot_info.loc[var_name, 'cmap']
    if ~np.isnan(plot_info.loc[var_name].vmin):
        vmin = plot_info.loc[var_name, 'vmin']
        vmax = plot_info.loc[var_name, 'vmax']
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
plt.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.1, hspace=0.2)
fig.suptitle(site)

im = ax.pcolormesh(
               ds.time.expand_dims(dim={"level": ds.level.shape[0]+1}).transpose(), 
               np.hstack([ds.surface_height.values.reshape([-1,1]),
                              ds.depth.values]),
               ds[var_name].isel(time=slice(1,None)),
               shading='flat',
               cmap = cmap, vmin=vmin, vmax=vmax
               )

if not zero_surf:
    ax.plot(ds.time, ds.surface_height, linewidth=2, color="k")
    
plt.colorbar(im, label=label, ax=ax)
ax.invert_yaxis()
if ylim:
    if len(ylim)==1: ax.set_ylim(ylim, ax.get_ylim()[1])
    if len(ylim)==2: ax.set_ylim(np.max(ylim), np.min(ylim))
# %%
ax.set_ylabel("Depth (m)")col = lines(4);
year = 1977
[~, ind_1977] =  min(abs(data_surf{1}.time - datenum(1977,7,1)));
[~, ind_2010] =  min(abs(data_surf{1}.time - datenum(2010,7,1)));
axis tight


f = figure('Visible',vis,'OuterPosition',[0 0 20 13]);
ha = tight_subplot(1,2,0.1,[0.15 0.1],0.1);
% ha(2).Visible = 'off';
% ha(6).Visible = 'off';
% ha(1).Position(3) = 0.8;
% ha(5).Position(3) = 0.8;

% set(f,'CurrentAxes',ha(1))
% hold on
% for ii =1:1
% %     x = data_surf{ii}.time;
% %     y1 = depth_annual_layers{ii}(find(ind_mod_in_77,1,'first'),:)';
% %     y2 = depth_annual_layers{ii}(find(ind_mod_in_77,1,'last'),:)';
% %     x2 = [x, fliplr(x)];
% %     ind1 = min(find(~isnan(y1),1,'first'), find(~isnan(y2),1,'first'))+10;
% %     ind2 = min(find(~isnan(y1),1,'last'), find(~isnan(y2),1,'last'));
% %     y1=y1(ind1:ind2);
% %     y2=y2(ind1:ind2);
% %     y1(isnan(y1))=0;
% %     y2(isnan(y2))=0;
% %     inBetween = [y1, fliplr(y2)];
% %     x2=x2(ind1:ind2,:);
% %     fill(x2, inBetween, RGB('light light green'));
% 
% plot(data_surf{ii}.time,...
%         depth_annual_layers{ii}([find(ind_mod_in_77,1,'first') ...
%             find(ind_mod_in_77,1,'last')],:),...
%             'Color',RGB('light green'));
%         
% plot(data_surf{ii}.time,...
%         depth_annual_layers{ii}([find(ind_mod_in_10,1,'first') ...
%             find(ind_mod_in_10,1,'last')],:),...
%             'Color',RGB('light blue'));
% end
% plot(data_surf{ii}.time(ind_1977)*ones(size(CC77_1.depth(ind_77_in_mod))),...
%     CC77_1.depth(ind_77_in_mod),'ok','MarkerFaceColor','k')
% plot(data_surf{ii}.time(ind_2010)*ones(size(CC10.depth(ind_10_in_mod))),...
%     CC10.depth(ind_10_in_mod),'ok','MarkerFaceColor','k')
% 
% set(gca,'YDir','reverse','XAxisLocation','top')
% xlabel('Year'); ylabel('Depth (m)');
% set_monthly_tick(data_surf{ii}.time);
% xlim(data_surf{ii}.time([1 end]))
    
set(f,'CurrentAxes',ha(1))
hold on
h = [];
for ii =4:-1:1
    h(ii) = plot(CC77_1.depth(ind_77_in_mod),...
        flipud(depth_annual_layers{ii}(ind_mod_in_77,ind_1977)),'o',...
        'Color',col(ii,:),'MarkerFaceColor',col(ii,:));
    h(ii) = plot(CC77_2.depth(ind_77_in_mod),...
        flipud(depth_annual_layers{ii}(ind_mod_in_77,ind_1977)),'x',...
        'Color',col(ii,:),'MarkerFaceColor',col(ii,:));
end
    b(1) = plot(NaN,NaN,'xk', 'MarkerFaceColor','k');
    b(2) = plot(NaN,NaN,'ok', 'MarkerFaceColor','k');
    text(CC77_1.depth(ind_77_in_mod),...
        flipud(depth_annual_layers{ii}(ind_mod_in_77,ind_1977))+2,...
        num2str(CC77_1.year(ind_77_in_mod)),'Rotation',90);
    plot([0 40],[0 40],'k')
    axis tight; box on; grid on
    xlim([3 10])
    ylim([3 14])
    xlabel('Observed depth (m)')
    ylabel('Simulated depth (m)')
title('Annual layers in the 1977 cores','FontSize',13)
% legend(h,model_list,'location','eastoutside','interpreter','none')
h_now = gca;
h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
    h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.08, 'A');
h_text.FontWeight='bold';
h_text.FontSize=16;

set(f,'CurrentAxes',ha(2))
hold on
h = [];
for ii =4:-1:1
    ind1 = find(ismember(CC10.year, years_all));
    ind2 = find(ismember(years_all, CC10.year));
    h(ii) = plot(CC10.depth(ind1),...
        flipud(depth_annual_layers{ii}(ind2,ind_2010)),'o',...
        'Color',col(ii,:),'MarkerFaceColor',col(ii,:));
end
    text(CC10.depth(ind1(1:2:end)),...
        flipud(depth_annual_layers{ii}(ind2((1:2:end)),ind_2010))+2,...
        num2str(CC10.year(ind1((1:2:end)))),'Rotation',90);
    plot([0 40],[0 40],'k')
    axis tight; box on; grid on
    xlim([2 30])
    ylim([2 45])
title('Annual layers in the 2010 core','FontSize',13)
% legend(h,model_list2,'location','eastoutside','interpreter','none')
    xlabel('Observed depth (m)')
    ylabel('Simulated depth (m)')
%     ha(4).YAxisLocation = 'right';
%     ha(3).XMinorTick = 'on';
%     ha(3).XMinorTick = 'on';

    h_now = gca;
h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
    h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.08, 'B');
h_text.FontWeight='bold';
h_text.FontSize=16;
print(f, sprintf('%s/9-tracking annual layers',OutputFolder), '-djpeg') 

%% 10- Tracking debris horizon and water percolation
    f = figure('Visible',vis,'OuterPosition',[0 0 28 13]);
h = [];
g=[];
ha= tight_subplot(1,2,0.08,[0.3 0.18] ,0.12 );                      
% ha(2).Position(4) = 0.4;
col = lines(4);
set(f,'CurrentAxes',ha(1))
hold on
for ii =2:4
    h(ii-1) = plot(data_surf{ii}.time,data_surf{ii}.depth_17,...
        'Color',col(ii,:),'LineWidth',2);
%        plot(data_surf{ii}.time,data_surf{ii}.depth_66,...
%         'Color',col(ii,:),'LineWidth',2);
%        plot(data_surf{ii}.time,data_surf{ii}.depth_12,...
%         'Color',col(ii,:),'LineWidth',2);
       
    depth_1 = data_subsurf{ii}.depth_act;
    depth_1(data_subsurf{ii}.slwc<1e-6) = 0;
    perc_depth = max(depth_1,[],1);
    
%     g(ii-1) = plot(data_surf{ii}.time,...
%         perc_depth,'-',...
%         'Color',col(ii,:),'LineWidth',2);
    set(f,'CurrentAxes',ha(2))
    hold on
    g(ii-1) = plot(data_surf{ii}.time,...
        perc_depth,'-',...
        'Color',col(ii,:),'LineWidth',1);
    disp(model_list{ii})
    disp(max(perc_depth))
    disp(data_surf{ii}.depth_17(end))
    set(f,'CurrentAxes',ha(1))
end
[~, ind_2017] =  min(abs(data_surf{1}.time - datenum(2017,7,1)));
axis tight

% h(5) = errorbar([1 1].*data_surf{1}.time(ind_2017), [6 33.5], [0.1, 2], ...
%     'ok','MarkerFaceColor','k');
h(4) = plot([1 1].*data_surf{1}.time(ind_2017), [32.5 32.5],...
    'ok','MarkerFaceColor','k');
set_monthly_tick(data_surf{4}.time)
% xlim(data_surf{2}.time([ind_2017 end]))
ylim([0 80])
set(gca,'YDir','reverse','XTickLabelRotation',45)
xlabel('Year')
grid on
ylabel('Depth (m)')
d1 = plot(NaN,NaN,'Color',RGB('gray'),'LineWidth',2);
d2 = plot(NaN,NaN,'--','Color',RGB('gray'),'LineWidth',2);
legendflex(h,{model_list2{2:4}},...
    'ref',f,'anchor',{'n' 'n'},...
    'ncol',3,...
    'Location','northoutside','Interpreter','none')
    h_now = gca;
h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
    h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.08, 'A');
h_text.FontWeight='bold';
h_text.FontSize=16;

title('Top of the debris layer','FontSize',14)
set(f,'CurrentAxes',ha(2))
title({'Meltwater percolation depth'},'FontSize',14)

hold on
[~, ind_2085] =  min(abs(data_surf{2}.time - datenum(2085,1,1)));
set_monthly_tick(data_surf{4}.time(ind_2085:end))
xlim(data_surf{2}.time([ind_2085 end]))
% ylim([0 80])
set(gca,'YDir','reverse','XTickLabelRotation',45,'YAxisLocation','right')
xlabel('Year')
grid on
ylabel('Depth (m)')

% Create rectangle
h_ann = annotation(f,'rectangle',...
    [ha(1).Position(1)  + ha(1).Position(3)*0.8...
    ha(1).Position(2)+ha(1).Position(4)*0.95 0.07 0.03]);
annotation(f,'line',[h_ann.Position(1)  + h_ann.Position(3)...
    ha(2).Position(1)],...
    [h_ann.Position(2) ha(2).Position(2)]);
annotation(f,'line',[h_ann.Position(1)  + h_ann.Position(3)...
    ha(2).Position(1)],...
    [h_ann.Position(2)+ h_ann.Position(4) ...
    ha(2).Position(2)+ha(2).Position(4)]);
    h_now = gca;
h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
    h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.08, 'B');
h_text.FontWeight='bold';
h_text.FontSize=16;
print(f, sprintf('%s/_10-tracking_debris_horizon',OutputFolder), '-djpeg') 
