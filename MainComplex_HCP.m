%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

ROOTDIR = 'J:\HCP';
MIPDIR = fullfile( ROOTDIR, 'AtlasBasedAnalsys\MI_MIP' );
PROJECTDIR = fullfile( ROOTDIR, 'MainComplex_HCP' );
ATLASDIR = fullfile( PROJECTDIR, 'atlas');

addpath(genpath( fullfile( PROJECTDIR, 'Toolbox' )));

% load suject ID
cd( PROJECTDIR );
load( 'id_common.mat' ); % id_common: 957 x 1, id_common_half1: 479 x 1, id_common_half2: 478 x 1

sessions = {'task_WM','task_Gambling','task_Motor','task_Language','task_Social','task_Relational','task_EM','rest'};
session_names = {'Working memory','Gambling','Motor','Language','Social Cognition','Relational Processing','Emotion Processing','Rest'};

Atlas = {'Yeo_7','Yeo_17'};
N_Networks = [7 17];

ITR = 10000;

Yeo_2011_17networks_names_abr = {'VN-c','VN-p','SMN-a','SMN-b','DAN-a','DAN-b','VAN','SN','LN-a','LN-b','ECN-c','ECN-a','ECN-b','TPN','DMN-c','DMN-a','DMN-b'};
Yeo_2011_7networks_names_abr = {'VN','SMN','DAN','SN/VAN','LN','ECN','DMN'};

%%
% for atlas = 1:length(Atlas)
%     for ii = 1:2
%         if ii ==1, tmp_id = id_common_half1; else, tmp_id = id_common_half2;end
%         for sess = 1:8
%             cd( fullfile( MIPDIR, Atlas{atlas}, sessions{sess} ));
%             datalist = dir('*.mat');
%
%             % initialazation
%             Imip = zeros([length(tmp_id) N_Networks(atlas)]);
%             Complex = zeros([ length(tmp_id) N_Networks(atlas)]);
%
%             m = 0;
%             for nn = 1:length(datalist)
%                 load( datalist(nn).name ); % CMP
%                 [~,id_str,~]=fileparts(datalist(nn).name);
%                 if sum(tmp_id == str2double(id_str))
%                     fprintf('%s\n', datalist(nn).name);
%                     m = m + 1;
%                     mi_mip = zeros(N_Networks(atlas));
%                     if length(CMP.complexes)==1
%                         mi_mip(CMP.complexes{1,1},CMP.complexes{1,1}) = CMP.phis_complexes(1);
%                     else
%                         for mm = length(CMP.complexes):-1:1
%                             mi_mip(CMP.complexes{mm,1},CMP.complexes{mm,1}) = CMP.phis_complexes(mm);
%                         end
%                     end
%                     Imip(m,:) = max(mi_mip);
%                     Complex(m,CMP.complexes{1}) = 1;
%
%                 end
%
%             end
%             cd( PROJECTDIR ); cd('data');
%             save( fullfile( strcat('Complex_', Atlas{atlas},'_',sessions{sess}, '_',sprintf('%0.3d',N_Networks(atlas)),'_half',int2str(ii),'.mat') ),'Imip','Complex');
%
%         end
%     end
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Permutation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for atlas = 1:length(Atlas)
    for ii = 1:2
        for sess = 1:8
            cd( PROJECTDIR ); cd('data');
            f_name = fullfile( strcat('Complex_', Atlas{atlas},'_',sessions{sess}, '_',sprintf('%0.3d',N_Networks(atlas)),'_half',int2str(ii),'.mat') ) ;
            load( f_name ); % Imip,Complex: N_sub x N_networks
            fprintf( 'Permutation test: Group %d, session %d, %s\n', ii, sess, f_name );
            
            [mu_Imip,p_Imip]       = permutation_test(Imip,ITR);
            [mu_Complex,p_Complex] = permutation_test(Complex,ITR);
            
            cd( PROJECTDIR ); cd('permutatoin_test');
            save( fullfile( strcat('Complex_', Atlas{atlas},'_',sessions{sess}, '_',sprintf('%0.3d',N_Networks(atlas)),'_half',int2str(ii),'.mat') ),'mu_Imip','p_Imip','mu_Complex','p_Complex');
            
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  make spider plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for atlas = 1:length(Atlas)
    if atlas == 1
        Order_for_figure = 1:7;
        label_names = Yeo_2011_7networks_names_abr;
        AxesLimits = [ zeros(1,N_Networks(atlas)); ones(1,N_Networks(atlas))*1.0];
    elseif atlas ==2
        Order_for_figure = 1:17;
        label_names = Yeo_2011_17networks_names_abr;
        AxesLimits = [ zeros(1,N_Networks(atlas)); ones(1,N_Networks(atlas))*1.2];
    end
    ParticipationRatio1 = zeros([N_Networks(atlas) length(sessions)]);
    ParticipationRatio2 = zeros([N_Networks(atlas) length(sessions)]);
    MutualInformation1 = zeros([N_Networks(atlas) length(sessions)]);
    MutualInformation2 = zeros([N_Networks(atlas) length(sessions)]);
    for sess = 1:8
        cd( PROJECTDIR ); cd('data');
        f_name = fullfile( strcat('Complex_', Atlas{atlas},'_',sessions{sess}, '_',sprintf('%0.3d',N_Networks(atlas)),'_half1.mat') ) ;
        load( f_name ); % Imip,Complex: N_sub x N_networks
        ParticipationRatio1(:,sess) = mean(Complex)';
        MutualInformation1(:,sess) = mean(Imip)';
        f_name = fullfile( strcat('Complex_', Atlas{atlas},'_',sessions{sess}, '_',sprintf('%0.3d',N_Networks(atlas)),'_half2.mat') ) ;
        load( f_name ); % Imip,Complex: N_sub x N_networks
        ParticipationRatio2(:,sess) = mean(Complex)';
        MutualInformation2(:,sess) = mean(Imip)';
    end
    
    figsize = [200 200 400 250];
    option.AxesLimits=AxesLimits;
    option.AxesLabels=label_names;
    option.AxesInterval=3;
    option.AxesPrecision=2;
    option.AxesDisplay='one';
    option.AxesColor=[0.9 0.9 0.9];
    option.AxesOffset=0;
    option.AxesLabelsEdge='none';
    option.AxesFontSize=6;
    option.LineWidth=0.75;
    option.MarkerSize=5;
    option.LabelFontSize=8;
    option.Direction='clockwise';
    
    figure('Position',figsize);
    subplot(1,2,1);
    make_spkder_plot(MutualInformation1(Order_for_figure,:)',option)
    legend(session_names, 'Location','best'); legend('boxoff');
    subplot(1,2,2);
    make_spkder_plot(MutualInformation2(Order_for_figure,:)',option)
    
    figure('Position',figsize);
    option.AxesLimits=[ zeros(1,N_Networks(atlas)); ones(1,N_Networks(atlas))*1];
    
    subplot(1,2,1);
    make_spkder_plot(ParticipationRatio1(Order_for_figure,:)',option)
    legend(session_names, 'Location','best'); legend('boxoff');
    
    subplot(1,2,2);
    make_spkder_plot(ParticipationRatio2(Order_for_figure,:)',option)
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for making brain map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% each session
for atlas = 1:length(Atlas)
    for ii =1:2
        
        cd( fullfile( ATLASDIR, Atlas{atlas} ));
        atlasfilename = dir('*dscalar.nii');
        
        atlas_nii = ft_read_cifti(atlasfilename(1).name);
        roi_atlas = atlas_nii.dscalar;
        no_rois = max(roi_atlas);
        
        for sess = 1:8
            
            cd( PROJECTDIR ); cd('permutatoin_test');
            load( fullfile( strcat('Complex_', Atlas{atlas},'_',sessions{sess}, '_',sprintf('%0.3d',N_Networks(atlas)),'_half',int2str(ii),'.mat') ));
            
            p_Imip(p_Imip==0)=0.00001;
            pt = FDR(p_Imip,0.05);
            Imip_sig = p_Imip<=pt;
            Imip_voxel = zeros( size(roi_atlas) );
            
            p_Complex(p_Complex==0)=0.00001;
            pt = FDR(p_Complex,0.05);
            Complex_sig = p_Complex<=pt;
            Complex_voxel = zeros( size(roi_atlas) );
            
            for jj = 1:N_Networks(atlas)
                Ivoxel = roi_atlas==jj ;
                if Imip_sig( jj )
                    Imip_voxel(Ivoxel) = mu_Imip( jj );
                end
                if Complex_sig( jj )
                    Complex_voxel(Ivoxel) = mu_Complex( jj );
                end
            end
            
            cd( PROJECTDIR ); cd('brain_map');
            
            Imip_nii = atlas_nii;
            Imip_nii.x1 = Imip_voxel;
            ft_write_cifti( strcat('Imip_', Atlas{atlas},'_',sessions{sess},int2str(ii)),Imip_nii,'parameter','x1');
            
            Complex_nii = atlas_nii;
            Complex_nii.x1 = Complex_voxel;
            ft_write_cifti( strcat('Complex_', Atlas{atlas},'_',sessions{sess},int2str(ii)),Complex_nii,'parameter','x1');
        end
    end
end
%%
% across conditions
for atlas = 1:length(Atlas)
    for ii =1:2
        
        cd( fullfile( ATLASDIR, Atlas{atlas} ));
        atlasfilename = dir('*dscalar.nii');
        
        atlas_nii = ft_read_cifti(atlasfilename(1).name);
        roi_atlas = atlas_nii.dscalar;
        no_rois = max(roi_atlas);
        
        p_Imip_all = zeros([N_Networks(atlas) 1]);
        p_Complex_all = zeros([N_Networks(atlas) 1]);
        
        for sess = 1:8
            cd( PROJECTDIR ); cd('permutatoin_test');
            load( fullfile( strcat('Complex_', Atlas{atlas},'_',sessions{sess}, '_',sprintf('%0.3d',N_Networks(atlas)),'_half',int2str(ii),'.mat') ));
            p_Imip(p_Imip==0)=0.00001;
            pt = FDR(p_Imip,0.05);
            p_Imip_all = p_Imip_all + ( p_Imip <= pt );
            p_Complex(p_Complex==0)=0.00001;
            pt = FDR(p_Complex,0.05);
            p_Complex_all = p_Complex_all + ( p_Complex <= pt);
        end
        
        Imip_voxel = zeros( size(roi_atlas) );
        Complex_voxel = zeros( size(roi_atlas) );
        for jj = 1:N_Networks(atlas)
            Ivoxel = roi_atlas==jj ;
            Imip_voxel(Ivoxel) = p_Imip_all( jj );
            Complex_voxel(Ivoxel) = p_Complex_all( jj );
        end
        
        cd( PROJECTDIR ); cd('brain_map');
        
        Imip_nii = atlas_nii;
        Imip_nii.x1 = Imip_voxel;
        ft_write_cifti( strcat('Imip_', Atlas{atlas},'_all_',int2str(ii)),Imip_nii,'parameter','x1');
        
        Complex_nii = atlas_nii;
        Complex_nii.x1 = Complex_voxel;
        ft_write_cifti( strcat('Complex_', Atlas{atlas},'_all_',int2str(ii)),Complex_nii,'parameter','x1');
    end
end
% dscalar.nii files are visualized in Workbench.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ICC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for atlas = 1:length(Atlas)
    ParticipationRatio1 = zeros([N_Networks(atlas) length(sessions)]);
    ParticipationRatio2 = zeros([N_Networks(atlas) length(sessions)]);
    MutualInformation1 = zeros([N_Networks(atlas) length(sessions)]);
    MutualInformation2 = zeros([N_Networks(atlas) length(sessions)]);
    for sess = 1:8
        cd( PROJECTDIR ); cd('data');
        f_name = fullfile( strcat('Complex_', Atlas{atlas},'_',sessions{sess}, '_',sprintf('%0.3d',N_Networks(atlas)),'_half1.mat') ) ;
        load( f_name ); % Imip,Complex: N_sub x N_networks
        ParticipationRatio1(:,sess) = mean(Complex)';
        MutualInformation1(:,sess) = mean(Imip)';
        f_name = fullfile( strcat('Complex_', Atlas{atlas},'_',sessions{sess}, '_',sprintf('%0.3d',N_Networks(atlas)),'_half2.mat') ) ;
        load( f_name ); % Imip,Complex: N_sub x N_networks
        ParticipationRatio2(:,sess) = mean(Complex)';
        MutualInformation2(:,sess) = mean(Imip)';
    end
    
    cd( PROJECTDIR ); cd('icc');
    
    [icc.MI.r, icc.MI.LB, icc.MI.UB, icc.MI.F, icc.MI.df1, icc.MI.df2, icc.MI.p] = ICC(MutualInformation1, 'C-k') ;
    [icc.PR.r, icc.PR.LB, icc.PR.UB, icc.PR.F, icc.PR.df1, icc.PR.df2, icc.PR.p] = ICC(ParticipationRatio1, 'C-k') ;
    save('ICC_group1.mat','icc');
    iccR_MI(1,atlas) = icc.MI.r;    iccR_PR(1,atlas) = icc.PR.r;
    iccU_MI(1,atlas) = icc.MI.UB;   iccU_PR(1,atlas) = icc.PR.UB;
    
    [icc.MI.r, icc.MI.LB, icc.MI.UB, icc.MI.F, icc.MI.df1, icc.MI.df2, icc.MI.p] = ICC(MutualInformation2, 'C-k');
    [icc.PR.r, icc.PR.LB, icc.PR.UB, icc.PR.F, icc.PR.df1, icc.PR.df2, icc.PR.p] = ICC(ParticipationRatio2, 'C-k') ;
    save('ICC_group2.mat','icc');
    iccR_MI(2,atlas) = icc.MI.r;    iccR_PR(2,atlas) = icc.PR.r;
    iccU_MI(2,atlas) = icc.MI.UB;   iccU_PR(2,atlas) = icc.PR.UB;
    
end
figure('position',[200 200 300 200]);
for ii = 1:2
    subplot(1,2,ii);
    if ii == 1, tmp_data = iccR_MI; else, tmp_data = iccR_PR; end
    if ii == 1, tmp_error = iccU_MI-iccR_MI; else, tmp_error = iccU_PR-iccR_PR; end
    hold on;
    b = bar(tmp_data); 
    b(1).FaceColor = [.5 .5 .5];    
    b(2).FaceColor = [.7 .7 .7];

    [ngroups,nbars] = size(tmp_data);
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x',tmp_data,tmp_error,'k','linestyle','none')
    
    box off;
    xticklabels({'group1','group2'});
    xticks([]);
    yticks([0.0 0.2 0.4 0.6 0.8 1.0]);
    if ii == 1, xticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'}); else, yticklabels({}); end
    if ii == 1, title('IMIP'); else, title('Participation Rate'); end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Permutation test funtion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mu,p] = permutation_test(data,itr)
s = size( data );
data_m = mean( data );
data_m_perm = zeros( [ itr s(2) ]);
p = zeros([s(2) 1]);
for ii = 1:itr
    data_tmp = data;
    for jj = 1:s(1)
        data_tmp(jj,:) = data_tmp(jj, randperm(s(2)));
    end
    data_m_perm(ii,:) = mean( data_tmp );
end
for ii = 1:s(2)
    p(ii) = mean( data_m_perm(:,ii) > data_m(ii) );
end
mu = data_m(:);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  make_spkder_plot function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_spkder_plot(data,opt)
option = opt;
disp(option)
spider_plot(data,'AxesLimits',option.AxesLimits,...
    'AxesLabels',option.AxesLabels,...
    'AxesInterval',option.AxesInterval,...
    'AxesPrecision',option.AxesPrecision,...
    'AxesDisplay',option.AxesDisplay,...
    'AxesColor',option.AxesColor,...
    'AxesOffset',option.AxesOffset,...
    'AxesLabelsEdge',option.AxesLabelsEdge,...
    'AxesFontSize',option.AxesFontSize,...
    'LineWidth',option.LineWidth,...
    'MarkerSize',option.MarkerSize,...
    'LabelFontSize',option.LabelFontSize,...
    'Direction',option.Direction);
end