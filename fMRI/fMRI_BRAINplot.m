clear all; clc; close all;

addpath(genpath('../support/brainPlotUtils'));


flag_type = 'RS';
%flag_type = 'TS';

%flag_dir = 'fw';
flag_dir = 'bw';

dataFilename = ['results/', flag_type, '_', flag_dir, '_list.mat'];
load(dataFilename);

[returnCode, hostName] = system('hostname');

% if running on NAS
if(strcmp(deblank(hostName),'miplabsrv3'))
   CodeBook=load('Glasser360_2mm_codebook');
   CodeBook=CodeBook.codeBook;

   CM = zeros(360,360);
   T_conn=0.5;
   Factor_Col=1;
   Exp_Sphere=1;
   View=1;

   Colormap_edges='jet';
   Colormap_nodes='winter';
   Gamma=0.5;
   LinearWeight=0.7;
   CA=[-1 1];

   CM2=zeros(size(CM));
   
   if(flag_dir == 'fw')
      UM = zeros(1,360); UM(ROI_fw_list) = 100;
   else
      UM = zeros(1,360); UM(ROI_bw_list) = 100;
   end
   
   filename = ['images/', flag_type, '_brain_ROI_', flag_dir, '.fig'];
   
   PlotBrainGraph(CM2,1*UM,0.8*UM,CodeBook,T_conn,max(UM),...
         Factor_Col,Exp_Sphere,View,Colormap_nodes,Colormap_edges,Gamma,...
         LinearWeight,CA, filename)
   exit
else
   % if running on laptop
   fig = openfig(['images/', flag_type, '_brain_ROI_', flag_dir, '.fig']);
   print(fig, '-depsc2', ['images/', flag_type, '_brain_ROI_', flag_dir]);
end


