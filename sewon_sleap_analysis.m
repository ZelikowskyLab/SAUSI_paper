% 2023 Dec 23 Sewon Park
% Find mean distance of each animal from SLEAP data
% change list below
% mother_root/Mice_list/File/animal body points

clc;
clear all;
close all;

%% Input

addpath('Z:\Jordan\elife')

% %Input
mother_root = 'Z:\Jordan\sausi_3F_rd3\sleap_csv';
Mice_list = {'16'}; %{'16','18','23','30','36','39','40','57','61','63','65','77', ...
    %'easi_A','easi_B','easi_C','easi_D','easi_E','easi_F','easi_G', ...
    %'easi_H','easi_I','easi_J','easi_K','easi_L','sausi_A','sausi_B', ...
   % 'sausi_C','sausi_D','sausi_E','sausi_F','sausi_G','sausi_H','sausi_I', ...
    %'sausi_J','sausi_K','sausi_L'};
% %

% % Output -- change the output path followed by directory
save_root = [mother_root '\Matlab results'];

if ~exist([save_root], 'dir')
    mkdir([save_root]);
end

% %

% % make excel file

cd(mother_root);
%outputCSV = fopen('Results_mean_distance_diff.csv', 'w');
%header = 'nn, nb, nt, bb, bt, tt\n';
%fprintf(outputCSV, header);
% %

% % Batch analysis
for animal_iter=1:length(Mice_list)
   
    % % csv read - SLEAP
    animal_id=Mice_list{animal_iter};
    
    File=[mother_root '\' animal_id, '.csv'];
    [Data, animal_pre]=xlsread(File);
    
    for i=1:length(Data)
        if char(animal_pre(i+1,1)) == 'track_0'
            animal(i,1)=0;
        else
            animal(i,1)=1;
        end
    end

    % % animal body points

    NumberofFrame = Data(length(Data),1)+1; % start frame=0
    
    frame_1=[]; Iscore_1=[]; nose.x_1=[]; nose.y_1=[]; neck.x_1=[]; neck.y_1=[]; earL.x_1=[]; earL.y_1=[]; earR.x_1=[]; earR.y_1=[]; bodyL.x_1=[]; bodyL.y_1=[]; bodyC.x_1=[]; bodyC.y_1=[]; bodyR.x_1=[]; bodyR.y_1=[]; tail.x_1=[]; tail.y_1=[]; tailE.x_1=[]; tailE.y_1=[];
    frame_2=[]; Iscore_2=[]; nose.x_2=[]; nose.y_2=[]; neck.x_2=[]; neck.y_2=[]; earL.x_2=[]; earL.y_2=[]; earR.x_2=[]; earR.y_2=[]; bodyL.x_2=[]; bodyL.y_2=[]; bodyC.x_2=[]; bodyC.y_2=[]; bodyR.x_2=[]; bodyR.y_2=[]; tail.x_2=[]; tail.y_2=[]; tailE.x_2=[]; tailE.y_2=[];
    
    for i=1:length(Data)
        if animal(i,1) == 0 % animal #1 (female)
            frame_1(length(frame_1)+1,1)=Data(i,1)+1; % frame (B)
            Iscore_1(length(Iscore_1)+1,1)=Data(i,2);
            nose.x_1(length(nose.x_1)+1,1)=Data(i,3); % nose x (D)
            nose.y_1(length(nose.y_1)+1,1)=Data(i,4); % nose y (E)
%             earL.x_1(length(earL.x_1)+1,1)=Data(i,9);
%             earL.y_1(length(earL.y_1)+1,1)=Data(i,10);
%             earR.x_1(length(earR.x_1)+1,1)=Data(i,12);
%             earR.y_1(length(earR.y_1)+1,1)=Data(i,13);
%             bodyL.x_1(length(bodyL.x_1)+1,1)=Data(i,15); 
%             bodyL.y_1(length(bodyL.y_1)+1,1)=Data(i,16); 
            bodyC.x_1(length(bodyC.x_1)+1,1)=Data(i,12); % body C x
            bodyC.y_1(length(bodyC.y_1)+1,1)=Data(i,13); % body C y
%             bodyR.x_1(length(bodyR.x_1)+1,1)=Data(i,21); 
%             bodyR.y_1(length(bodyR.y_1)+1,1)=Data(i,22);
            tail.x_1(length(tail.x_1)+1,1)=Data(i,21); % tail x
            tail.y_1(length(tail.y_1)+1,1)=Data(i,22); % tail y
%             tailE.x_1(length(tailE.x_1)+1,1)=Data(i,27); % tail x
%             tailE.y_1(length(tailE.y_1)+1,1)=Data(i,28); % tail y
        elseif animal(i,1) == 1 % animal #2 (male)
            frame_2(length(frame_2)+1,1)=Data(i,1)+1; % frame (B)
            Iscore_2(length(Iscore_2)+1,1)=Data(i,2);
            nose.x_2(length(nose.x_2)+1,1)=Data(i,3); % nose x (D)
            nose.y_2(length(nose.y_2)+1,1)=Data(i,4); % nose y (E)
%             earL.x_2(length(earL.x_2)+1,1)=Data(i,9);
%             earL.y_2(length(earL.y_2)+1,1)=Data(i,10);
%             earR.x_2(length(earR.x_2)+1,1)=Data(i,12);
%             earR.y_2(length(earR.y_2)+1,1)=Data(i,13);
%             bodyL.x_2(length(bodyL.x_2)+1,1)=Data(i,15); 
%             bodyL.y_2(length(bodyL.y_2)+1,1)=Data(i,16); 
            bodyC.x_2(length(bodyC.x_2)+1,1)=Data(i,12); % body C x
            bodyC.y_2(length(bodyC.y_2)+1,1)=Data(i,13); % body C y
%             bodyR.x_2(length(bodyR.x_2)+1,1)=Data(i,21); 
%             bodyR.y_2(length(bodyR.y_2)+1,1)=Data(i,22);
            tail.x_2(length(tail.x_2)+1,1)=Data(i,21); % tail x
            tail.y_2(length(tail.y_2)+1,1)=Data(i,22); % tail y
%             tailE.x_2(length(tailE.x_2)+1,1)=Data(i,27); % tail x
%             tailE.y_2(length(tailE.y_2)+1,1)=Data(i,28); % tail y
        end
    end

    % % interpolation
    % Bring location information if body points exist within 1s through 5s

    % nose
    [nose.x_1_missing, nose.y_1_missing, nose.x_1_smoothing, nose.y_1_smoothing]=smoothing_bodypoint(frame_1, nose.x_1, nose.y_1, 30, 150, NumberofFrame);
    [nose.x_2_missing, nose.y_2_missing, nose.x_2_smoothing, nose.y_2_smoothing]=smoothing_bodypoint(frame_2, nose.x_2, nose.y_2, 30, 150,  NumberofFrame);
    % body center
    [bodyC.x_1_missing, bodyC.y_1_missing, bodyC.x_1_smoothing, bodyC.y_1_smoothing]=smoothing_bodypoint(frame_1, bodyC.x_1, bodyC.y_1, 30, 150, NumberofFrame);
    [bodyC.x_2_missing, bodyC.y_2_missing, bodyC.x_2_smoothing, bodyC.y_2_smoothing]=smoothing_bodypoint(frame_2, bodyC.x_2, bodyC.y_2, 30, 150, NumberofFrame);
    % tail
    [tail.x_1_missing, tail.y_1_missing, tail.x_1_smoothing, tail.y_1_smoothing]=smoothing_bodypoint(frame_1, tail.x_1, tail.y_1, 30, 150, NumberofFrame);
    [tail.x_2_missing, tail.y_2_missing, tail.x_2_smoothing, tail.y_2_smoothing]=smoothing_bodypoint(frame_2, tail.x_2, tail.y_2, 30, 150, NumberofFrame);

    % distance (based on original data)
    for i=1:NumberofFrame
        frame_1_temp=find(frame_1==i); frame_2_temp=find(frame_2==i); 
        % frame information does exist!
        if length(frame_1_temp) == 1 & length(frame_2_temp) == 1 
            % each position information also exist! (female nose-male nose)
            if isnan(nose.x_1_smoothing(frame_1_temp,1)) == 0 & isnan(nose.y_1_smoothing(frame_1_temp,1)) == 0 & isnan(nose.x_2_smoothing(frame_2_temp,1)) == 0 & isnan(nose.y_2_smoothing(frame_2_temp,1)) == 0
                distance.nn(i,1) = sqrt((nose.x_1_smoothing(frame_1_temp,1)-nose.x_2_smoothing(frame_2_temp,1))^2+(nose.y_1_smoothing(frame_1_temp,1)-nose.y_2_smoothing(frame_2_temp,1))^2);
            else
                distance.nn(i,1) = NaN;
            end
            % female nose-male body
            if isnan(nose.x_1_smoothing(frame_1_temp,1)) == 0 & isnan(nose.y_1_smoothing(frame_1_temp,1)) == 0 & isnan(bodyC.x_2_smoothing(frame_2_temp,1)) == 0 & isnan(bodyC.y_2_smoothing(frame_2_temp,1)) == 0
                distance.nb(i,1) = sqrt((nose.x_1_smoothing(frame_1_temp,1)-bodyC.x_2_smoothing(frame_2_temp,1))^2+(nose.y_1_smoothing(frame_1_temp,1)-bodyC.y_2_smoothing(frame_2_temp,1))^2);
            else
                distance.nb(i,1) = NaN;
            end
            % female nose-male tail
            if isnan(nose.x_1_smoothing(frame_1_temp,1)) == 0 & isnan(nose.y_1_smoothing(frame_1_temp,1)) == 0 & isnan(tail.x_2_smoothing(frame_2_temp,1)) == 0 & isnan(tail.y_2_smoothing(frame_2_temp,1)) == 0
                distance.nt(i,1) = sqrt((nose.x_1_smoothing(frame_1_temp,1)-tail.x_2_smoothing(frame_2_temp,1))^2+(nose.y_1_smoothing(frame_1_temp,1)-tail.y_2_smoothing(frame_2_temp,1))^2);
            else
                distance.nt(i,1) = NaN;
            end
            % female body-male body
            if isnan(bodyC.x_1_smoothing(frame_1_temp,1)) == 0 & isnan(bodyC.y_1_smoothing(frame_1_temp,1)) == 0 & isnan(bodyC.x_2_smoothing(frame_2_temp,1)) == 0 & isnan(bodyC.y_2_smoothing(frame_2_temp,1)) == 0
                distance.bb(i,1) = sqrt((bodyC.x_1_smoothing(frame_1_temp,1)-bodyC.x_2_smoothing(frame_2_temp,1))^2+(bodyC.y_1_smoothing(frame_1_temp,1)-bodyC.y_2_smoothing(frame_2_temp,1))^2);
            else
                distance.bb(i,1) = NaN;
            end
            % female body-male tail
            if isnan(bodyC.x_1_smoothing(frame_1_temp,1)) == 0 & isnan(bodyC.y_1_smoothing(frame_1_temp,1)) == 0 & isnan(tail.x_2_smoothing(frame_2_temp,1)) == 0 & isnan(tail.y_2_smoothing(frame_2_temp,1)) == 0
                distance.bt(i,1) = sqrt((bodyC.x_1_smoothing(frame_1_temp,1)-tail.x_2_smoothing(frame_2_temp,1))^2+(bodyC.y_1_smoothing(frame_1_temp,1)-tail.y_2_smoothing(frame_2_temp,1))^2);
            else
                distance.bt(i,1) = NaN;
            end
            % female tail-male tail
            if isnan(tail.x_1_smoothing(frame_1_temp,1)) == 0 & isnan(tail.y_1_smoothing(frame_1_temp,1)) == 0 & isnan(tail.x_2_smoothing(frame_2_temp,1)) == 0 & isnan(tail.y_2_smoothing(frame_2_temp,1)) == 0
                distance.tt(i,1) = sqrt((tail.x_1_smoothing(frame_1_temp,1)-tail.x_2_smoothing(frame_2_temp,1))^2+(tail.y_1_smoothing(frame_1_temp,1)-tail.y_2_smoothing(frame_2_temp,1))^2);
            else
                distance.tt(i,1) = NaN;
            end
        else
            distance.nn(i,1) = NaN; distance.nb(i,1) = NaN; distance.nt(i,1) = NaN; distance.bb(i,1) = NaN; distance.bt(i,1) = NaN; distance.tt(i,1) = NaN;
        end
    end

    % Save distance.nn to a CSV file
    csvwrite('frame_distance_nn.csv', distance.nn); %jordan added to get nn dist data for each frame.

    %distance.mean_nn=nanmean(distance.nn,1);
    %distance.mean_nb=nanmean(distance.nb,1);
    %distance.mean_nt=nanmean(distance.nt,1);
    %distance.mean_bb=nanmean(distance.bb,1);
    %distance.mean_bt=nanmean(distance.bt,1);
    %distance.mean_tt=nanmean(distance.tt,1);

    % graph

    %plot(1:NumberofFrame, distance.nn, '.r'); %to change graph color: https://www.mathworks.com/help/matlab/creating_plots/specify-plot-colors.html
    %xlabel('time(frame)');
    %ylabel('nose distance difference(AU)');

    %cd([mother_root]);
    %print(['nose distance difference ' animal_id '.jpeg'], '-djpeg', '-r0');
    %close all;

    % save mat file

    save(['Z:\Jordan\sausi_3F_rd3\sleap_csv\' animal_id '.mat'], ...
    'nose','bodyC', 'tail');

    % save excel file
    %fprintf(outputCSV, '%d,%d,%d,%d,%d,%d\n', ...
    %    distance.mean_nn, distance.mean_nb, distance.mean_nt, distance.mean_bb, distance.mean_bt, distance.mean_tt);
    
end
% %

fclose all;

% load mat file

for animal_iter=1:length(Mice_list)

    animal_id=Mice_list{animal_iter};
    load(['Z:\Jordan\sausi_3F_rd3\sleap_csv\' animal_id '.mat']);

    %mouse 1
    %to add another bp, just increase to next colum number and change the
    %bp name to the correct one for each animal
    Interpolate_location{animal_iter}.A1(:,1)=nose.x_1_smoothing;
    Interpolate_location{animal_iter}.A1(:,2)=nose.y_1_smoothing;
    Interpolate_location{animal_iter}.A1(:,3)=bodyC.x_1_smoothing;
    Interpolate_location{animal_iter}.A1(:,4)=bodyC.y_1_smoothing;

    %mouse 2
    Interpolate_location{animal_iter}.A2(:,1)=nose.x_2_smoothing;
    Interpolate_location{animal_iter}.A2(:,2)=nose.y_2_smoothing;

    %on right panel, double click "Interpolar_location" variable. then
    %above in the variable tab, double click the first cell. then double
    %click th cell for A1 (animal 1) or A2 (animal 2) and this with give
    %the interpolated x & y body points. Ther are in the order of the
    %columns. Then just copy and paste into an excel file.

end

