function BlockOutsCombination(BlockOutFileLocation,BlockOutFileName,SaveName)
%% Generate Combined BlockOuts
load(fullfile(BlockOutFileLocation,BlockOutFileName))
%% NaN Matrices
MatrixVariableList = {'Chi2','AOD','SSA','Cv','f1','reff1','reff2','mr','mi','PSD'};
MatrixVariableLength = [1,6,6,1,1,1,1,6,6,22];
StructVariableList = {'DataBRF','ModelBRF'};
VariableList = [MatrixVariableList,StructVariableList];

for i = 1:length(MatrixVariableList)
    for j = 1:3
        eval( [MatrixVariableList{i},'_',num2str(j),' = ones(',num2str(MatrixVariableLength(i)),',size(Blocks_Info.latlon_Norm,1))*nan;' ] );
    end
end

for i = 1:length(StructVariableList)
    for j = 1:3
        eval( [StructVariableList{i},'_',num2str(j),' = cell(1,size(Blocks_Info.latlon_Norm,1));' ] );
    end
end

%% Combination
BlockOut_1 = Blocks_Info_1.BlockOut;
for i = 1:size(Blocks_Info.IsValid,1) %2:2%
    for j = 1:size(Blocks_Info.IsValid,2)%7:7%
%         if isfield(Blocks_Info_1,'BlockOut_1')

            if ~isempty(BlockOut_1{i,j})
                if isfield(BlockOut_1{i,j},'AOD')
                    CurrentIndexList = reshape(Blocks_Info.InitialPixelIndex{i,j}',1,[]);
                    IsValidPixel = reshape(BlockOut_1{i,j}.Patch_valid_mat',1,[]);
                    TempValidIndex = find(IsValidPixel);
                    Chi2_1(:,CurrentIndexList(TempValidIndex)) = BlockOut_1{i,j}.Chi_sq;
                    AOD_1(:,CurrentIndexList(TempValidIndex)) = BlockOut_1{i,j}.AOD;
                    SSA_1(:,CurrentIndexList(TempValidIndex)) = BlockOut_1{i,j}.SSA;
                    Cv_1(:,CurrentIndexList(TempValidIndex)) = exp(cell2mat(BlockOut_1{i,j}.pout(1,:)));
                    f1_1(:,CurrentIndexList(TempValidIndex)) = BlockOut_1{i,j}.f1out;
                    reff1_1(:,CurrentIndexList(TempValidIndex)) = exp(cell2mat(BlockOut_1{i,j}.pout(5,:)));
                    reff2_1(:,CurrentIndexList(TempValidIndex)) = exp(cell2mat(BlockOut_1{i,j}.pout(7,:)));
                    mr_1(:,CurrentIndexList(TempValidIndex)) = BlockOut_1{i,j}.mr1out;
                    mi_1(:,CurrentIndexList(TempValidIndex)) = BlockOut_1{i,j}.mi1out;
                    PSD_1(:,CurrentIndexList(TempValidIndex)) = BlockOut_1{i,j}.PSD;
                    for i_Temp = 1:length(TempValidIndex)
                        DataBRF_1{CurrentIndexList(i_Temp)} = BlockOut_1{i,j}.dataBRF{i_Temp};
                        ModelBRF_1{CurrentIndexList(i_Temp)} = BlockOut_1{i,j}.modelBRF{i_Temp};
                    end
                    
                end
            end
%         end

%         if isfield(Blocks_Info_2,'BlockOut_2')
        if exist('Blocks_Info_2',"var")
            if ~isempty(Blocks_Info_2.BlockOut{i,j})
                if isfield(Blocks_Info_2.BlockOut{i,j},'AOD')
                    CurrentIndexList = reshape(Blocks_Info.InitialPixelIndex{i,j}',1,[]);
                    IsValidPixel = reshape(Blocks_Info_2.BlockOut{i,j}.Patch_valid_mat',1,[]);
                    TempValidIndex = find(IsValidPixel);
                    Chi2_2(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_2.BlockOut{i,j}.Chi_sq;
                    AOD_2(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_2.BlockOut{i,j}.AOD;
                    SSA_2(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_2.BlockOut{i,j}.SSA;
                    Cv_2(:,CurrentIndexList(TempValidIndex)) = exp(cell2mat(Blocks_Info_2.BlockOut{i,j}.pout(1,:)));
                    f1_2(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_2.BlockOut{i,j}.f1out;
                    reff1_2(:,CurrentIndexList(TempValidIndex)) = exp(cell2mat(Blocks_Info_2.BlockOut{i,j}.pout(5,:)));
                    reff2_2(:,CurrentIndexList(TempValidIndex)) = exp(cell2mat(Blocks_Info_2.BlockOut{i,j}.pout(7,:)));
                    mr_2(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_2.BlockOut{i,j}.mr1out;
                    mi_2(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_2.BlockOut{i,j}.mi1out;
                    PSD_2(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_2.BlockOut{i,j}.PSD;
                    for i_Temp = 1:length(TempValidIndex)
                        DataBRF_2{CurrentIndexList(i_Temp)} = Blocks_Info_2.BlockOut{i,j}.dataBRF{i_Temp};
                        ModelBRF_2{CurrentIndexList(i_Temp)} = Blocks_Info_2.BlockOut{i,j}.modelBRF{i_Temp};
                    end
                end
            end
        end

%         if isfield(Blocks_Info_3,'BlockOut_3')
        if exist('Blocks_Info_3',"var")
            if ~isempty(Blocks_Info_3.BlockOut{i,j})
                if isfield(Blocks_Info_3.BlockOut{i,j},'AOD')
                    CurrentIndexList = reshape(Blocks_Info.InitialPixelIndex{i,j}',1,[]);
                    IsValidPixel = reshape(Blocks_Info_3.BlockOut{i,j}.Patch_valid_mat',1,[]);
                    TempValidIndex = find(IsValidPixel);
                    Chi2_3(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_3.BlockOut{i,j}.Chi_sq;
                    AOD_3(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_3.BlockOut{i,j}.AOD;
                    SSA_3(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_3.BlockOut{i,j}.SSA;
                    Cv_3(:,CurrentIndexList(TempValidIndex)) = exp(cell2mat(Blocks_Info_3.BlockOut{i,j}.pout(1,:)));
                    f1_3(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_3.BlockOut{i,j}.f1out;
                    reff1_3(:,CurrentIndexList(TempValidIndex)) = exp(cell2mat(Blocks_Info_3.BlockOut{i,j}.pout(5,:)));
                    reff2_3(:,CurrentIndexList(TempValidIndex)) = exp(cell2mat(Blocks_Info_3.BlockOut{i,j}.pout(7,:)));
                    mr_3(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_3.BlockOut{i,j}.mr1out;
                    mi_3(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_3.BlockOut{i,j}.mi1out;
                    PSD_3(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_3.BlockOut{i,j}.PSD;
                    PSD_3(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_3.BlockOut{i,j}.PSD;
                    PSD_3(:,CurrentIndexList(TempValidIndex)) = Blocks_Info_3.BlockOut{i,j}.PSD;
                    for i_Temp = 1:length(TempValidIndex)
                        DataBRF_3{CurrentIndexList(i_Temp)} = Blocks_Info_3.BlockOut{i,j}.dataBRF{i_Temp};
                        ModelBRF_3{CurrentIndexList(i_Temp)} = Blocks_Info_3.BlockOut{i,j}.modelBRF{i_Temp};
                    end
                end
            end
        end

    end
end
% Cv_1_all = exp(Cv_1_all);
% Cv_2_all = exp(Cv_2_all);
% Cv_3_all = exp(Cv_3_all);
%% Generate Odd / Even index lists
IndexList_even = [];
IndexList_odd = [];

for i = 1:size(Blocks_Info.IsValid,1)
    if rem(i,2) == 1
        for j = 1:2:size(Blocks_Info.IsValid,2)
            CurrentIndexList = reshape(Blocks_Info.InitialPixelIndex{i,j}',1,[]);
            CurrentIndexList = CurrentIndexList(~isnan(CurrentIndexList));
            IndexList_even = [IndexList_even,CurrentIndexList];
        end
        for j = 2:2:size(Blocks_Info.IsValid,2)
            CurrentIndexList = reshape(Blocks_Info.InitialPixelIndex{i,j}',1,[]);
            CurrentIndexList = CurrentIndexList(~isnan(CurrentIndexList));
            IndexList_odd = [IndexList_odd,CurrentIndexList];
        end
    else
        for j = 2:2:size(Blocks_Info.IsValid,2)
            CurrentIndexList = reshape(Blocks_Info.InitialPixelIndex{i,j}',1,[]);
            CurrentIndexList = CurrentIndexList(~isnan(CurrentIndexList));
            IndexList_even = [IndexList_even,CurrentIndexList];
        end
        for j = 1:2:size(Blocks_Info.IsValid,2)
            CurrentIndexList = reshape(Blocks_Info.InitialPixelIndex{i,j}',1,[]);
            CurrentIndexList = CurrentIndexList(~isnan(CurrentIndexList));
            IndexList_odd = [IndexList_odd,CurrentIndexList];
        end
    end
end

%% Save .mat of Ordered Variables

% FieldNameList = {'Chi2','AOD','SSA','Cv','f1','reff1','reff2','mr','mi','PSD','DataBRF','ModelBRF'};
for i = 1:length(VariableList)
    for j = 1:3
        xuout.([VariableList{i},'_',num2str(j)]) = eval( [VariableList{i},'_',num2str(j),';' ] );
    end
end

save(fullfile(BlockOutFileLocation,SaveName),'xuout','VariableList','IndexList_even','IndexList_odd')
