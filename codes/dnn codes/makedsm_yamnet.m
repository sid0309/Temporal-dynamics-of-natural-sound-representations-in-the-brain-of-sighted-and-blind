%% Adapted from Giordano et al., 2023 Nature Neuroscience
%% convert dnn activations to common models format
%% and compute distances.
%% Run this code after having fit the DNN models
%% with analyze_04[a,b,c]*.py code
clear
rootmain='/mnt/DATA/siddharth/myenv/workspace/';

clc
restoredefaultpath
addpath(genpath([rootmain,'code/']))
addpath(genpath('/media/siddharth/DATA/CPP/Projects/Aud_Cat/codes/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind'))


% dataset='giordano';
% dataset='formisano';
% whichmodel=1;%kell
% whichmodel=2;%vggish
% whichmodel=3;%yamnet

dataset=1;
for whichmodel=3
    model_name={'kell' 'vggish' 'yamnet'}';
    model_name_nicer={'Kell' 'VGGish' 'Yamnet'};
    % out_dist_dir=[datdir,dataset,'_dnns/'];
    dnn_dir= '/media/siddharth/DATA/CPP/Projects/Aud_Cat/codes/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind/others/dnn/data/';
    main_out_fn='/media/siddharth/DATA/CPP/Projects/Aud_Cat/codes/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind/others/dnn/models/';
    
    %dnn_dir=dnn_dirs{whichmodel};
    d=dir([dnn_dir,'*.hdf5']);
    d=struct2cell(d);
    fns=d(1,:)';
    fns=celfun(@(x) [dnn_dir,x],fns);    
    
    if whichmodel==1
        layer_nams={'conv1' %1
            'rnorm1' %2
            'pool1' %3
            'conv2' %4
            'rnorm2' %5
            'pool2' %6
            'conv3' %7
            'conv4_W' %8
            'conv5_W' %9
            'pool5_W' %10
            'pool5_flat_W' %11
            'fc6_W' %12
            'conv_logits_W' %13
            'conv4_G' %14
            'conv5_G' %15
            'pool5_G' %16
            'pool5_flat_G' %17
            'fc6_G' %18
            'conv_logits_G'}; %19
        in_layers=1:length(layer_nams); %use to discard unwanted layers
        in_layers=[3 6 7 8 11 12 13 14 17 18 19];
        in_layers=[3 6 7 8 11 12 14 17 18]; %no embeddings
    elseif whichmodel==2
        %             layer_nams={'input' %1
        %                 'conv1relu' %2
        %                 'pool1mp' %3
        %                 'conv2relu' %4
        %                 'pool2mp' %5
        %                 'conv3_1relu' %6
        %                 'conv3_2relu' %7
        %                 'pool3mp' %8
        %                 'conv4_1relu' %9
        %                 'conv4_2relu' %10
        %                 'pool4mp' %11
        %                 'fc1_1matm' %12
        %                 'fc1_1relu' %13
        %                 'fc1_2matm' %14
        %                 'fc1_2relu' %15
        %                 'fc_2matm' %16
        %                 'embedding'}; %17
        
        layer_nams={'input_3' %1
            'conv1' %2
            'pool1'%3
            'conv2'%4
            'pool2'%5
            'conv3_1'%6
            'conv3_2'%7
            'pool3'%8
            'conv4_1'%9
            'conv4_2'%10
            'pool4'%11
            'flatten'%12
            'fc1_1'%13
            'fc1_2'%14
            'fc2'};%15
        in_layers=[3 5 8 11 13 14 15];
        
        
        
    elseif whichmodel==3
        layer_nams={'layer01relu' %1
            'layer02relu' %2
            'layer03relu' %3
            'layer04relu' %4
            'layer05relu' %5
            'layer06relu' %6
            'layer07relu' %7
            'layer08relu' %8
            'layer09relu' %9
            'layer10relu' %10
            'layer11relu' %11
            'layer12relu' %12
            'layer13relu' %13
            'layer14relu' %14
            'embedding' %15
            'logits' %16
            'predictions'}; %17
        in_layers=1:14; %no embeddings etc.
    else
        error('not a supported model')
    end
    layer_nams=layer_nams(in_layers);
    n_layers=length(layer_nams);
    
    
    
    RDMsEuc=cell([n_layers,1]);
    RDMsCos=cell([n_layers,1]);
    Components=[];
    ndims=[];
    disp(repmat(dataset,[2 5]))
    for j=1:length(layer_nams)
        for i=1:length(d)
            tmp=hdf5read(fns{i},layer_nams{j});
            if rem(i,5)==0
                str=['sound: ',num2str(i)];
                %                     disp(str)
            end
            %         str=[layer_nams{j},' size: ',num2str(size(tmp))];
            %         disp(str)
            if whichmodel==1 %for Kell (one single frame analyzed)
                %put a singletone time dimension first
                s=size(tmp,1,2,3,4);
                idx=find(s>1);
                otherdims=setxor(1:4,idx);
                permdims=[otherdims(1) idx otherdims(2:end)];
                tmp=permute(tmp,permdims);
            elseif whichmodel==2 %put time first
                s=size(tmp,1,2,3,4);
                
                if strcmp(dataset,'formisano')
                    %put all non-singleton dimensions first
                    idx=find(s>1);
                    permdims=[idx setxor(1:4,idx)];
                    tmp=permute(tmp,permdims);
                    
                    idx=(find(s==1,1,'first')); %first singleton is time
                    otherdims=setxor(1:4,idx);
                    permdims=[idx otherdims];
                    tmp=permute(tmp,permdims);
                elseif strcmp(dataset,'giordano')
                    idx=(find(s>1,1,'last')); %last_nonsingleton is time
                    otherdims=setxor(1:4,idx);
                    permdims=[idx otherdims];
                    tmp=permute(tmp,permdims);
                end
                %             str=[layer_nams{j},' size: ',num2str(size(tmp))];
            elseif whichmodel==3
                %s=size(tmp,1,2,3,4);
                s = size(tmp);
                idx=(find(s>1,1,'last')); %last_nonsingleton is time
                otherdims=setxor(1:4,idx);
                permdims=[idx otherdims];
                tmp=permute(tmp,permdims);
                tmp=tmp(1:end-1,:,:,:); %discard last analysis frame
                %             str=[layer_nams{j},' size: ',num2str(size(tmp))];
            end
            if i==1
                dat_tmp=zeros(size(repmat(tmp,[1 1 1 1 length(d)])));
                %             disp(num2str(size(dat_tmp)))
                str=[layer_nams{j},' size: ',num2str(size(dat_tmp))];
                disp(str)
            end
            dat_tmp(:,:,:,:,i)=tmp;
            
            
            
        end
        
        dat_tmp=squeeze(mean(dat_tmp,1)); %average across analysis frames
        dat_tmp=reshape(dat_tmp,[size(dat_tmp,1)*size(dat_tmp,2)*size(dat_tmp,3)],24);
        for ii = 1:24
            for jj = 1:24
                                    dsm(ii,jj,j) = pdist2(dat_tmp(:,ii)',dat_tmp(:,jj)','cosine');
            end
        end
        
        dsm = dsm./max(max(dsm)); % normalize to 0 to 1
        
        
        %s=size(dat_tmp);
        %thisdat_tmp=reshape(dat_tmp,[prod(s(1:4)) s(5)]);
        %RDMsEuc{end+1,1}=BLG_EucDistND(thisdat_tmp);
        %RDMsCos{end+1,1}=BLG_CosDistND(thisdat_tmp);
        %Components{end+1,1}=[layer_nams{j}];
        %ndims=cat(1,ndims,size(thisdat_tmp,1));
    end
    
    
    % Model=model_name_nicer(whichmodel);
    %Components=Components(:)';
    
    %%% save the distances
    %D=cell2mat(RDMsEuc');
    %fnout=[main_out_fn,'dist_euc.mat'];
    %         save(fnout,'Model','Components','D','ndims','-v7.3')
    %disp(fnout)
    %D=cell2mat(RDMsCos');
    %fnout=[main_out_fn,'dist_cos.mat'];
    %         save(fnout,'Model','Components','D','ndims','-v7.3')
    %disp(fnout)
    %disp('all done for this model!')
    %disp(repmat([dataset,' ',Model{1},' '],[2 10]))
    
    
end
save('/media/siddharth/DATA/CPP/Projects/Aud_Cat/codes/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind/others/dnn/models/dsm_dsm.mat',dsm);

%%

% giordano Kell
% pool1 size: 1  96  43  43  80
% pool2 size: 1  256   11   11   80
% conv3 size: 1  512   11   11   80
% conv4_W size: 1  1024    11    11    80
% pool5_flat_W size: 1  18432      1      1     80
% fc6_W size: 1  1024     1     1    80
% conv4_G size: 1  1024    11    11    80
% pool5_flat_G size: 1  18432      1      1     80
% fc6_G size: 1  1024     1     1    80
%
%
% giordano VGGish
% pool1 size: 2  64  32  48  80
% pool2 size: 2  128   16   24   80
% pool3 size: 2  256    8   12   80
% pool4 size: 2  512    4    6   80
% fc1_1 size: 2  4096     1     1    80
% fc1_2 size: 2  4096     1     1    80
% fc2 size: 2  128    1    1   80
%
%
% giordano Yamnet
% layer01relu size: 2  32  32  48  80
% layer02relu size: 2  64  32  48  80
% layer03relu size: 2  128   16   24   80
% layer04relu size: 2  128   16   24   80
% layer05relu size: 2  256    8   12   80
% layer06relu size: 2  256    8   12   80
% layer07relu size: 2  512    4    6   80
% layer08relu size: 2  512    4    6   80
% layer09relu size: 2  512    4    6   80
% layer10relu size: 2  512    4    6   80
% layer11relu size: 2  512    4    6   80
% layer12relu size: 2  512    4    6   80
% layer13relu size: 2  1024     2     3    80
% layer14relu size: 2  1024     2     3    80
%
%
% formisano Kell
% pool1 size: 1   96   43   43  288
% pool2 size: 1  256   11   11  288
% conv3 size: 1  512   11   11  288
% conv4_W size: 1  1024    11    11   288
% pool5_flat_W size: 1  18432      1      1    288
% fc6_W size: 1  1024     1     1   288
% conv4_G size: 1  1024    11    11   288
% pool5_flat_G size: 1  18432      1      1    288
% fc6_G size: 1  1024     1     1   288
%
%
% formisano VGGish
% pool1 size: 1   64   32   48  288
% pool2 size: 1  128   16   24  288
% pool3 size: 1  256    8   12  288
% pool4 size: 1  512    4    6  288
% fc1_1 size: 1  4096     1     1   288
% fc1_2 size: 1  4096     1     1   288
% fc2 size: 1  128    1    1  288
%
%
% formisano Yamnet
% layer01relu size: 1   32   32   48  288
% layer02relu size: 1   64   32   48  288
% layer03relu size: 1  128   16   24  288
% layer04relu size: 1  128   16   24  288
% layer05relu size: 1  256    8   12  288
% layer06relu size: 1  256    8   12  288
% layer07relu size: 1  512    4    6  288
% layer08relu size: 1  512    4    6  288
% layer09relu size: 1  512    4    6  288
% layer10relu size: 1  512    4    6  288
% layer11relu size: 1  512    4    6  288
% layer12relu size: 1  512    4    6  288
% layer13relu size: 1  1024     2     3   288
% layer14relu size: 1  1024     2     3   320