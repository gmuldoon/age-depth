function [paramRange,nparam] = setParams2(accumFlag,H,pik)
    lp = length(pik);

    if accumFlag == 7
        
        % set prior distributions on these parameters
        paramRange(1,1:2) =  [0.0 1.0];              % Schwander s
        
%         for i = 2:40
%             paramRange(i,1:2) = [0.05 0.15];
%         end
%         paramRange(2,1:2) =  [0.05 0.20];           % accum 1294m < depth < 2191m
%         paramRange(3,1:2) =  [0.05 0.20];           % accum 1024m < depth < 1294m
%         paramRange(4,1:2) =  [0.03 0.20];          % accum 150 m < depth < 1024m
%         paramRange(5,1:2) =  [0.05 0.20];           % accum depth < 150m
%         paramRange(6,1:2) =  [0.05 0.20];
%         paramRange(7,1:2) =  [0.05 0.14];
%         paramRange(8,1:2) =  [0.05 0.20];
%         paramRange(9,1:2) =  [0.05 0.20];
%         paramRange(10,1:2) = [0.05 0.14];
%         paramRange(11,1:2) = [0.05 0.14];
        paramRange(2,1:2) =  [0.01 0.25];           % bottom of the ice
        paramRange(3,1:2) =  [0.01 0.25];           
        paramRange(4,1:2) =  [0.01 0.25];          
        paramRange(5,1:2) =  [0.01 0.25];           
        paramRange(6,1:2) =  [0.01 0.25];
        paramRange(7,1:2) =  [0.01 0.25];
        paramRange(8,1:2) =  [0.01 0.25];
        paramRange(9,1:2) =  [0.01 0.25];
        paramRange(10,1:2) = [0.01 0.25];
        paramRange(11,1:2) = [0.01 0.25];           % top of the ice
        
        paramRange(12,1:2)=[0.01 0.5];            % Nye h
        paramRange(13,1:2)=[1.68*10^8 1.695*10^8]; %velocity in glacial ice
        paramRange(14,1:2)=[8.6572 8.6574];             %firn correction (diff in ice thickness between ice column with/out firn
        
         % determine number of parameters
        nparam = length(paramRange(:,1))+lp;
        n_schwander_params = nparam - lp; % Dr, S, Schwander s, h, accum levels
        
         for i = 1:lp
            %paramRange(i+n_schwander_params,1:2) = [0 H];
            v_est = 1.68*10^8;
            pikTWTT_sec = pik(i)*10^-6;
            depth_est = v_est*pikTWTT_sec/2;
            paramRange(i+n_schwander_params,1:2) = [depth_est-0.2*H depth_est+.1*H];
        end
        
               
    end
 
end
%     prevParam(2,1) = 0.14; %deepest part
%     prevParam(3,1) = 0.14;
%     prevParam(4,1) = 0.06;
%     prevParam(5,1) = 0.0805;
%     prevParam(6,1) = 0.12;
%     prevParam(7,1) = 0.1275;
%     prevParam(8,1) = 0.125;
%     prevParam(9:11,1)= 0.14; %shallowest part