function [paramRange,nparam] = setParams2(accumFlag,H,pik)
    lp = length(pik);

    if accumFlag == 7
        
        % set prior distributions on these parameters
%         paramRange(1,1:2)=[0.5 1.];              % Schwander s
%         paramRange(2,1:2)=[0.07 0.005];           % accum 1294m < depth < 2191m
%         paramRange(3,1:2)=[0.07 0.005];           % accum 1024m < depth < 1294m
%         paramRange(4,1:2)=[0.07 0.005];           % accum 150 m < depth < 1024m
%         paramRange(5,1:2)=[0.07 0.005];           % accum depth < 150m    %mean and std
%         paramRange(6,1:2) =  [0.10 0.005];
%         paramRange(7,1:2) =  [0.10 0.005];
%         paramRange(8,1:2) =  [0.10 0.005];
%         paramRange(9,1:2) =  [0.10 0.005];
%         paramRange(10,1:2) = [0.10 0.005];
        paramRange(1,1:2) =  [0.5 1.0];              % Schwander s
        
        for i = 2:40
            paramRange(i,1:2) = [0.05 0.20];
        end
%         paramRange(2,1:2) =  [0.01 0.15];           % accum 1294m < depth < 2191m
%         paramRange(3,1:2) =  [0.01 0.15];           % accum 1024m < depth < 1294m
%         paramRange(4,1:2) =  [0.01 0.15];          % accum 150 m < depth < 1024m
%         paramRange(5,1:2) =  [0.01 0.15];           % accum depth < 150m
%         paramRange(6,1:2) =  [0.01 0.15];
%         paramRange(7,1:2) =  [0.01 0.15];
%         paramRange(8,1:2) =  [0.01 0.15];
%         paramRange(9,1:2) =  [0.01 0.15];
%         paramRange(10,1:2) = [0.01 0.15];
%         paramRange(11,1:2) = [0.01 0.20];
%         paramRange(12,1:2) = [0.01 0.20];
%         paramRange(13,1:2) = [0.01 0.20];
%         paramRange(14,1:2) = [0.01 0.20];
%         paramRange(15,1:2) = [0.01 0.20];
%         paramRange(16,1:2) = [0.01 0.20];
%         paramRange(17,1:2) = [0.01 0.20];
%         paramRange(18,1:2) = [0.01 0.20];
%         paramRange(19,1:2) = [0.01 0.20];
%         paramRange(20,1:2) = [0.01 0.20];
%         paramRange(21,1:2) = [0.01 0.20];
%         paramRange(22,1:2) = [0.01 0.20];
%         paramRange(23,1:2) = [0.01 0.20];
%         paramRange(24,1:2) = [0.01 0.20];
%         paramRange(25,1:2) = [0.01 0.20];
%         paramRange(26,1:2) = [0.01 0.20];
%         paramRange(27,1:2) = [0.01 0.20];
%         paramRange(28,1:2) = [0.01 0.20];
%         paramRange(29,1:2) = [0.01 0.20];
%         paramRange(30,1:2) = [0.01 0.20];
%         paramRange(31,1:2) = [0.01 0.20];
%         paramRange(32,1:2) = [0.01 0.20];
%         paramRange(33,1:2) = [0.01 0.20];
%         paramRange(34,1:2) = [0.01 0.20];
%         paramRange(35,1:2) = [0.01 0.20];
%         paramRange(36,1:2) = [0.01 0.20];
%         paramRange(37,1:2) = [0.01 0.20];
%         paramRange(38,1:2) = [0.01 0.20];
%         paramRange(39,1:2) = [0.01 0.20];
%         paramRange(40,1:2) = [0.01 0.20];

        


        paramRange(41,1:2)=[0.0 0.5];            % Nye h
        paramRange(42,1:2)=[1.68*10^8 1.695*10^8]; %velocity in glacial ice
        paramRange(43,1:2)=[1 65];             %firn correction (diff in ice thickness between ice column with/out firn
        

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
