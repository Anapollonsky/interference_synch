classdef wififrame
    % wififrame: Simulates the properties of a frame.
    %   Detailed explanation goes here
    
    properties
        OFDMLen = 48;
        DataLength        % Length of data in binary
        MaxDataLength = 2304*8;     % Maximum binary length of data
        CyclicPrefixLength = 16;
        FrameType = 'Data';
        TrainingTimeLen = 12;
        TransmissionRate  % 1-8
        Mod               % Number of points in constellation
        CodeRate          % Code Rate
        Training          % Parallel
        Preamble          % Series
        Data
        TrainingMod
        PreambleMod
        DataMod
        SignalParallel
    end
    
    methods
        function frame = wififrame(data, tr, training)
            frame.TransmissionRate = tr;      
            
            %% Basic Allocation
            if tr < 3
                frame.Mod = 2;
            elseif tr < 5
                frame.Mod = 4;
            elseif tr < 7
                frame.Mod = 16;
            else
                frame.Mod = 64;
            end
          
            if tr == 1 || tr == 3 || tr == 5
                frame.CodeRate = 1/2;
            elseif tr == 7
                frame.CodeRate = 2/3;
            else
                frame.CodeRate = 3/4;
            end
%             frame.FFTPoints = (1728*8/48 + 13) * frame.CodeRate;
            
            %% Signal Construction
            frame.DataLength = length(data);
            frame.Training = training;
            frame.Preamble = [dec2bin(tr-1, 3) - 48 ...
                dec2bin(frame.DataLength/8,16)-48 randi([0,1], 1, 5 + 48)] + 0;
            
            frame.Data = data;
%             frame.DataMod = wifimod(frame.Data, tr);
%             frame.HeaderMod = [reshape(wifimod(reshape(frame.Training, ...
%                 [], 1), 1), [], frame.OFDMLen); wifimod(frame.Preamble.', 1).'];
%             frame.SignalParallel = [frame.HeaderMod; reshape(frame.DataMod,...
%                 [], frame.OFDMLen)];

            frame.TrainingMod = reshape(wifimod(reshape(frame.Training, [], 1), 1),...
                [], frame.OFDMLen);
            frame.PreambleMod = reshape(wifimod(wifint(wifienc(frame.Preamble.', ...
                1)), 1), [], frame.OFDMLen);
            frame.DataMod = reshape(wifimod(wifint(wifienc(frame.Data, ...
                tr)), tr), [], frame.OFDMLen);
            frame.SignalParallel = ...
                [frame.TrainingMod; frame.PreambleMod; frame.DataMod];
            
        end
    end
end


