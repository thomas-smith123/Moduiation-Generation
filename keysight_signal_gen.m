classdef keysight_signal_gen
    properties
        fs
        resample_fs
        max_signal_per_frame
        resolution
        window_length
        overlap
        noverlap
        path
        snr_
        number
        digital_phase_modulate_class1
        pulse_modulate_class1
        pulse_modulate_class2
        others
        all_modulation
        framePoints
        check
        seq
    end

    methods
        %% 构造函数
        function self = keysight_signal_gen(fs,resample_fs,max_signal_per_frame,resolution,window_length,overlap,path,snr_,number,check)
            self.fs = fs;
            self.resample_fs = resample_fs;
            self.max_signal_per_frame = max_signal_per_frame;
            self.resolution = resolution;
            self.window_length = window_length;
            self.overlap = overlap;
            self.path = path;
            self.snr_ = snr_;
            self.number = number;
            self.check = check;

            self.digital_phase_modulate_class1=["BPSK","BPSK_X","QPSK","OQPSK","8-PSK","QAM16","QAM32","QAM64","QAM128","QAM256","QAM512","APSK16","APSK32","PAM4","QAM8","OOK"];
            % self.digital_phase_modulate_class1=[];
            % self.pulse_modulate_class1= [];
            self.pulse_modulate_class1 = ["none","increasing","decreasing","v-shape","inverted v","barker-2 +-","barker-2 ++","barker-3","barker-4 ++-+","barker-4 +++-","barker-5","barker-7","barker-11","barker-13","frank-4","frank-6"];
            self.pulse_modulate_class2 = ["2FSK","costas","frequency_hop","nlfm_sine",...
            %"nlfm_atan", 还没整，但是可以搞
            ];
            % self.others = [ ];
            self.others = [ "adsb","dsb" ];
            self.all_modulation = [self.digital_phase_modulate_class1,self.pulse_modulate_class1,self.pulse_modulate_class2,self.others];
            self.noverlap = resolution(1)-overlap;
            self.framePoints = (window_length-overlap)*(resolution(2)-1)+self.window_length;
        end
%% 判断文件夹是否存在
        function check_path(self)
            if ~exist([self.path,'raw_complex/'], 'dir')
                mkdir([self.path,'raw_complex/']);
            end
            if ~exist([self.path,'labels/'], 'dir')
                mkdir([self.path,'labels/']);
            end
            if ~exist([self.path,'images/'], 'dir')
                mkdir([self.path,'images/']);
            end
            if ~exist([self.path,'check/'], 'dir')
                mkdir([self.path,'check/']);
            end
        end
%% 生成随机数
        function [rand_u] = rand_uniform(self,low,high)
            rand_u = low + (high-low)*rand;
        end
%% 数组移位
        function array_shifted = shift_with_zero(self,array, shiftAmount, cval)
            % 创建一个和原数组同样大小的0数组
            zeroArray = zeros(size(array));
            
            % 如果左移的位数大于数组的长度，返回0数组
            if abs(shiftAmount) >= length(array)
                array_shifted = zeroArray;
            else
                % 左移数组
                array_shifted = circshift(array, -shiftAmount);
                
                % 用0填充空位
                if shiftAmount > 0
                    array_shifted(end-shiftAmount+1:end) = 0;
                else
                    array_shifted(1:abs(shiftAmount)) = 0;
                end
            end
        end
 %% 生成
        function [time_domain_signal,label] = get_frame(self)
            
            time_domain_signal = zeros(self.framePoints,1);
            label = []; % start, stop, bandwidth, centerFreq, modType
            % 抽取调制方式
            num = randi([1,self.max_signal_per_frame]);
            if length(self.all_modulation)<num
                num=length(self.all_modulation);
            end
            index = randperm(length(self.all_modulation),num);
            choosed_signal = self.all_modulation(index);
            % 开始产生
            for tmp_i = 1:length(choosed_signal)
                i = char(choosed_signal(tmp_i));
                startPoint = 1;
                endPoint = self.framePoints;
                flag = 2;
                if ismember(i,self.digital_phase_modulate_class1) 
                    
                    numSymbol = randi([5,10])*1e2;                    
                    centerFreq = self.rand_uniform(-self.resample_fs/2*0.9,self.resample_fs/2*0.9);
                    
                    symbolrate = randi([1,20])*1e7;
                    oversampling = ceil(self.fs/symbolrate);
                    position = randi([0,3]); % 0: 前不足，至最后; 1：全满; 2: 后不足，至最前
                    
                    bandwidth = (1+0.35)*symbolrate;
                    % 产生信号
                    
                    [iqdata, a_,a_,a_,a_] = iqmod('sampleRate', self.fs, 'numSymbols', numSymbol, 'data', 'Random', 'modType', i,...
                        'oversampling', oversampling, 'filterType', 'Root Raised Cosine', 'filterNsym', 80.0, 'filterBeta', 0.35,...
                        'carrierOffset', centerFreq, 'quadErr', 0.0, 'iqskew', 0.0, 'gainImbalance', 0.0, 'channelMapping', [1,0]);

                    iqdata = resample(iqdata, self.resample_fs, self.fs); % 这个是我们需要的
                    
                    tmp_time_domain = zeros(self.framePoints,1);
                    if size(iqdata,1) >= self.framePoints % 如果信号长度大于一帧的长度
                        if position == 0
                            startPoint = randi([0,ceil(self.framePoints*0.5)]);
                            tmp_time_domain(startPoint:self.framePoints,:) = iqdata(startPoint:self.framePoints,:);
                            iqdata = tmp_time_domain;
                            
                        elseif position == 1 %中间截取，填满
                            startPoint = randi([0,(size(iqdata,1)-self.framePoints)]); % 信号截取的位置
                            iqdata = iqdata(startPoint:self.framePoints+startPoint-1,:);
                            startPoint = 1;endPoint = self.framePoints;
                        else
                            endPoint = randi([ceil(self.framePoints*0.3),ceil(self.framePoints)]); % 信号截取的位置
                            tmp_time_domain(1:endPoint,1:size(tmp_time_domain,2)) = iqdata(size(iqdata,1)-endPoint+1:size(iqdata,1),1:size(iqdata,2));
                            iqdata = tmp_time_domain;
                        end

                    else % 如果信号长度小于一帧的长度
                        tmp = floor((self.framePoints-size(iqdata,1))/2);
                        startPoint = ceil(self.rand_uniform(0,tmp)); % 信号截取的位置
                        tmp = zeros(self.framePoints,1);
                        if startPoint+size(iqdata,1) <= self.framePoints
                            tmp(startPoint:startPoint+size(iqdata,1)-1,:) = iqdata;
                            endPoint = startPoint+size(iqdata,1);
                        else
                            tmp(startPoint:self.framePoints-1,1:size(tmp,2)) = iqdata(1:self.framePoints-startPoint);
                            endPoint = self.framePoints;
                        end
                        iqdata = tmp;
                    end

                    % 频率: signalFreq; 起始：startPoint; 终止：endPoint; 带宽：bandwidth; 调制方式：i
                elseif ismember(i,self.pulse_modulate_class1)% #@audit 这里可能有问题
                    t = self.framePoints/self.resample_fs;
                    % 随机产生一些参数
                    centerFreq = self.rand_uniform(-self.resample_fs/2*0.9,self.resample_fs/2*0.9);
                    bandwidth = self.rand_uniform(5e7,5e8);
                    amp = self.rand_uniform(-10,10);
                    
                    pw = self.rand_uniform(self.framePoints/self.resample_fs*0.25,...
                        self.framePoints/self.resample_fs*0.8);
                    pri = self.rand_uniform(1.5,3)*pw;
                    % position = 1 ## 0: 前移; 1：后移
                    position = randi([0,1]); %# 0: 前移; 1：后移
                    
                    if ismember(i,["increasing","decreasing","v-shape","inverted v"])
                        if self.resample_fs/2-abs(centerFreq)<bandwidth/2
                            centerFreq=(self.resample_fs/2-bandwidth/2)*centerFreq/abs(centerFreq);
                        end
                        [iqdata, marker, numRepeats, chMap] = iqpulse('sampleRate', 64e9, 'PRI', pri, 'PW', pw, 'pulseShape',...
                            'Raised Cosine', 'span', bandwidth, 'offset', centerFreq, 'amplitude', amp, 'modulationType', i,...
                            'correction', 0, 'channelMapping', [1,0], 'normalize', 1);
                        %iqdata = conj(iqdata);
                    
                    else % 脉冲相位调制
                        pri = self.rand_uniform(1e-6,10e-6);
                        pw = self.rand_uniform(0.4,0.8)*pri;
                        bandwidth = 5e7;
                        [iqdata, marker, numRepeats, chMap] = iqpulse('sampleRate', 64e9, 'PRI', pri, 'PW', pw, 'pulseShape', 'Raised Cosine',...
                            'span', bandwidth, 'offset', centerFreq, 'amplitude', amp, 'modulationType', i, 'correction', 0, 'phase', 0,...
                            'channelMapping', [1,0], 'normalize', 1);
                    end
                    iqdata = resample(iqdata, self.resample_fs, self.fs); %# 这个是我们需要的
                    
                    %# 长度补齐
                    if max(size(iqdata))-self.resample_fs*(pri-pw) < self.framePoints % 如果信号长度小于一帧的长度
                        repeat_times = ceil(ceil(self.framePoints/size(iqdata,1))); % 必须重复
                        iqdata = repmat(iqdata,repeat_times+1,1);
                        %iqdata = np.tile(iqdata[:,0],repeat_times); %注意一下维度
                        flag=1;
                    else
                        flag=0;
                        if rand()<0.5
                            iqdata = conj(iqdata(1:self.framePoints))';
                        else
                            iqdata = conj(iqdata(max(size(iqdata))-self.framePoints+1:max(size(iqdata))))';
                        end
                    end

                    % #@audit 时移
                    if flag == 0
                        if position == 0 %左移
                            startPoint = randi([0,ceil(self.framePoints*0.7)]); % 信号截取的位置
                            if length(size(iqdata))==1
                                iqdata = self.shift_with_zero(iqdata, -startPoint, 0);
                            else 
                                iqdata = self.shift_with_zero(reshape(iqdata,1,[]), -startPoint, 0);
                            end
                            if ceil(startPoint+max(size(iqdata))-(pri-pw)*self.resample_fs) <= self.framePoints
                                endPoint = startPoint+max(size(iqdata))-(pri-pw)*self.resample_fs;
                            else
                                endPoint = self.framePoints;
                            end
                        elseif position == 1 %右移
                            endPoint = randi([0,ceil(self.resample_fs*pw/2+(pri-pw)*self.resample_fs)]); % 信号截取的位置
                            if length(size(iqdata))==1
                                iqdata = self.shift_with_zero(iqdata, endPoint, 0);
                            else 
                                iqdata = self.shift_with_zero(reshape(iqdata,1,[]), endPoint, 0);
                            end
                            endPoint = self.framePoints-(endPoint+ceil((pri-pw)*self.resample_fs));
                        end
                        iqdata = conj(iqdata');
                    else
                        if position == 0
                            startPoint = randi([0,ceil(self.resample_fs*pw/2)]); % 信号截取的位置
                            if length(size(iqdata))==1
                                iqdata = self.shift_with_zero(iqdata, -startPoint, 0);
                            else 
                                iqdata = self.shift_with_zero(reshape(iqdata,1,[]), -startPoint, 0);
                            end
                            
                            tmp_ = (self.framePoints-startPoint)/(pri*self.resample_fs);
                            if (tmp_-floor(tmp_))*(pri) > pw %移动后为空
                                endPoint = startPoint+floor(tmp_)*self.resample_fs*pri+pw*self.resample_fs;
                            end
                            iqdata = iqdata(1:self.framePoints);
                        else
                            endPoint = randi([0,ceil(self.framePoints/2)]); % 信号截取的位置)
                            iqdata = iqdata(1:max(size(iqdata))-ceil(self.resample_fs*(pri-pw))-1);
                            if max(size(iqdata))-ceil(self.resample_fs*(pri-pw))<self.framePoints
                                if length(size(iqdata))==1
                                    iqdata = self.shift_with_zero(iqdata, endPoint, 0); 
                                else 
                                    iqdata = self.shift_with_zero(reshape(iqdata,1,[]), endPoint, 0);
                                end
                                endPoint = self.framePoints-(endPoint);
                            else
                                %iqdata = iqdata(1:max(size(iqdata))-ceil(self.resample_fs*(pri-pw))-1);
                                if length(size(iqdata))==1
                                    iqdata = self.shift_with_zero(iqdata, endPoint, 0);
                                else 
                                    iqdata = self.shift_with_zero(reshape(iqdata,1,[]), endPoint, 0);
                                end
                                tmp_ = (self.framePoints-endPoint)/(pri*self.resample_fs);
                                if (tmp_-floor(tmp_))*(pri) > pw %移动后为空
                                    startPoint = self.framePoints-(endPoint+floor(tmp_)*self.resample_fs*pri+pw*self.resample_fs);
                                end
                                endPoint = self.framePoints-endPoint;
                            end
                            iqdata = iqdata(max(size(iqdata))-self.framePoints+1:max(size(iqdata)));
                        end
                    end
                    % endPoint_ = self.framePoints-startPoint;
                    % startPoint = self.framePoints-endPoint;
                    % endPoint = endPoint_;
                    if startPoint>endPoint
                        disp('error')
                    end
                elseif ismember(i,self.pulse_modulate_class2)
                    tmp_time_domain = zeros(self.framePoints,1);
                    position = randi([0,3]);
                    pri = self.rand_uniform(1e-6,50e-5);
                    centerFreq = self.rand_uniform(-self.resample_fs/2,self.resample_fs/2)*1.0;
                    amp = self.rand_uniform(-10,10);
                    
                    if strcmp(i,"2FSK")
                        %#@audit 注意带宽/码率
            
                        bandwidth = self.rand_uniform(1e7,1e8)*1.0;
                        repeat_rate = randi([1000,50000])*10;
                        symbolrate = 64e9/repeat_rate; % 可能有问题
                        tmp_exp = sprintf('repelem(((randi(2,70,1)-1)/(2-1)-0.5)*2, %f)',repeat_rate);
                        
                        [iqdata, marker, numRepeats, chMap] = iqpulse('sampleRate', 64e9, 'PRI', pri, 'PW', pri,...
                            'pulseShape', 'Raised Cosine', 'span', bandwidth, 'offset', centerFreq, 'fmFormula', tmp_exp,...
                            'amplitude', amp, 'modulationType', 'User defined', 'correction', 0, 'channelMapping', [1,0], 'normalize', 1);
                    end
                    if strcmp(i,"frequency_hop")
                        %#@audit 注意带宽/码率
                        
                        bandwidth = self.rand_uniform(1e7,5e8)*1.0;
                        repeat_rate = randi([8000,80000]);
                        symbolrate = 64e9/repeat_rate; %# 可能有问题
                        n = randi([4,9]);
                        tmp_exp = sprintf('repelem(((randi(%d,70,1)-1)/(%d-1)-0.5)*2,%f)',n,n,repeat_rate);
                        
                        [iqdata, marker, numRepeats, chMap] = iqpulse('sampleRate', 64e9, 'PRI', pri, 'PW', pri, 'pulseShape', 'Raised Cosine', 'span', bandwidth, 'offset', centerFreq, 'fmFormula', tmp_exp,'pmFormula', '360*floor(x*4)/4','amplitude', amp, 'modulationType', 'User defined', 'correction', 0, 'channelMapping', [1,0], 'normalize', 1);
                    end
                    if strcmp(i,"costas")
                        bandwidth = self.rand_uniform(1e7,5e8)*1.0;
                        n = randi([9,15]);
                        repeat_rate = randi([20000,50000]);
                        symbolrate = 64e9/repeat_rate; %# 可能有问题
                        tmp_exp = sprintf('repelem((costas_sequence(%d)/(%d)-0.5)*2,%f)', n,n,repeat_rate);
                        [iqdata, marker, numRepeats, chMap] = iqpulse('sampleRate', 64e9, 'PRI', pri, 'PW', pri,...
                            'pulseShape', 'Raised Cosine', 'span', bandwidth, 'offset', centerFreq, 'fmFormula', tmp_exp,...
                            'amplitude', 0, 'modulationType', 'User defined', 'channelMapping', [1,0], 'normalize', 1);
                    end
                    if strcmp(i,'fmcw')
                        bandwidth = self.rand_uniform(1e7,5e8)*1.0;
                        if self.resample_fs/2-abs(centerFreq)<bandwidth/2
                            centerFreq=(self.resample_fs/2-bandwidth/2)*centerFreq/abs(centerFreq);
                        end
                        pri = self.rand_uniform(1e-6,10e-6);
                        pw = self.rand_uniform(0.4,0.8)*pri;
                        rt = self.rand_uniform(1,7)*1e-9;
                        ft = rt;
                        if pw+rt+ft > pri
                            rt = (pri-pw)/2;
                            ft = rt;
                        end
                        [iqdata, marker, numRepeats, chMap] = iqpulse('sampleRate', 64e9, 'PRI', pri, 'PW', pw, 'riseTime', rt,...
                            'fallTime', ft, 'pulseShape', 'Raised Cosine', 'span', bandwidth, 'offset', centerFreq, 'amplitude', amp,...
                            'modulationType', i, 'correction', 0, 'channelMapping', [1,0], 'normalize', 1);
                        n = self.framePoints/length(iqdata);
                        if n>1
                            n = randi([1,ceil(n)+2]);
                            iqdata = repmat(iqdata,n,1);
                        end
                    end
                    if strcmp(i,"nlfm_sine")
                        bandwidth = self.rand_uniform(1e7,5e8)*1.0;
                        if 2.5e9-abs(centerFreq)<bandwidth/2
                            centerFreq=(self.resample_fs/2-bandwidth/2)*(centerFreq/abs(centerFreq));
                        end
                        f = self.rand_uniform(2,100);
                        tmp_exp = sprintf('sin(%f*pi*x)', f);
                        [iqdata, marker, numRepeats, chMap] = iqpulse('sampleRate', 64e9, 'PRI', pri, 'PW', pri, 'pulseShape', 'Raised Cosine', 'span', bandwidth, 'offset', centerFreq, 'amplitude', 0, 'fmFormula', tmp_exp, 'pmFormula', '0', 'exactPRI', 0, 'modulationType', 'User defined', 'correction', 0, 'channelMapping', [1,0], 'normalize', 1);
                    end

                    iqdata = resample(iqdata, self.resample_fs, self.fs); %# 这个是我们需要的
                    %iqdata = conj(iqdata);
                    if max(size(iqdata)) >= self.framePoints %# 如果信号长度大于一帧的长度
                        if position == 0
                            startPoint = ceil(self.rand_uniform(floor(self.framePoints*0.1),floor(self.framePoints*0.7))); % 信号截取的位置
                            tmp_time_domain(startPoint:self.framePoints,:) = iqdata(1:self.framePoints-startPoint+1,:);
                            iqdata = tmp_time_domain;
                            endPoint = self.framePoints;

                        elseif position == 1
                            startPoint = ceil(self.rand_uniform(0,size(iqdata,1)-self.framePoints)); % 信号截取的位置
                            iqdata = iqdata(startPoint:self.framePoints+startPoint-1,1:size(iqdata,2));
                            startPoint = 0;endPoint = self.framePoints;

                        else
                            endPoint = ceil(self.rand_uniform(floor(self.framePoints*0.1),floor(self.framePoints*0.7))); % 信号截取的位置
                            tmp_time_domain(1:endPoint,1:size(tmp_time_domain,2)) = iqdata(size(iqdata,1)-endPoint+1:size(iqdata,1),1:size(iqdata,2));
                            iqdata = tmp_time_domain;
                        end

                    else %# 如果信号长度小于一帧的长度
                        tmp = floor((self.framePoints-size(iqdata,1))/2);
                        startPoint = ceil(self.rand_uniform(0,tmp)); % 信号截取的位置
                        endPoint = startPoint+size(iqdata,1)-1;
                        tmp = zeros(self.framePoints,1);
                        tmp(startPoint:endPoint,:) = iqdata;
                        iqdata = tmp;
                    
                    end    
                elseif ismember(i,self.others)
                    centerFreq = self.rand_uniform(-self.resample_fs/2*0.9,self.resample_fs/2*0.9);
                    bandwidth = 4e8;%self.rand_uniform(8e7,3e8);
                    if strcmp(i,"adsb")
                        ads_b = ADSB(self.resample_fs,self.framePoints);
                        ads_b = ads_b.call();
                        iqdata = ads_b.iq_data;
                        centerFreq = 1090e6;
                        bandwidth = 1e6;
                        startPoint = ads_b.startPoint;
                        endPoint = ads_b.endPoint;
                    elseif strcmp(i,'dsb')%'dsb','fm','usb','lsb' 
                        if self.resample_fs/2-abs(centerFreq)<bandwidth/2
                            centerFreq=(self.resample_fs/2-bandwidth/2)*centerFreq/abs(centerFreq);
                        end
                        dsb = DSB(self.resample_fs,centerFreq,self.framePoints,bandwidth);
                        dsb = dsb.call();
                        iqdata = dsb.iq_data;
                        startPoint = dsb.startPoint;
                        endPoint = dsb.endPoint;
                        
                    elseif strcmp(i,'fm')
                        if self.resample_fs/2-abs(centerFreq)<bandwidth/2
                            centerFreq=(self.resample_fs/2-bandwidth/2)*centerFreq/abs(centerFreq);
                        end
                        fm = FM(self.resample_fs,centerFreq,self.framePoints,bandwidth);
                        fm = fm.call();
                        iqdata = fm.iq_data;
                        startPoint = fm.startPoint;
                        endPoint = fm.endPoint;
                    elseif strcmp(i,'usb')
                        if self.resample_fs/2-abs(centerFreq)<bandwidth/2
                            centerFreq=(self.resample_fs/2-bandwidth/2)*centerFreq/abs(centerFreq);
                        end
                        usb = USB(self.resample_fs,centerFreq,self.framePoints,bandwidth);
                        usb = usb.call();
                        iqdata = usb.iq_data;
                        startPoint = usb.startPoint;
                        endPoint = usb.endPoint;
                        if centerFreq>0
                            centerFreq = centerFreq+bandwidth/2;
                        else
                            centerFreq = centerFreq-bandwidth/2;
                        end
                    elseif strcmp(i,'lsb')
                        if self.resample_fs/2-abs(centerFreq)<bandwidth/2
                            centerFreq=(self.resample_fs/2-bandwidth/2)*centerFreq/abs(centerFreq);
                        end
                        lsb = LSB(self.resample_fs,centerFreq,self.framePoints,bandwidth);
                        lsb = lsb.call();
                        iqdata = lsb.iq_data;
                        startPoint = lsb.startPoint;
                        endPoint = lsb.endPoint;
                        if centerFreq<0
                            centerFreq = centerFreq+bandwidth/2;
                        else
                            centerFreq = centerFreq-bandwidth/2;
                        end
                    end
                    
                end
                if size(time_domain_signal) == size(iqdata)
                    time_domain_signal = iqdata+time_domain_signal;
                else
                    time_domain_signal = conj(iqdata')+time_domain_signal;
                end
                % if centerFreq<0
                    % centerFreq=self.resample_fs/2+centerFreq;
                % end
                
                label_ = [startPoint/self.framePoints,endPoint/self.framePoints,bandwidth/self.resample_fs,centerFreq/self.resample_fs,find(self.all_modulation==i)];
                if length(label_)~=5
                    disp('error')
                end
                label = [label;label_];
                
            end
        end
        
%% 其他
        function call(self)
            self.check_path();
            if exist([self.path,'seq.mat'],'file')
                pp = load([self.path,'seq.mat']);
                self.seq = pp.seq_;
            else
                self.seq = 0;
            end
            h = waitbar(0,'Wait until finish.');
            for i = 1:self.number
                if i == self.number
                    waitbar(i/self.number,h,'DONE.');
                else
                    waitbar(i/self.number,h,'Wait until finish.');
                end
                for tmp_j = 1:length(self.snr_)
                    j = self.snr_(tmp_j);
                    [time_domain,label] = self.get_frame();
                    noisy_time_domain = awgn(time_domain,j);
                    win = rectwin(self.resolution(1));
                    % figure(1);
                    % stft(noisy_time_domain,self.resample_fs,'Window',win,'OverlapLength',self.overlap,...
                    %     'FFTLength',self.resolution(1),'FrequencyRange','centered'); %#@audit 时频，可以直接替换
                    % figure(2);
                    [tf,f,t] = stft(noisy_time_domain,self.resample_fs,'Window',win,'OverlapLength',self.overlap,...
                    'FFTLength',self.resolution(1),'FrequencyRange','centered');
                    heatmapshow = self.NormMinandMax(abs(tf));
                    imagesc(heatmapshow);
                    cmp=[parula,turbo,hsv,hot,cool,spring,summer,autumn,winter,gray,bone,copper,pink,sky ,jet];
                    sz_cmp=size(cmp,2)/3;cz=randi([0,sz_cmp-1]);
                    colormap(cmp(:,3*cz+1:3*(cz+1)));
                    axis off; % 关闭坐标轴
                    colorbar off; % 关闭颜色条
                    % 捕获图形的内容
                    frame = getframe(gca);
                    im = frame2im(frame);
                    % 保存图像
                    resizedImage = imresize(im,[self.resolution(1),self.resolution(2)]);
                    imwrite(resizedImage, char([self.path,'images/',char(string(self.seq)), '_', char(num2str(j)),'.png']));
                    for index=1:size(label,1)
                        sub=label(index,:);
                        time_duration = sub(2)-sub(1);
                        freq_duration = sub(3);
                        cent_fre = sub(4);
                        % if cent_fre>0.5
                        %     y_start = 1-cent_fre;
                        % else
                            y_start = cent_fre+0.5;
                        % end
                        if time_duration<0 || freq_duration<0
                            disp('error');
                        end
                        if self.check
                            rectangle('position',[ceil(sub(1)*self.resolution(1)) ceil((((y_start)-0.5*freq_duration)*self.resolution(2)))-3 ...
                                ceil(time_duration*self.resolution(1)) ceil(freq_duration*self.resolution(2))+6], 'EdgeColor', 'r');
                            text(ceil(sub(1)*self.resolution(1)),ceil((((y_start)-0.5*freq_duration)*self.resolution(2)))-3,self.all_modulation(sub(5)));
                        end

                   
                        label(index,1) = sub(1)+0.5*time_duration;
                        label(index,2) = y_start;
                        label(index,3) = time_duration;
                        label(index,4) = freq_duration+6/self.resolution(2);

                    end
                    if self.check
                        frame = getframe(gca);
                        im = frame2im(frame);
                        resizedImage = imresize(im,[self.resolution(1),self.resolution(2)]);
                        imwrite(resizedImage, [self.path,'check/',char(string(self.seq)), '_', char(num2str(j)),'.png']);
                    end
                    fid = fopen([self.path, 'labels/', char(num2str(self.seq)), '_', char(num2str(j)), '.txt'], 'w');
                    tf = tf';
                    save([self.path,'raw_complex/',char(string(self.seq)), '_', char(num2str(j))],'tf');
                    for n = 1:size(label, 1)
                        fprintf(fid, '%f %f %f %f %f %f\n', label(n,1), label(n, 2:5), j); % class, xywh
                    end
                    fclose(fid);
                    close all;
                    self.seq = self.seq+1;seq_=self.seq;
                    save([self.path,'seq.mat'],'seq_');
                end
            end
        end
 %% 辅助函数
        function last = NormMinandMax(self,arr)
            min_ = 0;
            max_ = 255;
            Ymax = max(arr,[],'all');  % 计算最大值
            Ymin = min(arr,[],'all');  % 计算最小值
            k = (max_ - min_) / (Ymax - Ymin);
            last = min_ + k * (arr - Ymin);
            
        end

    end
end