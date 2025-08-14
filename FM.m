classdef FM
    properties
        
        carrier_freq
        sample_rate
        duration
        num_iq_samples
        iq_data
        endPoint
        startPoint
        frame
        bandwidth
        beta
        rand_
    end
    methods
        %% 构造函数
        function self = FM(sample_rate,carrier_freq,num_iq_samples,bandwidth,beta)
            self.sample_rate = sample_rate;
            self.carrier_freq = carrier_freq;

            self.num_iq_samples = num_iq_samples;
            self.duration = num_iq_samples/self.sample_rate;
            self.bandwidth = bandwidth;
            self.beta = self.bandwidth/1;
            self.rand_ = true;
        end
        %% 生成随机数
        function rand_u = rand_uniform(self,low,high)
            rand_u = low + (high-low)*rand;
        end
        %% 生成随机数
        function Hd = getFilter(self)
            
            Fstop1 = abs(self.carrier_freq)-self.bandwidth/2*1.2;    % First Stopband Frequency
            if Fstop1<0
                Fstop1 = 0;
            end
            Fpass1 = abs(self.carrier_freq)-self.bandwidth/2;    % First Passband Frequency
            if Fpass1<0
                Fpass1 = 0;
            end
            Fpass2 = abs(self.carrier_freq)+self.bandwidth/2;   % Second Passband Frequency
            Fstop2 = abs(self.carrier_freq)+self.bandwidth/2*1.2;   % Second Stopband Frequency
            Astop1 = 100;          % First Stopband Attenuation (dB)
            Apass  = 1;           % Passband Ripple (dB)
            Astop2 = 100;          % Second Stopband Attenuation (dB)
            Fs     = self.sample_rate;  % Sampling Frequency
            
            h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
                Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);
            
            Hd = design(h, 'cheby1', ...
                'MatchExactly', 'passband', ...
                'SystemObject', true,...
                UseLegacyBiquadFilter=true);
            %GETFILTER Returns a discrete-time filter System object.

        end
       
       %% 主函数
       function self = call(self)
            t = linspace(0,self.num_iq_samples/self.sample_rate,self.num_iq_samples);
            tmp_time_domain = zeros(self.num_iq_samples,1);
            % carrier = exp(1i*2*pi*self.carrier_freq*t);
            %message-rand
            nn = randi([4,10]);
            message_ = zeros(1,self.num_iq_samples);
            for i=1:nn
                tmp = randi([0,2]);
                if tmp==0
                    message = chirp(t, rand()*6e5, self.num_iq_samples/self.sample_rate, rand()*6e5, 'linear');
                elseif tmp==1
                    message = chirp(t, rand()*6e5, self.num_iq_samples/self.sample_rate, rand()*6e5, 'quadratic');
                elseif tmp==2
                    message = chirp(t, rand()*6e5, self.num_iq_samples/self.sample_rate, rand()*6e5, 'logarithmic');
                end
                message_ = (message+message_);
            end
            message_=message_/nn;
            % message_ = message_';

            message_ = fmmod(message_, abs(self.carrier_freq), self.sample_rate,self.bandwidth/2);
            message_=hilbert(message_);
            
            Hd = self.getFilter();
            self.frame = Hd(message_);

            if self.carrier_freq>0
                message_ = conj(message_');
            end

            self.iq_data = zeros(1,self.num_iq_samples);
            self.frame = (message_)';

            if size(self.frame)~=size(tmp_time_domain)
                self.frame = (self.frame)';
            end
            position = randi([0,2]);
            
            if max(size(self.frame)) >= self.num_iq_samples % 如果信号长度大于一帧的长度
                if position == 0
                    self.startPoint = randi([0,ceil(self.num_iq_samples*0.5)]);
                    tmp_time_domain(self.startPoint:self.num_iq_samples,:) = self.frame(self.startPoint:self.num_iq_samples,:);
                    self.frame = tmp_time_domain;
                    self.endPoint = self.num_iq_samples;
                elseif position == 1
                    % self.startPoint = randi([0,(size(self.frame,1)-self.num_iq_samples)]); % 信号截取的位置
                    % self.frame = self.frame(self.startPoint:self.num_iq_samples+self.startPoint-1,:);
                    self.startPoint = 1;self.endPoint = self.num_iq_samples;
                else
                    self.endPoint = randi([0,ceil(self.num_iq_samples*0.5)]); % 信号截取的位置
                    tmp_time_domain(1:self.endPoint,:) = self.frame(size(self.frame,1)-self.endPoint+1:size(self.frame,1),:);
                    self.frame = tmp_time_domain;
                    self.startPoint = 1;
                end

            else % 如果信号长度小于一帧的长度
                tmp = floor((self.num_iq_samples-size(self.frame,1))/2);
                self.startPoint = ceil(self.rand_uniform(0,tmp)); % 信号截取的位置
                tmp = zeros(self.num_iq_samples,1);
                if self.startPoint+size(self.frame,1) <= self.num_iq_samples
                    tmp(self.startPoint:self.startPoint+max(size(self.frame))-1,:) = self.frame;
                    self.endPoint = self.startPoint+size(self.frame,1);
                else
                    tmp(self.startPoint:self.num_iq_samples-1,1:size(tmp,2)) = iqdata(1:self.num_iq_samples-self.startPoint);
                    self.endPoint = self.num_iq_samples;
                end
                self.frame = tmp;
            end
            self.iq_data = self.frame;
       end
    end
end