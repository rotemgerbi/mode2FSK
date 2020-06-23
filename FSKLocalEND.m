%%Reset
close all;
clear all;
clc;

%%
%"Transmitter"
fileID = fopen('test.txt');
FileBits=fread(fileID,'ubit1');
currentCharacterEncoding = slCharacterEncoding()
% fof = fopen('binfile.bin','w');
% fwrite(fof,FileBits);
% fclose(fof);
% FSK Modulationmyfilebin and Demodulation in AWGN
% Modulate and demodulate a signal using 8-FSK modulation with a frequency
% separation of 25kHz.


% Set the modulation order and frequency separation parameters.
M =2;%Number of frequencies
freqSep = 25e3; %Distance between frequencies

% Create FSK modulator and demodulator System objects(TM) with modulation
% order 2 and 5kHz frequency separation.

SamplesPerSymbol=10;%samples per Symbol*
BaudRate=20e3;%symbols/Second
Fsampling=SamplesPerSymbol*BaudRate;
%Modulation
hMod = comm.FSKModulator (M,freqSep,'SamplesPerSymbol',SamplesPerSymbol,'SymbolRate',BaudRate);

%%
%Setting Up variables and matrices
%build preamble
PreambleSize=100;
preamble=zeros(PreambleSize,1);
preamble(1:20)=1;
preamble(21:80)=0;
preamble(81:100)=1;
%  preamble=randi([0,1],1,100).';
%Setting Up variables and matrices
SerialNumberSize=10; %Serial Number
SerialNumber=zeros(SerialNumberSize,1);
BITSfileSize=20; % File Size`
BITSNumberfileSize=zeros(BITSfileSize,1) ;
PayLoadSize=200; %bits

%%
TotalFrames=ceil(size(FileBits,1)/PayLoadSize);
packetsize=BITSfileSize+PayLoadSize+SerialNumberSize+PreambleSize;
packetSizeSamples=packetsize*SamplesPerSymbol;
%%
%Create three main matrices:1.modSignal-Information with extra noise.
%2.packets-Information is divided into messages.
%3.Quiet_breaks-Information with silent times.
modSignal=zeros(packetSizeSamples,TotalFrames);
packets= zeros(PreambleSize+BITSfileSize+SerialNumberSize+PayLoadSize,TotalFrames);

%This loop fills the packets matrix with preamble and PayLoad.
%Until the last line, where added the missing and adds zeros.
%Packets is=[preamble SerialNumberSize BITSfileSize PayLoadSize].
for counter = 1:TotalFrames
    
    if counter~=TotalFrames
        packets(:,counter)= [preamble;...
            de2bi(counter,SerialNumberSize,'left-msb').';....
            de2bi(length(FileBits),BITSfileSize,'left-msb').';...
            FileBits((counter-1)*PayLoadSize+1:(counter-1)*PayLoadSize+PayLoadSize)];
        
    else
        len_last=length(FileBits((counter-1)*PayLoadSize+1:end));
        packets (:,counter) =[preamble;...
            de2bi(counter,SerialNumberSize,'left-msb').';...
            de2bi(length(FileBits),BITSfileSize,'left-msb').';...
            FileBits((counter-1)*PayLoadSize+1:end); zeros(PayLoadSize-len_last,1)];
        
    end
    
    modSignal(:,counter)=step(hMod,packets(:,counter));
    
end

%Insert random zeros at the end and at the beginning.
%To understand that the receiver does not know when they will start broadcasting
%and when they will stop.
modSignal=modSignal(:);
% [bw1,flo1,fhi1,power1] = obw(modSignal);

%%
% Test SNR before attenuation and noise
% SNR_Before=snr(real(modSignal));
% [bw,flo,fhi,power] = obw(modSignal);
%%
%%The channel

modSignalNoise_Atten=[ zeros(packetSizeSamples*TotalFrames,1); modSignal ;zeros(packetSizeSamples*TotalFrames,1) ];

%Call to the Attenuate_dB function
% attenFactordB=3;
% modSignal_attenuated=AttenuatedB(modSignal,attenFactordB);
% 
% %In order to create a graph
% % BER_arr=[];
% % SNR_Arr=(1:1:50);
% %
% %Four times to make an average
% % for z=1:4
% %     for i=1:length(SNR_Arr)
% %
% % Add noise
% SNR=1000;
% modSignalNoise_Atten=awgn(modSignal_attenuated,SNR);
% NoiseVector=BuildNoiseBySNR(modSignalNoise_Atten,SNR);
% % [bw1,flo1,fhi1,power1] = obw(NoiseVector);
% Spectrum = dsp.SpectrumAnalyzer('SampleRate',Fsampling);
% while 1
%     step(Spectrum, modSignalNoise_Atten);
% end
%
%%
%Receiver
%Test SNR after attenuation and noise
% SNR_After=snr(real(modSignalNoise_Attent))
% figure;
% stem(abs(modSignalNoise_Atten));

%Object deciphering
hDeMod = comm.FSKDemodulator (M,freqSep,'SamplesPerSymbol',SamplesPerSymbol,'SymbolRate',BaudRate);

%%Setting variables
% PreambleSize=100;
% Preamble=zeros(SPreambleSize,1);
% Preamble(1:20)=1;
% Preamble(21:80)=0;
% Preamble(81:100)=1;
PayloadeReceiver=200;
SerialSize=10;
BITSSize=20;

%First bit modulation, for correlation
release(hMod);
SPreamble_hMod=step(hMod,preamble);

length_modSignalNoise=length(modSignalNoise_Atten);
zerost=length_modSignalNoise-length(SPreamble_hMod);
Message_size=PayloadeReceiver+SerialSize+BITSSize+PreambleSize;
packetSizeSamples=SamplesPerSymbol*Message_size;
zerost_matrix=zeros(1,zerost);

%In order to normalize should SPreamble_hMod and modSignalNoise_Attent be
%equal.
SPreamble_hMod_zeros=[SPreamble_hMod.' zerost_matrix];

%Correlation between modSignalNoise_Attent and SPreamble_hMod_zeros+plot.
% [acor lag] = xcorrEND(modSignalNoise_Atten,SPreamble_hMod_zeros,'coeff');
% [acor lag] = xcorr(modSignalNoise_Atten,SPreamble_hMod);
 [acor lag] = xcorr(modSignalNoise_Atten,SPreamble_hMod_zeros,'coeff');

% %Plot correlation

% figure;
% acorAbs=abs(acor);
% figure;
% plot(acorAbs);

% spectrum = dsp.SpectrumAnalyzer('SampleRate', Fsampling);
% % info(spectrum)
% % % % show the real frequency and not 0
% while 1
% step(spectrum,acor)
% end
% [Msg_deModp , suspectST ]=CorrelationP(acor,lag,modSignalNoise_Atten, M,freqSep,BaudRate,SamplesPerSymbol,preamble);
% % Discovery of the first maximum
% % [maxVal maxIdx]=max(acorAbs);
% % suspectST=lag(find(acorAbs==max(acorAbs)));
[~,I] = max(abs(acor));
suspectST = lag(I);
suspectST=suspectST(1);
suspectPacket=modSignalNoise_Atten(suspectST:suspectST+packetSizeSamples-1);
% % 
% % Decrypt the entire message

Msg_deMod=step(hDeMod,suspectPacket);

stem(abs(Msg_deMod));
% % Xor_Preamble=sum(xor(Msg_deMod(1:100),preamble));
% % if Xor_Preamble~=0
% %     XOR_PREAMBLE(acorAbs);
% % end

%Discovery file size+packet one
FileSizeIs=Msg_deMod(111:130);
FileSizeIs=bi2de(FileSizeIs.',2,'left-msb');
NumMsg_One=Msg_deMod(101:110);
NumMsg_One=bi2de(NumMsg_One.',2,'left-msb');

%Discovery number packet
TotalFrames= ceil(FileSizeIs/PayloadeReceiver);

%the number of left and right of remaining
Cut_Right=TotalFrames-NumMsg_One;
Cut_Left=TotalFrames-Cut_Right;

%%Treatment in extreme case + Expected case.
%Reveals the file size, according to the started, Of the first packet.

if NumMsg_One==1 %Input the one packet first
    Msg_Total=modSignalNoise_Atten(suspectST:suspectST+packetSizeSamples*TotalFrames-1);
end

if NumMsg_One==TotalFrames %Input the last packet first
    Msg_Total=modSignalNoise_Atten(suspectST-packetSizeSamples*TotalFrames:(suspectST+packetSizeSamples)-1);
end

if NumMsg_One~=TotalFrames &&NumMsg_One~=1 %Input the packet between the last and the first
    Msg_Total=modSignalNoise_Atten(((suspectST-(packetSizeSamples*(Cut_Left-1)))):(suspectST+(packetSizeSamples*(Cut_Right+1))-1));
end

%Decrypt the message
release(hDeMod);
Msg_Total_deMod=step(hDeMod,Msg_Total);


%build resullt file
for i=1:Message_size:length(Msg_Total_deMod)
    
    NumberMsg=Msg_Total_deMod(i+100:i+109);
    NumMsg_Received=bi2de(NumberMsg.',2,'left-msb');
    result((NumMsg_Received-1)*200+1:NumMsg_Received*200)=Msg_Total_deMod(i+130:i+329);
    
end

%%Download zeros from the end file.
TotalResWithZeros=(length(Msg_Total_deMod)/330);
TotalReciver=FileSizeIs/PayloadeReceiver;
Detection_Information=PayloadeReceiver*((TotalResWithZeros-TotalReciver));
result=result(:);

if Detection_Information>0
    result=result(1:end-Detection_Information);
end

%         %If the length is bigger or smaller,subtract
%         Difference=length(result)-length(FileBits);
%         Biterros=sum(xor(FileBits, result(1:end-Difference)));
%         BER =(Biterros/length(FileBits))*100;
%         BER_arr=[BER_arr BER];

% %%Writing to a file
fileIDNew = fopen('myNewPoem.txt','w');
fwrite(fileIDNew,result,'ubit1');
fclose(fileIDNew);

% fid = fopen('myNewPoem.txt');
% if fid < 0, error('Cannot open file'); end
% C = textscan(fid, '%s', 'Delimiter', '\n');
% fclose(fid);



%     end
% end

%Average
% BER_arr1=BER_arr(1:50);
% BER_arr2=BER_arr(51:100);
% BER_arr3=BER_arr(101:150);
% BER_arr4=BER_arr(151:200);
% BER_arr_Total=BER_arr1+BER_arr2+BER_arr3+BER_arr4;
% BER_arr_Total=BER_arr_Total/4;
%
% %If no error is placed 0.00001,doing log10 on BER average,and brings out.
% for i=1:50
%     if BER_arr_Total(i)==0
%         BER_arr_Total(i)=0.00001;
%     end
%     BER_arr_Total(i)=log10(BER_arr_Total(i));
% end
%
% %Plot SNR_Arr VS BER_arr_Total(10^)
% figure;
% plot(SNR_Arr,BER_arr_Total);
% grid on
% title('SNR logarithmic function of BER');
% xlabel('SNR');
% ylabel('BER(%)__10^');