function varargout = MMAT(varargin)
% MMAT M-file for MMAT.fig
%      MMAT, by itself, creates a new MMAT or raises the existing
%      singleton*.
%
%      H = MMAT returns the handle to a new MMAT or the handle to
%      the existing singleton*.
%
%      MMAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MMAT.M with the given input arguments.
%
%      MMAT('Property','Value',...) creates a new MMAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MMAT_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MMAT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MMAT

% Begin initialization code - DO NOT EDIT
%global variables to allow sharing between functions
global curr_path data frames;


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MMAT_OpeningFcn, ...
                   'gui_OutputFcn',  @MMAT_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MMAT is made visible.
function MMAT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MMAT (see VARARGIN)

% Choose default command line output for MMAT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MMAT wait for user response (see UIRESUME)
% uiwait(handles.application);


% --- Outputs from this function are returned to the command line.
function varargout = MMAT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectfolder.
function selectfolder_Callback(hObject, eventdata, handles)
% hObject    handle to selectfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global curr_path;

curr_path = uigetdir(pwd, 'Choose directory containing images to analyze');
set(handles.beginanalysis,'Enable','on');


function results_Callback(hObject, eventdata, handles)
% hObject    handle to results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global curr_path data frames;
figure()
%this plots the extracted features values
plot(data(:,1),data(:,2),'*');



% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global curr_path data frames currentfigure;
%changes the displayed image on the GUI to previous if it exists
if(currentfigure>1)
    currentfigure = currentfigure-1;
   imgholder=imshow(frames(currentfigure).image);
end


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global curr_path data frames currentfigure;
%advances the displayed image on the GUI to next if it exists
if(currentfigure<length(frames))
    currentfigure = currentfigure+1;
   imgholder=imshow(frames(currentfigure).image);
end


% --- Executes on button press in openfigure.
function openfigure_Callback(hObject, eventdata, handles)
% hObject    handle to openfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global curr_path data frames currentfigure;
%opens the original saved figure file of currently displayed image
hgload(['figures/' frames(currentfigure).title]);


% --- Executes on button press in beginanalysis.
function beginanalysis_Callback(hObject, eventdata, handles)
% hObject    handle to beginanalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global curr_path data frames currentfigure;

%disable all buttons on GUI during this process
set(handles.beginanalysis,'Enable','off');
set(handles.selectfolder,'Enable','off');
%disable all buttons on GUI during this process

%read the images in the choosed folder
cd(curr_path);
files_in_path=dir(curr_path)
images_list=[];
for (i=1:length(files_in_path))
    if (files_in_path(i).isdir==0)
        this_file = files_in_path(i).name;
			[pathname,filename,ext,ver] = fileparts(this_file);
            if(ismember(upper(ext), '.TIF'))			
            images_list=strvcat(images_list, this_file);            
            end
    end
end



%initialize the frames verctor, where snapshoot images will be saved
frames=[]; 
%read in the first 2 images in tiff file from choosen file
ch1 = imread(images_list(1,:));
ch2 = imread(images_list(2,:));
%creating separate folders to hold images, figures, and data
mkdir('images');
mkdir('figures');
mkdir('m-files');
%create a temp rgb image to help create colored versions of images
temp=ch1;
%this creates a red image of the first channel
temp(:,:,2:3)=0;
redlayer=temp;
%this creates a red * green image of the both channels
temp(:,:,2)=ch2;
redgreenlayer=temp;
%finally creating a green image of the second channel
temp(:,:,1)=0;
greenlayer=temp;

%To improve on gridding accuracy we combine the intensities of both
%channels , this is only to help finding dividing lines
v=ch1+ch2;

%we start by global gridding in the horizontal direction to find the lines
%that divide the wells. first we create a profile of intensities
XProfile = mean(v);
%this find the period of occurance of spots in the horizontal direction 
ac = xcov(XProfile);                        %unbiased autocorrelation
s1 = diff(ac([1 1:end]));                   %left slopes
s2 = diff(ac([1:end end]));                 %right slopes
maxima = find(s1>0 & s2<0);                 %peaks
estPeriod = round(median(diff(maxima))) ;
seLine = strel('line',estPeriod,0);
%this is to create a better profile of intensities by removing the noise
%morphologically using imtophat
XProfile2 = imtophat(XProfile,seLine);
level = graythresh(XProfile2/255)*255;
bw = im2bw(XProfile2/255,level/255);
L = bwlabel(bw);
Xraisingedges=[];
Xfallingedges=[];

%removing in singular noise value
for i=2:length(L)
   if (L(i-1)<L(i))
      Xraisingedges =[Xraisingedges i];
      else if (L(i-1)>L(i))
              Xfallingedges =[Xfallingedges i];
          end
    end
 
end

%populating the vertical lines, image edges included
XGrid=1;
for r=11:size(Xraisingedges,2)-11 %ignoring first and last 10
    if(Xraisingedges(r)-Xfallingedges(r-1)>round(estPeriod*3/4))
    XGrid = [XGrid round((Xraisingedges(r)+Xfallingedges(r-1))/2)];
    end
end
XGrid=[XGrid size(v,2)];

% now we average subgrids' widths, if we find subgrid width below the
% average, we remove it because it's most likely a mistake and repeat again
badpoints = [];
keeplooping =true;
while(keeplooping)
gaps = diff(XGrid);

for i=1:length(gaps)-1
    if (gaps(i)<median(gaps)/2)       
        badpoints=[badpoints i];
    end
end

if(~(isempty(badpoints)))
XGrid(badpoints)=[];
badpoints = [];
else
    keeplooping =false;
end

end

%This will repeat the same global gridding process on the vertical spacing,
%we will simply use the inverse of the image vector and proceed as before.
YProfile = mean(v');
ac = xcov(YProfile);                        %unbiased autocorrelation
s1 = diff(ac([1 1:end]));                   %left slopes
s2 = diff(ac([1:end end]));                 %right slopes
maxima = find(s1>0 & s2<0);                 %peaks
estPeriod = round(median(diff(maxima))); 
seLine = strel('line',estPeriod,0);
YProfile2 = imtophat(YProfile,seLine);
level = graythresh(XProfile2/255)*255;
bw = im2bw(YProfile2/255,level/255); 
L = bwlabel(bw);
Yraisingedges=[];
Yfallingedges=[];

%removing in singular noise value
for i=2:length(L)
  if (L(i-1)<L(i))
      Yraisingedges =[Yraisingedges i];
      else if (L(i-1)>L(i))
              Yfallingedges =[Yfallingedges i];
          end
      end
  
end
YGrid=1;

for r=11:size(Yraisingedges,2)-11 %ignoring first and last 10
    if(Yraisingedges(r)-Yfallingedges(r-1)>=round(estPeriod*3/4))
    YGrid = [YGrid round((Yraisingedges(r)+Yfallingedges(r-1))/2)];
    end
end
YGrid=[YGrid size(v,1)];

badpoints = [];
keeplooping =true;
while(keeplooping)
gaps = diff(YGrid);
for i=1:length(gaps)-1
    if (gaps(i)<median(gaps)/2)
        badpoints=[badpoints i];
    end
end

if(~(isempty(badpoints)))
YGrid(badpoints)=[];
badpoints = [];
else
    keeplooping =false;
end

end


%the following will display the outcome of the global gridding as well as
%the subgrid numbers on green, red, combined layers.
figure(1);
imshow(redlayer)
line(XGrid'*[1 1],YGrid([1 end]),'color','w')
line(XGrid([1 end]),YGrid'*[1 1],'color','w')
for s=1:size(YGrid,2)-1
  for j=1:size(XGrid,2)-1
text((XGrid(j)+XGrid(j+1))/2,(YGrid(s)+YGrid(s+1))/2,sprintf(...
    'Well %dx%d',s,j),'color','w','HorizontalAlignment','center')
  end 
end
title('Red channel wells outline');
%we capture a snapshot of the current figure
f1=getframe(figure(1));
frames(length(frames)+1).image=f1.cdata;
frames(length(frames)).title='Red-wells-outline';
%next we save the snapshot image as well as the original figure
imwrite(f1.cdata,'images/Red-wells-outline.jpg','jpg','Quality',100);
hgsave('figures/Red-wells-outline');

figure(1);
imshow(greenlayer)
line(XGrid'*[1 1],YGrid([1 end]),'color','w')
line(XGrid([1 end]),YGrid'*[1 1],'color','w')
for s=1:size(YGrid,2)-1
  for j=1:size(XGrid,2)-1
text((XGrid(j)+XGrid(j+1))/2,(YGrid(s)+YGrid(s+1))/2,sprintf(...
    'Well %dx%d',s,j),'color','w','HorizontalAlignment','center')
  end 
end
title('Green channel wells outline');
f1=getframe(figure(1));
frames(length(frames)+1).image=f1.cdata;
frames(length(frames)).title='Green-wells-outline';
imwrite(f1.cdata,'images/Green-wells-outline.jpg','jpg','Quality',100);
hgsave('figures/Green-wells-outline');

figure(1);
imshow(redgreenlayer)
line(XGrid'*[1 1],YGrid([1 end]),'color','w')
line(XGrid([1 end]),YGrid'*[1 1],'color','w')
for s=1:size(YGrid,2)-1
  for j=1:size(XGrid,2)-1
text((XGrid(j)+XGrid(j+1))/2,(YGrid(s)+YGrid(s+1))/2,sprintf(...
    'Well %dx%d',s,j),'color','w','HorizontalAlignment','center')
  end 
end
title('Combined channels wells outline');
f1=getframe(figure(1));
frames(length(frames)+1).image=f1.cdata;
frames(length(frames)).title='Combined-wells-outline';
imwrite(f1.cdata,'images/Combined-wells-outline.jpg','jpg','Quality',100);
hgsave('figures/Combined-wells-outline');
close(1);

intensities=[];
intensities2=[];

%the next process is a local gridding on each subgrid found from the
%previous process by cropping the combined channel to subgrid dimensions
for s=1:length(YGrid)-1
  for j=1:length(XGrid)-1
%crop combined channel only to find dividing lines
z = imcrop(v,[XGrid(j) YGrid(s) XGrid(j+1)-XGrid(j) YGrid(s+1)-YGrid(s)]);
%crops the red channel for feature extraction after gridding is done
z1 = imcrop(ch1,[XGrid(j) YGrid(s) XGrid(j+1)-XGrid(j) YGrid(s+1)-YGrid(s)]);
%crops the green channel for feature extraction after gridding is done
z2 = imcrop(ch2,[XGrid(j) YGrid(s) XGrid(j+1)-XGrid(j) YGrid(s+1)-YGrid(s)]);
%crops the red & green channels for illustration purposes only
zboth = imcrop(redgreenlayer,[XGrid(j) YGrid(s) XGrid(j+1)-XGrid(j) YGrid(s+1)-YGrid(s)]);

%similar process to the previous gridding but on a subgrid level
xProfile = mean(z);
ac = xcov(xProfile);                        %unbiased autocorrelation
s1 = diff(ac([1 1:end]));                   %left slopes
s2 = diff(ac([1:end end]));                 %right slopes
maxima = find(s1>0 & s2<0);                 %peaks
estPeriod = round(median(diff(maxima))) ;
seLine = strel('line',estPeriod,0);
xProfile2 = imtophat(xProfile,seLine);
level = graythresh(xProfile2/255)*255;
bw = im2bw(xProfile2/255,level/255);
L = bwlabel(bw);
xraisingedges=[];
xfallingedges=[];
%removing in singular noise value
for i=2:length(L)
   if (L(i-1)<L(i))
      xraisingedges =[xraisingedges i];
      else if (L(i-1)>L(i))
              xfallingedges =[xfallingedges i];
          end
    end
 
end

badpoints = [];
for i=1:min([length(xraisingedges) length(xfallingedges)])
   if (xfallingedges(i)<xraisingedges(i))
       xfallingedges(i)=[];
   else   if (xfallingedges(i)-xraisingedges(i)<10)
    badpoints=[badpoints i];
       end
   end 
end
 xraisingedges(badpoints)=[];
 xfallingedges(badpoints)=[];

 %adding padding to the edges of the grid
xGrid=xraisingedges(1)-4;
for r=2:size(xraisingedges,2) %ignoring first and last 10
    
    xGrid = [xGrid round((xraisingedges(r)+xfallingedges(r-1))/2)];
end


yProfile = mean(z');
ac = xcov(yProfile);                        %unbiased autocorrelation
s1 = diff(ac([1 1:end]));                   %left slopes
s2 = diff(ac([1:end end]));                 %right slopes
maxima = find(s1>0 & s2<0);                 %peaks
estPeriod = round(median(diff(maxima))); 
seLine = strel('line',estPeriod,0);
yProfile2 = imtophat(yProfile,seLine);
level = graythresh(xProfile2/255)*255;
bw = im2bw(yProfile2/255,level/255); 
L = bwlabel(bw);
yraisingedges=[];
yfallingedges=[];
%removing in singular noise value

for i=2:length(L)
  if (L(i-1)<L(i))
      yraisingedges =[yraisingedges i];
      else if (L(i-1)>L(i))
              yfallingedges =[yfallingedges i];
          end
      end
  
end


badpoints = [];
for i=1:min([length(yraisingedges) length(yfallingedges)])
   if (yfallingedges(i)<yraisingedges(i))
       yfallingedges(i)=[];
   else   if (yfallingedges(i)-yraisingedges(i)<10)
    badpoints=[badpoints i];
       end
   end 
end
 yraisingedges(badpoints)=[];
 yfallingedges(badpoints)=[];
 
%adding padding to the edges of the grid
yGrid=yraisingedges(1)-4;
for r=2:size(yraisingedges,2) %ignoring first and last 10
   
    yGrid = [yGrid round((yraisingedges(r)+yfallingedges(r-1))/2)];
end


if(length(yGrid)>length(xGrid))
    yGrid(end-(length(yGrid)-length(xGrid))+1:end)=[];
else if (length(xGrid)>length(yGrid))
         xGrid(end-(length(xGrid)-length(yGrid))+1:end)=[];
    end
end

%adding padding to the edges of the grid
xGrid=[xGrid xfallingedges(end)+4];
yGrid=[yGrid yfallingedges(end)+4];
   
%showing the result of gridding on
figure(2);
imshow(zboth);
line(xGrid'*[1 1],yGrid([1 end]),'color','b')
line(xGrid([1 end]),yGrid'*[1 1],'color','b')
[X,Y] = meshgrid(xGrid(1:end-1),yGrid(1:end-1));
[dX,dY] = meshgrid(diff(xGrid),diff(yGrid));
ROI = [X(:) Y(:) dX(:) dY(:)];

%segmentation process:
%the following line of code mimics the im2bw function but with a very low
%threshild to recognize even the most faded spots, which essentially
%creates logical map of which pixels have intensity informaion of interest
BW2 = z>=4000;
%we then fill those shapes
BW_filled = imfill(BW2,'holes');
%we then filter the noise produced from over segmenting with low threshild
bw = medfilt2(BW_filled,[10 10]);
%tracing the boundaries of the remaning segmented shapes
boundaries = bwboundaries(bw);
for k=1:length(boundaries)
   b = boundaries{k};
   if(length(b)>2)
   hold on;
   plot(b(:,2),b(:,1),'b','LineWidth',1);
   end
end

%display and save the result of local gridding as both a jpg & figure format
title(sprintf('Well %dx%d',s,j));
f1=getframe(figure(2));
frames(length(frames)+1).image=f1.cdata;
frames(length(frames)).title=sprintf('Well%dx%d',s,j);
imwrite(f1.cdata,sprintf('images/Well%dx%d.jpg',s,j),'jpg','Quality',100);
hgsave(sprintf('figures/Well%dx%d',s,j));
close(2);

%creating masks from every shape found to be used for feature extraction
L = zeros(size(bw));
for i=1:length(ROI)
  rows = ROI(i,2)+[0:(ROI(i,4)-1)];
  cols = ROI(i,1)+[0:(ROI(i,3)-1)];
  rectMask = L(rows,cols);
  spotMask = bw(rows,cols);
  rectMask(spotMask) = i;
  L(rows,cols) = rectMask;
end

%loop through every spot in the subgrid, extracting the feature value of
%each channel by virtue of the previously created mask if it exists
for g=1:length(ROI)
 rect = ROI(g,:);                            %[X Y dX dY]
spot1 = imcrop(z1,rect);                      %region around spot
spot2 = imcrop(z2,rect); 
mask = imcrop(L,rect)==g;
if(~isempty(find(mask)))
 intensities = [intensities double(median(spot1(mask)))];
 intensities2 = [intensities2 double(median(spot2(mask)))];
else
intensities = [intensities 0];
 intensities2 = [intensities2 0];
end
end


  end
end

%saving the vectors of each channel intensities as well as their difference
data = [intensities' intensities2' intensities2'-intensities'];
save('m-files/analysis-data.mat', 'data');
save('m-files/analysis-images.mat', 'frames');

%process is done, we enable all buttons on the GUI
set(handles.results,'Enable','on');
set(handles.selectfolder,'Enable','on');
set(handles.beginanalysis,'Enable','off');
set(handles.previous,'Visible','on');
set(handles.next,'Visible','on');
set(handles.openfigure,'Visible','on');
set(handles.imgholder,'Visible','on');
%display the first saved image of red channel global gridding
currentfigure=1;
imgholder=imshow(frames(currentfigure).image);



