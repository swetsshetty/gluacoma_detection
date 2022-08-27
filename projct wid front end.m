function varargout = projectnew(varargin)
% PROJECTNEW MATLAB code for projectnew.fig
%      PROJECTNEW, by itself, creates a new PROJECTNEW or raises the existing
%      singleton*.
%
%      H = PROJECTNEW returns the handle to a new PROJECTNEW or the handle to
%      the existing singleton*.
%
%      PROJECTNEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECTNEW.M with the given input arguments.
%
%      PROJECTNEW('Property','Value',...) creates a new PROJECTNEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before projectnew_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to projectnew_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help projectnew

% Last Modified by GUIDE v2.5 04-May-2013 01:39:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @projectnew_OpeningFcn, ...
                   'gui_OutputFcn',  @projectnew_OutputFcn, ...
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


% --- Executes just before projectnew is made visible.
function projectnew_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to projectnew (see VARARGIN)

% Choose default command line output for projectnew
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes projectnew wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = projectnew_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    str=get(handles.edit1,'string');
    F=imread(str);
    % resize the image
    F= imresize(F,[256,256]);      
    
    B = strel('disk',15);         %structural element
    he=imclose(F,B);
axes(handles.axes2)
imshow(he);
    
  
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    str=get(handles.edit1,'string');
    F=imread(str);
    % resize the image
    F= imresize(F,[256,256]);      
    
    B = strel('disk',15);         %structural element
    he=imclose(F,B);
    cform = makecform('srgb2lab');
    % Applying above color transform to the sRGB image 
    lab_he = applycform(he,cform);
    % Converting into double 
    ab = double(lab_he(:,:,2:3));
    % obtaining rows and columns of transformed image 
    nrows = size(ab,1);
    ncols = size(ab,2);
    % Reshaping image taking each value column wise 
    ab = reshape(ab,nrows*ncols,2);
    % No of clusters to be created with five iterations 
    nColors = 5;
    [cluster_idx cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean',...
                                          'Replicates',5);
    % Reshaping and showing the clusters 
    pixel_labels = reshape(cluster_idx,nrows,ncols);
    %subplot(2,6,3), imshow(pixel_labels,[]),
    title('image labeled by cluster index');
    % creating five element array 
    segmented_images = cell(1,5);
    % Creating tiles for three different colors 
    rgb_label = repmat(pixel_labels,[1 1 3]); 
    % Assigning clustered objects to array(segmented_image) 
    for k = 1:nColors
        color = he;
        color(rgb_label ~= k) = 0;
        segmented_images{k} = color;
    end
    % displaying different cluster objects 
    %subplot(2,6,4),imshow(segmented_images{1}), title('cluster 1');
    %subplot(2,6,5),imshow(segmented_images{2}), title('cluster 2');
    %subplot(2,6,6),imshow(segmented_images{3}), title('cluster 3');
    %subplot(2,6,7),imshow(segmented_images{4}), title('cluster 4');
    %subplot(2,6,8),imshow(segmented_images{5}), title('cluster 5'); 
    %automatic selection of cluster containing cup by conveting into lab
     lab_he1 = applycform(segmented_images{1},cform);
     lab_he2= applycform(segmented_images{2},cform);
     lab_he3= applycform(segmented_images{3},cform);
     lab_he4= applycform(segmented_images{4},cform);
     lab_he5= applycform(segmented_images{5},cform);
    % calculating max point of lab of segmented image   
     m1=max(lab_he1(:));
     m2=max(lab_he2(:));
     m3=max(lab_he3(:));
     m4=max(lab_he4(:));
     m5=max(lab_he5(:));
    %  image with greatest max point taken as cup
     a=[m1 m2 m3 m4 m5];
     temp=a(1);
     for pos=1:5
         if a(pos)>=temp
             temp=a(pos);
             j=pos;
         end
     end
    % convert segmented image into grey scale   
    cup = rgb2gray(segmented_images{j});
    % convert into black and white
    cup = im2bw(cup,0.5);
    %subplot(2,6,9),imshow(cup), title('cup ');
    %remove unwanted portion of segmented cup
    SE = strel('square',3);
    cup = imerode(cup,SE);
    %find edge of cup     
    cup = double(cup) + 1;
    cup = edge(cup,'canny');
    %     [B,cup]=bwboundaries(cup,'noholes');
    %     cup=bwmorph(cup,'bridge');
    cup =bwmorph(cup,'dilate');
%aa=imread(cup);
    axes(handles.axes3)
imshow(cup);
    %figure,subplot(2,6,10),imshow(cup), title(' cup boundary');
  
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    str=get(handles.edit1,'string');
    F=imread(str);
    % resize the image
    F= imresize(F,[256,256]);      
    
    B = strel('disk',15);         %structural element
    he=imclose(F,B);
axes(handles.axes2)
imshow(he);
    
  
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    str=get(handles.edit1,'string');
    F=imread(str);
    % resize the image
    F= imresize(F,[256,256]);      
    
    B = strel('disk',15);
    
%structural element
     F=imopen(F,B);

    he=imclose(F,B);
    cform = makecform('srgb2lab');
    % Applying above color transform to the sRGB image 
    lab_he = applycform(he,cform);
    % Converting into double 
    ab = double(lab_he(:,:,2:3));
    % obtaining rows and columns of transformed image 
    nrows = size(ab,1);
    ncols = size(ab,2);
    % Reshaping image taking each value column wise 
    ab = reshape(ab,nrows*ncols,2)
    % No of clusters to be created with five iterations 
    nColors = 5;
    [cluster_idx cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean',...
                                          'Replicates',5);
    % Reshaping and showing the clusters 
    pixel_labels = reshape(cluster_idx,nrows,ncols);
    %subplot(2,6,3), imshow(pixel_labels,[]),
    title('image labeled by cluster index');
    % creating five element array 
    segmented_images = cell(1,5);
    % Creating tiles for three different colors 
    rgb_label = repmat(pixel_labels,[1 1 3]); 
    % Assigning clustered objects to array(segmented_image) 
    for k = 1:nColors
        color = he;
        color(rgb_label ~= k) = 0;
        segmented_images{k} = color;
    end
    % displaying different cluster objects 
    %subplot(2,6,4),imshow(segmented_images{1}), title('cluster 1');
    %subplot(2,6,5),imshow(segmented_images{2}), title('cluster 2');
    %subplot(2,6,6),imshow(segmented_images{3}), title('cluster 3');
    %subplot(2,6,7),imshow(segmented_images{4}), title('cluster 4');
    %subplot(2,6,8),imshow(segmented_images{5}), title('cluster 5'); 
    %automatic selection of cluster containing cup by conveting into lab
     lab_he1 = applycform(segmented_images{1},cform);
     lab_he2= applycform(segmented_images{2},cform);
     lab_he3= applycform(segmented_images{3},cform);
     lab_he4= applycform(segmented_images{4},cform);
     lab_he5= applycform(segmented_images{5},cform);
    % calculating max point of lab of segmented image   
     m1=max(lab_he1(:));
     m2=max(lab_he2(:));
     m3=max(lab_he3(:));
     m4=max(lab_he4(:));
     m5=max(lab_he5(:));
    %  image with greatest max point taken as cup
     a=[m1 m2 m3 m4 m5];
     temp=a(1);
     for pos=1:5
         if a(pos)>=temp
             temp=a(pos);
             j=pos;
         end
     end
    % convert segmented image into grey scale   
    cup = rgb2gray(segmented_images{j});
    % convert into black and white
    cup = im2bw(cup,0.5);
    %subplot(2,6,9),imshow(cup), title('cup ');
    %remove unwanted portion of segmented cup
    SE = strel('square',3);
    cup = imerode(cup,SE);
    %find edge of cup     
    cup = double(cup) + 1;
    cup = edge(cup,'canny');
    %     [B,cup]=bwboundaries(cup,'noholes');
    %     cup=bwmorph(cup,'bridge');
    cup =bwmorph(cup,'dilate');
    %aa=imread(cup);
    %figure,subplot(2,6,10),imshow(cup), title(' cup boundary');



I=cup;
    subplot(2,6,11),imshow(I), title('  boundary tracing');
    %column scan to select initial point
    BW=cup;
    [m,n]=size(BW);
    loc=[];
    flag=false;
    for row=1:m
        for col=1:n
            if BW(row,col)
                x=row;
                y=col;
                loc=[x y];
                flag=true;%%flag for other loops
                  break;
             end
        end
     if flag,break,end;    
    end
    row=loc(1,1);
    col=loc(1,2);
    row;col;
    contour = bwtraceboundary(I, [row, col], 'W',8,inf,'clockwise')
    if(~isempty(contour))
    hold on;
    plot(contour(:,2),contour(:,1),'g','LineWidth',2);
    hold on;
    plot(col, row,'rx','LineWidth',5);
    else
    hold on; plot(col, row,'rx','LineWidth',5);
    end
    %Find the number of pixels that make the outline
    n=length(contour)
    % there are n-1 distinct pixels as the first and last pixels are the same
    % initialize variables to display slope
    data=zeros(1,2)
    % Find the number of steps to "walk around" cup with a step size, lambda.
    % Determine step-size, lambda
    for number=1:5
   if number==1
      lambda=10;
   elseif number==2
      lambda=8;
       elseif number==3
      lambda=6;
           elseif number==4
      lambda=4;
           elseif number==5
      lambda=2;
           
   else 
        lambda=1;
   end
        % initialize variables
         steps=0;
         nextpt=zeros(1);
         step_dist=0;
         distance=0;
         initialpt=zeros(1);
         delta=zeros(2);
        % find pixels exactly lambda away from each other and 
        % calculate step distance then sum the step distances
        initialpt(1)=contour(1,2);
        initialpt(2)=contour(1,1);
        initialpt;
        for i=(lambda+1):lambda:n
            nextpt(1)=contour(i,2);
            nextpt(2)=contour(i,1);
            nextpt;
            delta(1)=initialpt(1)-nextpt(1);
            delta(2)=initialpt(2)-nextpt(2);
            step_dist=sqrt((delta(1))^2+(delta(2))^2);
            steps=steps+1;
            distance= distance + step_dist;
            initialpt=nextpt;
            initialpt;
         end
         lambda;
         steps;
         initialpt;
         % Find the distance btwn last pixel used above and  the last pixel that completes the outline and add to the distance if there are left over  pixels
         remainder=n-steps*lambda;
         if remainder>0
              nextpt(1)=contour(n,2);
              nextpt(2)=contour(n,1);
              nextpt;
              delta(1)=initialpt(1)-nextpt(1);
              delta(2)=initialpt(2)-nextpt(2);
              delta;
              r_dist=sqrt((delta(1))^2+(delta(2))^2); %remainder distance
              distance=distance+r_dist;
        end
        %  Store data
       data(number,1)=lambda;
       data(number,2)=distance;
    end
    % Displays lambda size with corresponding distance
    data, disp(' lambda, distance') 
    % plot the log steps versus log of lambda
    % find log of values
    data
    log_data(:,1)=log10(data(:,1));
    log_data(:,2)=log10(data(:,2));
    %log_data
    % Graph points (ln(perimeter) vs ln(of pixels skipped))
    x=0; y=0;
    x=log_data(:,1);
    y=log_data(:,2);
    %axes(handles.axes4)
    %imshow(plot(x,y));
    %subplot(2,6,12),plot(x,y),
    %title('log10(distance) vs log10(lambda)')
    % fractal dimension, D, is the slope of plot
    P=polyfit(x,y,1);
    hold on
    message=sprintf('Fractal dimension is %2.4f',1-P(1,1));
    z=1-P(1,1)
   zz = num2str(z);
 if(z>=1.04)
     set(handles.edit3,'string','Not Affected by Glaucoma');
     set(handles.edit2,'string',zz); 
 end
    if(z<=1.04)
       set(handles.edit3,'string','Affected by Glaucoma');
     set(handles.edit2,'string',zz); 
    end 
    hold off

    
    
    
    
    
    
    
    
    
    
    
    

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=get(handles.edit1,'string');
a=imread(str);
axes(handles.axes1)
imshow(a);



F= imresize(a,[256,256]);      
    figure;
    subplot(2,5,1),imshow(F), title('original image');


%% PREPROCESSING


    he=F;
    B =   strel('disk',15)         %structural element
    he=imclose(he,B);               % to remove blood vessels
    subplot(2,6,2),imshow(he), 
    title('preprocessed image');


%% SEGMENTATION OF CUP USING KMEANS


    % Creating color transformation from sRGB to L*a*b 
    cform = makecform('srgb2lab');
    % Applying above color transform to the sRGB image 
    lab_he = applycform(he,cform);
    % Converting into double 
    ab = double(lab_he(:,:,2:3));
    % obtaining rows and columns of transformed image 
    nrows = size(ab,1);
    ncols = size(ab,2);
    % Reshaping image taking each value column wise 
    ab = reshape(ab,nrows*ncols,2);
    % No of clusters to be created with five iterations 
    nColors = 5;
    [cluster_idx cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean',...
                                          'Replicates',5);
    % Reshaping and showing the clusters 
    pixel_labels = reshape(cluster_idx,nrows,ncols);
    subplot(2,6,3), imshow(pixel_labels,[]),
    title('image labeled by cluster index');
    % creating five element array 
    segmented_images = cell(1,5);
    % Creating tiles for three different colors 
    rgb_label = repmat(pixel_labels,[1 1 3]); 
    % Assigning clustered objects to array(segmented_image) 
    for k = 1:nColors
        color = he;
        color(rgb_label ~= k) = 0;
        segmented_images{k} = color;
    end
    % displaying different cluster objects 
    subplot(2,6,4),imshow(segmented_images{1}), title('cluster 1');
    subplot(2,6,5),imshow(segmented_images{2}), title('cluster 2');
    subplot(2,6,6),imshow(segmented_images{3}), title('cluster 3');
    subplot(2,6,7),imshow(segmented_images{4}), title('cluster 4');
    subplot(2,6,8),imshow(segmented_images{5}), title('cluster 5'); 
    %automatic selection of cluster containing cup by conveting into lab
     lab_he1 = applycform(segmented_images{1},cform);
     lab_he2= applycform(segmented_images{2},cform);
     lab_he3= applycform(segmented_images{3},cform);
     lab_he4= applycform(segmented_images{4},cform);
     lab_he5= applycform(segmented_images{5},cform);
    % calculating max point of lab of segmented image   
     m1=max(lab_he1(:));
     m2=max(lab_he2(:));
     m3=max(lab_he3(:));
     m4=max(lab_he4(:));
     m5=max(lab_he5(:));
    %  image with greatest max point taken as cup
     a=[m1 m2 m3 m4 m5];
     temp=a(1);
     for pos=1:5
         if a(pos)>=temp
             temp=a(pos);
             j=pos;
         end
     end
    % convert segmented image into grey scale   
    cup = rgb2gray(segmented_images{j});
    % convert into black and white
    cup = im2bw(cup,0.5);
    subplot(2,6,9),imshow(cup), title('cup ');
    %remove unwanted portion of segmented cup
    SE = strel('square',3)
    cup = imerode(cup,SE);
    %find edge of cup     
    cup = double(cup) + 1;
    cup = edge(cup,'canny');
    %     [B,cup]=bwboundaries(cup,'noholes');
    %     cup=bwmorph(cup,'bridge');
    cup =bwmorph(cup,'dilate');
    subplot(2,6,10),imshow(cup), title(' cup boundary');
    
    
    
%%  FRACTAL DIMENSION FOR SEGMENTED CUP



    I=cup;
    subplot(2,6,11),imshow(I), title('  boundary tracing');
    %column scan to select initial point
    BW=cup;
    [m,n]=size(BW);
    loc=[];
     flag=false;
    for row=1:m
        for col=1:n
            if BW(row,col)
                x=row;
                y=col;
                loc=[x y];
                 flag=true;%%flag for other loops
                  break;
             end
        end
     if flag,break,end;    
    end
    row=loc(1,1);
    col=loc(1,2);
    row;col;
    contour = bwtraceboundary(I, [row, col], 'W',8,inf,'clockwise');
  
    if(~isempty(contour))
     hold on; plot(contour(:,2),contour(:,1),'g','LineWidth',2);
     hold on;plot(col, row,'rx','LineWidth',5);
    else
    hold on; plot(col, row,'rx','LineWidth',5);
    end
    %Find the number of pixels that make the outline
    n=length(contour)
    % there are n-1 distinct pixels as the first and last pixels are the same
    % initialize variables to display slope
    data=zeros(1,2)
    % Find the number of steps to "walk around" cup with a step size, lambda.
    % Determine step-size, lambda
    for number=1:5
   if number==1
      lambda=10;
   elseif number==2
      lambda=8;
       elseif number==3
      lambda=6;
           elseif number==4
      lambda=4;
           elseif number==5
      lambda=2;
           
   else 
        lambda=1;
   end
        % initialize variables
         steps=0;
         nextpt=zeros(2);
         step_dist=0;
         distance=0;
         initialpt=zeros(2);
         delta=zeros(2);
        % find pixels exactly lambda away from each other and 
        % calculate step distance then sum the step distances
        initialpt(1)=contour(1,2);
        initialpt(2)=contour(1,1);
        initialpt;
        for i=(lambda+1):lambda:n
            nextpt(1)=contour(i,2);
            nextpt(2)=contour(i,1);
            nextpt;
            delta(1)=initialpt(1)-nextpt(1);
            delta(2)=initialpt(2)-nextpt(2);
            step_dist=sqrt((delta(1))^2+(delta(2))^2);
            steps=steps+1;
            distance= distance + step_dist;
            initialpt=nextpt;
            initialpt;
         end
         lambda;
         steps;
         initialpt;
         % Find the distance btwn last pixel used above and  the last pixel that completes the outline and add to the distance if there are left over  pixels
         remainder=n-steps*lambda;
         if remainder>0
              nextpt(1)=contour(n,2);
              nextpt(2)=contour(n,1);
              nextpt;
              delta(1)=initialpt(1)-nextpt(1);
              delta(2)=initialpt(2)-nextpt(2);
              delta;
              r_dist=sqrt((delta(1))^2+(delta(2))^2); %remainder distance
              distance=distance+r_dist;
        end
        %  Store data
       data(number,1)=lambda;
       data(number,2)=distance;
    end
    % Displays lambda size with corresponding distance
    data, disp(' lambda, distance') 
    % plot the log steps versus log of lambda
    % find log of values
    data
    log_data(:,1)=log10(data(:,1));
    log_data(:,2)=log10(data(:,2));
    %log_data
    % Graph points (ln(perimeter) vs ln(of pixels skipped))
    x=0; y=0;
    x=log_data(:,1);
    y=log_data(:,2);
    subplot(2,6,12),plot(x,y),
    title('(distance) vs (lambda)')
    % fractal dimension, D, is the slope of plot
    P=polyfit(x,y,1);
    hold on
    message=sprintf('Fractal dimension is %2.4f',1-P(1,1));
    z=1-P(1,1)
    zz=num2str(z);
    fdzz=strcat('FD=',zz);
     set(handles.edit3,'string',zz);
      if(z>=1.04)
     set(handles.text1,'string','Healthy eye');
     set(handles.edit3,'string',fdzz); 
 end
    if(z<=1.04)
       set(handles.text1,'string','Affected by Glaucoma');
     set(handles.edit3,'string',fdzz); 
    end    
    hold off





    
    
    
    
    
    
    


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
