
%GENERATING THE VIDEO DISPLAYING PENDULUMN DISTURBANCES 

%actual pendulumn coordinates
pend_length_y = 3:1:7 ;
pend_length_x = repelem(5,length(pend_length_y));

%mapping the pendulumn length to make reference point(i.e point of rotation) (0,0)
for i = 1:length(pend_length_y)
    mapped_coord{i} =  [mapfun(pend_length_x(i), 1, 9, -4, 4),mapfun(pend_length_y(i), 0, 7, -3, 4)];
end
mapped_coord = mapped_coord';
mapped_val = cell2mat(mapped_coord);


%rotating the pendulumn length for a given angle (theta)
pend_length_y = 0:1:mapped_val(end,2);
pend_length_x = repelem(0,length(pend_length_y));
plot(pend_length_x,pend_length_y,"Color",'red')
hold on

theta = -90; %angle of rotation

v = [pend_length_x(end);pend_length_y(end)];% getting the coordinates of the end of the pendulum
rot_theta = [cosd(theta) sind(theta); -sind(theta) cosd(theta)]; %rotation matrix along Y axis
x_y_rotval = rot_theta*v; %endpoint coordinates of rotated endpoint

rot_ylim = 0:0.01:x_y_rotval(2,1); %ylimits
if x_y_rotval(1,1) > 0
    rot_xlim = 0:0.01:x_y_rotval(1,1);
    plot([0,rot_xlim(end)],[rot_ylim(1) rot_ylim(end)],"Color",'blue')%xlimits
    title(['rotation through: ',num2str(theta),'degrees'])

elseif x_y_rotval(1,1) < 0
    rot_xlim = x_y_rotval(1,1):0.01:0;   
    plot([rot_xlim(end),rot_xlim(1,1)],[rot_ylim(1) rot_ylim(end)],"Color",'blue')%xlimits
    title(['rotation through: ',num2str(theta),'degrees'])

end
hold off
axis([-4 4 -3 4])


%% 

%ROTATING THE ACTUAL PENDULUM OVER THE RANGE -90 to 89 DEGREES

delta = -90:1:89; % range of rotation

v = [pend_length_x(end);pend_length_y(end)]; %x,y coord of the endpoint of the inverted pendulumn

close all
vid = VideoWriter('pendulumn_swing3.avi'); %video recorder object
open(vid);

%adding frames to the video onject
for indx = 1:length(delta)
    rot_theta = [cosd(indx) sind(indx); -sind(indx) cosd(indx)]; %rotation along y matrix
    x_y_rotval = rot_theta*v; %endpoint coordinates of rotated line

    %mapping the obtained rotated coordinates back to original endpoint values
    remapped_x = mapfun(x_y_rotval(2,1), -4, 4, 1, 9);
    remapped_y = mapfun(x_y_rotval(1,1), -3, 4, 0, 7);

    %plotting the frame containing the actual coordinates of ref image(rotation about point (5,3) )
    rot_ylim = 3:0.01:remapped_y; %ylimits

    check_y{indx} = rot_ylim; % catches errors in the loop
    check_x{indx} = rot_xlim;

    if remapped_x > 5 % for values on the right of reference line
        rot_xlim = 5:0.1:remapped_x;
        plot([5,rot_xlim(end)],[rot_ylim(1) rot_ylim(end)],"Color",'green')%xlimits


    elseif remapped_x < 5 % for values on the left of reference line
        rot_xlim = remapped_x:0.1:5;   
        plot([5,rot_xlim(1,1)],[rot_ylim(1) rot_ylim(end)],"Color",'black')%xlimits
    end
    hold on

%plot of the whole image showing all possible disturbances of the pendulumn
    hor_lim = 4:0.05:6;
    plot([hor_lim(1) hor_lim(end)],[3,3],"color",'black')
    hold on
    plot([hor_lim(1) hor_lim(end)],[1,1],"color",'black')
    hold on
    vert_lim = 1:0.05:3;
    plot([4,4],[vert_lim(1) vert_lim(end)],"Color",'black')
    plot([6,6],[vert_lim(1) vert_lim(end)],"Color",'black')
    hold on
    circle(4.5,0.5,0.5)
    hold on
    circle(5.5,0.5,0.5)
    axis([0 10 0 10])
    % set(gca,'visible','off')
    hold on
    axis([0 10 0 10])
    hold off

%extracting frames for videos

    frame = getframe(gcf); %get current figure handle
    writeVideo(vid,frame)

end %end of get video frame

close(vid)

%% Reading the video file created
close all
video_obj = VideoReader('pendulumn_swing4.avi');%interface to the video


%determining the reference line

frame_ip = read(video_obj,90);%frame corresponnding to a zero degree displacement
% imshow(frame_ip)
frame_ip = rgb2gray(frame_ip);%grayscale conversion


%removing the axes and reducing image size for faster computation
r_fr = centerCropWindow2d(size(frame_ip),[315 400]);
frame_crop = imcrop(frame_ip,r_fr);
% imshow(frame_crop);

%getting the point of rotation of the pendulumn in pixel coordinates
[x_up,y_up] = find(frame_crop < 200,1,'first'); %extracting the row location of pendulumn rotation axis
[x_up_pend,y_up_pend] = find(frame_crop(1:217,:) < 200,1,'first'); %extracting the column location
%point of rotation = [x_up,y_up_pend] == [218,210]

%determining pendulumn workspace
ref_line = frame_crop(1:217,210);
imshow(ref_line)
[R_ref,C_ref] = find(ref_line < 255,1,'first');

%% Reading frames from the input video

for count = 1:video_obj.NumFrames %video_obj.hasFrame()
    frame_inp = read(video_obj,count); 
    input_img = frame_inp ;
    input_img = rgb2gray(input_img);
           
    %cropping the images to improve computation perfomance

    r1 = centerCropWindow2d(size(input_img),[315 400]); %crop to 315 by 400
    cropped_ip = imcrop(input_img,r1);
%   imshow(cropped_ip);
 
   % determining quadrant the pendulum lies in
    [x_lim,y_lim] = find(cropped_ip(1:x_up-1,:)<200,1,'first');
    
    if y_lim < y_up_pend
        flag = 1; %negative_quadrant 
        
    elseif y_lim >= y_up_pend %tuned value of variable y_up_pend
        
        flag = 2; %positive quadrant
        
    else
        disp('pendulumn outside workspace')
    end
    
    % for pendulum in the positive region
    if flag == 2
        if count == 1
            mapped_input = cropped_ip(75:218,210:400); %75 is the upper pixel row value of the reference line
            imshow(mapped_input)
        else
           mapped_input = cropped_ip(75:217,210:400); %75 is the upper pixel row value of the reference line
           imshow(mapped_input) 
        end
    % for pendulum in the negative region
    elseif flag == 1
        if count < 179
            mapped_input = cropped_ip(75:217,1:210);
            imshow(mapped_input)
        elseif count >= 179
            mapped_input = cropped_ip(75:218,1:210);
            imshow(mapped_input)
        end
    end
    
      %mathematical approach to determine offset angle
    if flag == 2    
        [lower_r,lower_c] = find(mapped_input < 245,1,'first');
        [upper_r,upper_c] = find(mapped_input < 200,1,'last');
        b = (abs(upper_c-1)/abs(143 - upper_r));
        if b == 1 %point along reference line
            theta_calc = 0 ;     
        else 
            %ref point = (143,1) i.e lower left corner (tuned value)
            theta_calc = round(atand(abs(upper_c-1)/abs(143 - upper_r)));
        end
    elseif flag == 1
        im = flip(mapped_input, 2);
        [lower_x,lower_y] = find(im < 245,1,'first');
        [upper_x,upper_y] = find(im < 200,1,'last');
        theta_calc = -round(atand(abs(upper_y-1)/abs(143 - upper_x)));
    end
%     disp(theta_calc)
    theta_theory(count) = theta_calc;
    
    %EDGE DETECTION AND HOUGH TRANSFORM
    
    if flag == 1
      rot_IP = flip(mapped_input, 2); 
%     imshow(rot_IP)
      BW = edge(rot_IP,'canny'); %edge detection
    elseif flag == 2 
        BW = edge(mapped_input,'canny'); %edge detection
    end
    
    %hough transform
    [H,Theta,Rho] = hough(BW,'RhoResolution',1,'Theta',-90:0.5:89);
    
    %getting the highest peak
    P = houghpeaks(H,5,'threshold',ceil(0.9*max(H(:)))); %round off
    x = Theta(P(:,2)); %theta corresponding to the longest line
    y = Rho(P(:,1)); %rho corresponding to the longest line 

    %Assigning the necessary sign to detected angle based on the section
    %the pendulum line is in
    
    if flag == 1
        angle_displ = -round(x);
    elseif flag == 2
        angle_displ =  round(x);
    else
        disp('error')
    end
theta_img(count) = angle_displ(1,1);

end %for loop end

%%
%plotting a comparison between calculated angle and detected angle
comp = [theta_theory', theta_img'];

figure,
 t = (1:1:video_obj.NumFrames);
 plot(t, comp(:,1),'DisplayName','Theoretical angle values','color','green',"LineWidth",3); %calculated angle
 title("comparison between theoretical angle and detected angle")
 legend
 hold on
  plot(t, comp(:,2),'DisplayName','angles detected from image processing','color','black',"LineWidth",1); %result after image processing
  axis([0 180 -90 90])
  ylabel('angle values')
  xlabel('frame index')
  legend
 hold off





%Mapping function
function output = mapfun(value,fromLow,fromHigh,toLow,toHigh)
narginchk(5,5)
nargoutchk(0,1)
output = (value - fromLow) .* (toHigh - toLow) ./ (fromHigh - fromLow) + toLow;
end

 %circle plotting function
function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit,yunit,"Color",'black');
hold off
end
