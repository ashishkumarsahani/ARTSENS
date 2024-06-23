%Scope_Link dll test

opengl info
niScope_Link(1,5195,'Dev3');
pause;
[frame,last_Time_Stamp] = niScope_Link(2,5195);
current_Time_Stamp = 0;
h=figure;
set(0, 'DefaultFigureRenderer', 'OpenGL');
set(gcf,'Renderer','OpenGL');
set(gcf,'RendererMode','manual');
axis manual
xlim([1,5195]); ylim([-4,4]); hold on;
tic
for i = 1:1005
    [frame,current_Time_Stamp] = niScope_Link(2,5195);
    int32(current_Time_Stamp-last_Time_Stamp)
    last_Time_Stamp = current_Time_Stamp;
    cla;
    plot((1:5195),frame);
    drawnow;
end
toc
niScope_Link(3,5195);