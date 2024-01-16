function p1=plt_network(B1,f,x_coord,y_coord,width_ratio,color_value)
% B1: incidence matrix
% f: edge flow on top
% x,y_coord: coordinates of the nodes
% width_ratio: edge thickness to indicate the edge value
% color_value = 1/0: indicate the edge values by colors meanwhile
for i = 1:size(B1,2)
   tb1(i,1)  = find(B1(:,i)==1);
   tb1(i,2) = find(B1(:,i)==-1);
end
tb1(:,3) = 1;
tb1(:,4) = f;

% if the flow is negative, we reverse the edge
temp1 = tb1(find(f<0),2);
temp2 = tb1(find(f<0),1);
tb1(find(f<0),1) = temp1; 
tb1(find(f<0),2) = temp2; 

 
% M = [tb1(:,1), tb1(:,2), tb1(:,3)];
% cHeader = {'source','target','weight'};
% textHeader = strjoin(cHeader, ',');
% fid = fopen('chicago_list.csv','w'); 
% fprintf(fid,'%s\n',textHeader);
% fclose(fid);
% dlmwrite('chicago_list.csv',M,'-append','delimiter',',')

G = digraph(tb1(:,1),tb1(:,2),abs(tb1(:,4)));
LWidths = width_ratio*G.Edges.Weight/max(G.Edges.Weight);
if nargin == 4 % if there is no coordinate input
    p1=plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths);
elseif nargin == 6 % if there are coordinate inputs
    p1=plot(G,'XData',x_coord,'YData',y_coord,'LineWidth',LWidths);    
end
if color_value == 1
    G.Edges.value = G.Edges.Weight;
    G.Edges.EdgeColors = G.Edges.value;
    p1.EdgeCData = G.Edges.EdgeColors;
    colormap  
    colorbar('East')
end
box off;
set(gca,'visible','off')
p1.NodeLabel = {};
p1.EdgeLabel = {};
end 
