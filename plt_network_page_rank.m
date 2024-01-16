function plt_network_page_rank(B1,pos,color,label,x_coord,y_coord)
% B1:incidence matrix
% pos: indices of the edges to be highlighted
% color: colors of the edges to be highlighted
% label: labels of the edges to be highlighted
for i = 1:size(B1,2)
   tb1(i,1)  = find(B1(:,i)==1);
   tb1(i,2) = find(B1(:,i)==-1);
end
tb1(:,3) = 1;

G = digraph(tb1(:,1),tb1(:,2),tb1(:,3));
%h = plot(G,'NodeColor','#0072BD','EdgeColor','#0072BD','LineWidth',0.5);

graycolor= gray(256);
c1 = graycolor(100,:);
c2 = graycolor(50,:);
%h = plot(G,'NodeColor',c,'EdgeColor',c,'LineWidth',0.5);
if nargin == 4
    h = plot(G,'NodeColor',c1,'EdgeColor',c2,'LineWidth',0.5);
elseif nargin == 6
    h = plot(G,'XData',x_coord,'YData',y_coord,...
        'NodeColor',c1,'EdgeColor',c2,'LineWidth',0.5);    
end
h.NodeLabel = {};
h.EdgeLabel = {};

% find the nodes and edges on the positions of pos
for j = 1:size(pos,2)
    node_label1 = tb1(pos(:,j),1);
    node_label2 = tb1(pos(:,j),2);
    for i = 1:length(pos)
        highlight(h,node_label1(i),node_label2(i),'EdgeColor',color(i,j),'LineWidth',10);
        % highlight the main edge index (the indicator edge), which is always the edge itself
        if label == ""
            labeledge(h,node_label1(i),node_label2(i),"");
        else
            labeledge(h,node_label1(i),node_label2(i),label(i));
        end

    end
end
h.EdgeFontSize=12;
box off;
set(gca,'visible','off')

end 
