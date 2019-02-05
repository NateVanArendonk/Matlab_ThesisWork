function [a,b,c,d,e]=Scatter(varargin)

%%Scatter Plot with the Polynomial fitting of degree 1.
%%Get rid of writing separate commands for curve fitting and statistics
%%too.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scatter is a part of TimeSaving ToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Ankur Kumar                        Wednesday; January 03, 2018
%                                              Version: 1.0
%
% National Institute of Techonology, Rourkela
% Rourkela, Odisa - 769008, India
% Department of Earth and Atmospheric Sciences
% Email: ankurk017@gmail.com
%        416AS2025@nitrkl.ac.in
%
%
% Function:
%           One can draw the Scatter Plot with the figure handling options
%           along with the polynomial fitting of degree 1. This function
%           gives the basis statistics i.e. BIAS, Correlation Coefficient,
%           Standard Deviation of both variables, Root Mean Square Error
%           and write on the graph too. If you dont want statistics on the
%           plot, then you can change the format in the function itself.
%
% Syntax:
%           Scatter([2,3,6,4,7],[1,2,5,10,89])
%           Scatter([2,3,6,4,7],[1,2,5,10,89],20)
%
%
%Inputs:
%           First and second input should be the datas of scatter plot,
%           which you want to plot. First variable will plot on the X axis
%           and second variable will plot on the Y axis.
%           If you want to specify the size of marker size, write its size
%           as in the third argument.
%           Use LineWidth, MarkerFaceColor, MarkerFaceAlpha, LineWidth,
%           TextColor and Format to modify the Scatter Plots.
%           If you don't want the Statistics to print on your Scatter,
%           then write 0 as the Format argument.
%           There are 2 formats available presently, 1 and 2. Format 1 is
%           for the supressed statistics and Format 2 is the detail
%           Statistics.
%
%           For format 1, the equation of the fitting line is written on
%           the first line. Standard Deviation of X and Y is on the second.
%           SD of Y is in bracket. Correlation is also on the second line
%           seperated by slant. Bias and RMSE is on the third line.
%
%           You can store the output in variables to set the properties of
%           these lines later in the program.
%
%
% Example:
%           Scatter([2,3,6,4,7],[1,2,5,10,89],20,'LineWidth',1)
%           Scatter([2,3,6,4,7],[1,2,5,10,89],20,'LineWidth',1,'MarkerFaceColor','r')
%           Scatter([2,3,6,4,7],[1,2,5,10,89],20,'LineWidth',1,'MarkerFaceColor','r','LineColor','r')
%           Scatter([2,3,6,4,7],[1,2,5,10,89],20,'Format',1)
%           [a b]=Scatter([2,3,6,4,7],[1,2,5,10,89],20,'Format',1)
%           [a,b,c,d,e,f]=Scatter([2,3,6,4,7],[1,2,5,10,89],20,'Format',1)
%
%
% Please send your suggestions to the email id: ankurk017@gmail.com or
%                                               416AS2025@nitrkl.ac.in
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scatter is a part of TimeSaving ToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
    error('Enter at least 2 arguments to plot the Scatter')
end
AWS=cell2mat(varargin(1));
sat=cell2mat(varargin(2));
third_argument=10;
if nargin>=3
    if isnumeric(cell2mat(varargin(3)))
        third_argument=cell2mat(varargin(3));
    end
end
AW=AWS(:);
sat=sat(:);
id=union(find(isnan(AW)),find(isnan(sat)));
AW(id)=[];
sat(id)=[];
c=polyfit(AW,sat,1);
n=2*length(sat);
x1=linspace(min(AW),max(AW),n);
y1=(x1*c(1))+c(2);

MarkerFaceAlphaValue=0.8;
MarkerFaceColorValue='r';
LineWidthValue=5;
LineColorValue='b';
if nargin>=4
    id=find(cellfun(@ischar,varargin));
    id1=id(find(contains(varargin(id),'MarkerFaceAlpha')));
    if ~isempty(id1)
        MarkerFaceAlphaValue=cell2mat(varargin(id1+1));
    end
    id2=id(find(contains(varargin(id),'MarkerFaceColor')));
    if ~isempty(id2)
        if ~isletter(varargin(id2+1));
            MarkerFaceColorValue=cell2mat(varargin(id2+1));
        else
            MarkerFaceColorValue=varargin(id2+1);
        end
    end
    id3=id(find(contains(varargin(id),'LineWidth')));
    if ~isempty(id3)
        LineWidthValue=cell2mat(varargin(id3+1));
    end
end
a=scatter(AW,sat,third_argument,'filled',...
    'MarkerFaceAlpha',MarkerFaceAlphaValue,...
    'MarkerFaceColor',MarkerFaceColorValue);
hold on
id4=id(find(contains(varargin(id),'LineColor')));
if ~isempty(id4)
    if ~isletter(varargin(id4+1));
        LineColorValue=cell2mat(varargin(id4+1));
    else
        LineColorValue=varargin(id4+1);
    end
end
b= plot(x1,y1,'col',LineColorValue,'LineWidth',LineWidthValue);


hold on
get(gca,'XTick');
set(gca,'FontSize', 8,'fontweight','bold');
stats=SM_statistics(AW,sat);
if c(2)<0
    Line1=sprintf('y= %.2f x - %.2f',c(1),abs(c(2)));
else
    Line1=sprintf('y= %.2f x + %.2f',c(1),c(2));
end
std=stats{1};
corr=stats{2};
rmse=stats{3};
bias=stats{4};
id5=id(find(contains(varargin(id),'Format')));
formatindex=cell2mat(varargin(id5+1));
id6=id(find(contains(varargin(id),'TextColor')));
TextColorValue='k';
if ~isempty(id6)
    TextColorValue=char(varargin(id6+1));
end
if isempty(id5)|| formatindex==1
    Line2=sprintf(' %.3f (%.3f) / %.3f',std(1),std(2),corr);
    Line3=sprintf('%.3f / %.3f',rmse,bias);
    c=text(.65, .21, Line1,'units', 'normalized','Color',TextColorValue,'fontsize',8,'fontweight','bold');
    d=text(.65, .15, Line2,'units', 'normalized','Color',TextColorValue,'fontsize',8,'fontweight','bold');
    e=text(.65, .09, Line3,'units', 'normalized','Color',TextColorValue,'fontsize',8,'fontweight','bold');
    
    
elseif formatindex==2
    c=text(.15, .89, Line1 ,'units', 'normalized','Color',TextColorValue,'fontsize',8,'fontweight','bold');
    d=text(.15, .81,sprintf('STD= %.3f / %.3f',std(1),std(2)) ,'units', 'normalized','Color',TextColorValue,'fontsize',8,'fontweight','bold');
    e=text(.15, .74, sprintf('CORR= %.3f ',corr),'units', 'normalized','Color',TextColorValue,'fontsize',8,'fontweight','bold');
    f=text(.15, .66,sprintf('RMSE= %.3f ',rmse) ,'units', 'normalized','Color',TextColorValue,'fontsize',8,'fontweight','bold');
    g=text(.15, .58,sprintf('BIAS= %.3f ',bias) ,'units', 'normalized','Color',TextColorValue,'fontsize',8,'fontweight','bold');
    
    
elseif formatindex==0
else
    error('No further formats available. Use 0, if you don''t want statistics. Only Format 1 and Format 2 are available.')
    
end
end
function [ AQ] = SM_statistics( s_d,g_d)
%takes two series input (first is satellite and second is ground) and gives
% S.D., CC, RMSE,Bias
format short
a=[std(s_d) std(g_d)];
b=corr2(s_d,g_d);
c=sqrt(sum((s_d(:)-g_d(:)).^2)/numel(s_d));
d=nanmean(s_d-g_d);
AQ={a,b,c,d};
end
