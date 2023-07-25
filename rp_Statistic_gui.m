function varargout = rp_Statistic_gui(varargin)
%   varargout = rp_Statistic_gui(varargin)
%   Statistical Analysis GUI.
%   By YAN Chao-Gan 100401.
%   State Key Laboratory of Cognitive Neuroscience and Learning, Beijing Normal University, China, 100875
%	http://www.restfmri.net
% 	Mail to Authors:  <a href="ycg.yan@gmail.com">YAN Chao-Gan</a>
%	Version=1.0;
%	Release=200100401;
%--------------------------------------------------------------------------

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rp_Statistic_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @rp_Statistic_gui_OutputFcn, ...
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

function rp_Statistic_gui_OpeningFcn(hObject, eventdata, handles, varargin)

% Make Display correct in linux - YAN Chao-Gan 111025 Added.
if ~ispc
    ZoomFactor=0.85;
    ObjectNames = fieldnames(handles);
    for i=1:length(ObjectNames);
        eval(['IsFontSizeProp=isprop(handles.',ObjectNames{i},',''FontSize'');']);
        if IsFontSizeProp
            eval(['PCFontSize=get(handles.',ObjectNames{i},',''FontSize'');']);
            FontSize=PCFontSize*ZoomFactor;
            eval(['set(handles.',ObjectNames{i},',''FontSize'',',num2str(FontSize),');']);
        end
    end
end

handles.output = hObject;

guidata(hObject, handles);

function varargout = rp_Statistic_gui_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function btnTtest1_Callback(hObject, eventdata, handles)
rp_ttest1_gui;

function btnTtestpaired_Callback(hObject, eventdata, handles)
rp_ttestpaired_gui;

function btnTtest2cov_Callback(hObject, eventdata, handles)
rp_ttest2cov_gui;

function btnPCA_Callback(hObject, eventdata, handles)% Add PCA module SUN Jia-Wei 200502
get_citation('pca');
rp_PCA_gui;

function btn_MixedEffectsAnalysis_Callback(hObject, eventdata, handles)% Add MixedEffectsAnalysis module SUN Jia-Wei 200717
rp_MixedEffectsAnalysis_gui;

function btnAncova_Callback(hObject, eventdata, handles)
rp_ancova1_gui;

function btnCorr_Callback(hObject, eventdata, handles)
rp_corr_gui;

% --- Executes on button press in normalize_pushbutton.
function normalize_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to normalize_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rp_NormalityTest_gui;

% --- Executes on button press in Alpha_sim_pushbutton.
function Alpha_sim_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Alpha_sim_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rp_AlphaSim_gui;

% --- Executes on button press in smooth_pushbutton.
function smooth_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rp_Smoothest_gui;

% --- Executes on button press in GRF_pushbutton.
function GRF_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to GRF_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rp_GRFCorrection_gui;



% --- Executes on button press in TwoSample_EffectSize.
function TwoSample_EffectSize_Callback(hObject, eventdata, handles)
% hObject    handle to TwoSample_EffectSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
es_twoSampleES_image_gui;

% --- Executes on button press in Paired_Effect_Size.
function Paired_Effect_Size_Callback(hObject, eventdata, handles)
% hObject    handle to Paired_Effect_Size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
es_PairedES_image_gui;

% --- Executes on button press in IntraclassCorrelation_pushbutton.
function IntraclassCorrelation_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to IntraclassCorrelation_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ICC_GUI;







