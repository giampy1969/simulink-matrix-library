function tbxStruct  = demos

%DEMOS Demo list for the Simulink Matrix Library (smxl).

% Version 2.3
% Giampiero Campa 13-Nov-2004

if nargout == 0, demo blockset; return; end

tbxStruct.Name='Matrix Library';
tbxStruct.Type='Blockset';

tbxStruct.Help={
   
   'The Simulink MatriX Library (SMXL) consists in a collection'
   'of blocks that allow (varying) matrices to be handled directly'
   'in Simulink, it was built for R11.1 but works fine in R12'
   ' '
   'The smxl library allows matrix operations like multiplication'
   'transposition, and pseudoinversion to be performed in simulink' 
   'without any decomposition or type conversion.'
   ' '
   'A library of ROTATION MATRICES which can be useful for'
   'simulations of rigid bodies in 3D space is also included.'
   'All blocks are compatible with Real Time Workshop.'
   
};

 tbxStruct.DemoList = { ' Real Matrices Composition',       'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;vrcmspex';
                        ' Complex Matrices Composition',    'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;vccmspex';
                        ' Complex Matrices Transposition',  'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;vctrspex';
                        ' Complex Matrices Multiplication', 'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;vcmultex';
                        ' Real Matrices Pseudoinversion',   'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;vrpinvex';
                        ' Real Matrices SVD',               'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;vrsvdex';
                        ' Real Matrices Rank',              'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;vrrankex';
                        ' Complex Matrices Pseudoinversion','smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;vcpinvex';
                        ' Complex Matrices SVD',            'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;vcpsvdex';
                        ' Complex Matrices Rank',           'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;vcrankex';
                        ' Adaptive Neural Network',         'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;nnex';
                        ' Discrete Variable State Space',   'smxlexamplespath=which(''smxl'');path(path,[ smxlexamplespath(1:length(smxlexamplespath)-8) ''examples'' ]);clear smxlexamplespath;dvss' };
