
                SMXL 2.3 : Simulink MatriX Library

                              G.Campa
                       (campa@cemr.wvu.edu)
                     West Virginia University             
         Mechanichal and Aerospace Engineering Department

                            October 2000  



  ABSTRACT  
  ========

  The Simulink MatriX Library (SMXL) consists in a collection of blocks 
  that allow (varying) matrices to be handled directly in Simulink.

  The signals that flow through Simulink connections, which in general are
  the only quantities allowed to vary with time, can be scalars or vectors.
  It is not permitted to carry a matrix through a connection, that is,
  we could certainly carry a p by n matrix as a vector of p*n entries, 
  but all the existing Simulink built-in blocks that will process this 
  p*n by 1 vector will still treat it just like a vector. 
  As a consequence, every time that we have to process the p by n
  matrix we have to decompose it into its elements. 
  These decompositions, and the successive compositions make it very 
  uncomfortable to deal with varying matrices in the Simulink environment.

  The smxl library solves this problem, allowing operations like 
  matrix multiplication, transposition, and pseudoinversion to be 
  performed without any decomposition or type conversion.

  The library contains several examples as well as a pdf paper that gives
  further details and explains several motivations for using the blocks.  
  A library of ROTATION MATRICES which can be very useful 
  for simulations of rigid bodies in 3D space, is included.
 
  All blocks are compatible with Real Time Workshop.
  Smxl has been developed under Matlab 5.3.1 and tested
  under several platforms.



  DOWNLOADING
  ===========

  SMXL files are all contained in the directory smxl, which 
  can be freely downloaded from the official mathworks fileexchange site:
  http://www.mathworks.com/matlabcentral/fileexchange/



  HOW TO INSTALL SMXL 
  ===================

  Note: the path-notations in this section are in Windows format.
  Replace this by the appropriate format for other operating systems.  

  1. Copy the folder smxl from the mathworks ftp site to a 
     suitable destination folder.
     NOTE: SMXL requires less than 1Mbyte of disk space!
  2. Now start Matlab. 
     The Matlab search path should now be enhanced with the SMXL folders.
     It is recommended to add the SMXL root-folder to the Matlab path 
     permanently, since this simplifies future use. It is also suggested
     to add the EXAMPLES folder permanently to the matlab path.
     To add new folders to the matlab path use the PATH BROWSER utility,
     in the Matlab toolbar, (or the utility EDITPATH.M).
  3. When the SMXL folder is in the path you can access it from the 
     Simulink Library Browser, or by typing "smxl" from the matlab 
     command window. If the EXAMPLES folder is in the path as well,
     then you can access them with the matlab utility "demo",
     or by typing "dmxldemo" from the matlab command window.



  CURRENT VERSION INFO 
  ====================

  The most current version of the library is SMXL 2.3 (November 2004).

  In version 2.3, a much more efficient transposition block, based on the 
  selector block, was implemented (should have been done long time ago).
  Also, a couple of new blocks for aerospace engineering usage 
  was added in the Rotation Matices / Euler section of the library.

  Version 2.2 has been organized as a standard simulink library. 

  In Version 2.1 the code has been slightly optimized and a block
  for complex SVD has been added. 

  In version 2.0, a new code for pseudoinversion and SVD (which is 
  completely Real Time Workshop compatible) was added, along with new 
  blocks for complex matrices handling. The code for multiplication and 
  transposition was rewieved, and lots of new  examples were added.
  The library is now fully compatible with Real Time Workshop.


  CONTENTS 
  ========

  The SMXL folder contains the following DLL files:

  vrmult.dll, vrpinv.dll, vrsvd.dll,
  they are needed to run the blocks under windows platforms. 

  It contains as well the MDL file smxl.mdl which is the Matrix Library
                  and contains all the blocks that performs operations
                  on matrices. These blocks can be dragged and dropped
                  into any Simulink scheme.

  Finally there are 2 m-files :

  - Contents.m    which is useful to display basic info from the matlab
                  command prompt when typing ver or help smxl.
  - Slblocks.m    which allows the Simulink Library Browser to display
                  the smxl library in it.

  
  The EXAMPLES folder contains the following files :

  - Demos.m     : This m-file allows the allows the the matlab utility
                  "demo" to display the smxl demos in it.

  - smxldemo.mdl: This is the main Smxl Demo library, it can be used 
                  to call any example.

  - nnex.mdl    : This is a nonlinear control of an underwater vehicle,
                  which uses a self adaptive sigmoidal neural network,
                  (contained in the Anna block, under the mask).
                  The network has been implemented using real transposition
                  and real multiplication blocks.
  - abcdkkf.mat : This file contains the parameters needed to run the
                  nnex scheme. It is loaded from the nnex scheme
                  by double clicking on the "Load Model" block.

  - vccmspex.mdl: Example on complex matrices composition and splitting.
  - vcmultex.mdl: Example on complex matrices multiplication.
  - vcpinvex.mdl: Example on complex matrices pseudoinversion.
  - vcpsvdex.mdl: Example on complex matrices singular values decomposition.
  - vcrankex.mdl: Example on complex matrices rank.
  - vctrspex.mdl: Example on complex matrices transposition.
  - vrcmspex.mdl: Example on real matrices composition and splitting.
  - vrpinvex.mdl: Example on real matrices pseudoinversion.
  - vrrankex.mdl: Example on real matrices rank.
  - vrsvdex.mdl : Example on real matrices singular values decomposition.


  The SOURCE folder contains the c source code of the 4 s-functions:

  - vrmult.c    : real multiplication
  - vrpinv.c    : real pseudoinversion
  - vrsvd.c     : real singular value decomposition


  The DOC folder contains the following files:

  - index.txt   : This file contains the text included in the 
                  mathworks ftp INDEX.simulink file.
  - readme.txt  : The file you are reading now.
  - license.txt : This file explains some restrictions for copying
                  and distribution of SMXL, meant to prevent the 
                  distribution of incomplete versions.
  - vblocks.pdf : A paper which explains the motivation that lead
                  to the implementation of the blocks, as well as the
                  various motivations for using them. This paper
                  refers to the first version of the library.



  PLATFORMS
  =========

  - note for WINDOWS users: by default dll files are hidden,
    so you must set the option "show all files" in your explorer
    configuration if you want to see them, this is not necessary
    to use the package anyway.

  - note for NON WINDOWS users: you need to recompile the code 
    in the source directory and put the resulting executables
    in the main folder in order to use SMXL blocks.



  COPYING 
  =======

  The SMXL files may be copied freely under the restrictions put 
  forward in LICENSE.TXT (subfolder DOC of the toolbox). In short, 
  these restrictions are meant to prevent the distribution of incomplete
  or incompatible versions of SMXL; they are NOT intended to limit 
  the public accessibility of the toolbox in any way!



  CONTACT ADDRESS
  ===============

  For more information about the SMXL package, please contact:
  Giampiero Campa, Ph.D., Research Assistant Professor
  West Virginia University, Mechanical & Aerospace Engineering
  College of Engineering and Mineral Resources,
  Morgantown WV 26506-6106  USA
  E-mail      : campa@cemr.wvu.edu
                gcampa@supereva.it
  Note: these addresses may be subject to change.

  Thanks to M.O.Rauw for allowing the use of the FDC 1.3 files 
  readme.txt and licence.txt as a template for the respective SMXL files. 

  ========================================================================
  The SMXL library. Copyright G.Campa, 1999-2000. All rights reserved.