function []=cufsmhelp(number)
%BWS
%October 2001 (last modified)
%Dec 2002 :last modified by Badri Hiriyur ( number 20)
%Almost all help messages are in this file, activated when appropriate button is pushed
%September 2010 updating for cufsm4
%
%1 = main help   
if number==1   
%commented out of cufsm.m in December 2002, added tooltips instead
msgbox(['      Load = Load a previous file (input only or completed anlysis)   ';...
        '      Save = Save a file                                              ';...
        '     Input = General input screen for analysis                        ';...
        'Properties = Section properties and tools for applying moment and load';...
        '   Analyze = Perform finite strip analysis                            ';...
        '      Post = View buckling curve and mode shapes from analysis        ';...
        '   Compare = Compare multiple analyses against one another            '],...
        'Help for Main Commands','help')
end
%main help in pre
if number==2   
%commented out of cufsm.m in December 2002, added tooltips instead
msgbox(['C and Z template = template for entering geometry of Cs and Zs';...
        'Double elem. = double the elements of the cross-section       ';...
        'Update and Replot = Update the variables for Input and replot '],...
        'Help for Input Commands','help')
end
%definition of prop
if number==3   
%CreateStruct.WindowStyle='';
%CreateStruct.Interpreter='tex';
msgbox(['Enter an identification number for the material (mat#) followed by the material properties',...
        'for the member, separate values by spaces. ',...
        'The material properties allow for orthotropic definition of E, v, and G. Simply enter in ',....
		'istoropic values for the needed quantities to ignore this.'],...
        'Help for material properties','help')
end
%definition of node
if number==4  
msgbox(['enter the node data, nodal coordinates, boundary conditions and stress, separate values by commas or spaces',...
        'node# = node number',...
        'x = x coordinate (left-right direction), ',...
        'z = z coordinate (up-down direction), ',...
        'xdof = degree of freedom in the x direction, 1 or 0,  1=free 0=fixed (generally leave as 1), ',...
        'zdof = degree of freedom in the z direction, 1 or 0, ',...
        'ydof = degree of freedom in the y (axial) direction, 1 or 0, ',...
        'qdof = degree of freedom in the q (q=theta, or twist) direction, 1 or 0, ',...
        'stress = stress at the node - use 1.0 if you want to ignore for now, and use the loading button to generate a stress distribution.'],...
        'Help for node entry','help')
end
%definition of elem
if number==5   
msgbox(['enter the element data, separate values by single spaces end with the mat# which ',...
        'refers back to the properties entered above.'],...
        'Help for element entry','help')
end
%definition of lengths
%this message superseded in cufsm4
if number==6   
msgbox(['enter the lengths to be considered, IN MATLAB ONLY, you may enter vectors as well as singular numbers ',...
		  'for example 10:10:100 is the same as entering 10 20 30 40 50 60 70 80 90 100 ',...
		  'you may also space things evenly in logspace, enter logspace(0,1,20) puts 20 points ',...
        'between 10^0 and 10^1 (evenly spaced in a semilogx plot). In the Standalone version in DOS enter the lengths separated by a space or a comma.'],...
        'Help for lengths entry','help')
end
%help on calculation of loads and moments using fy
if number==7   
msgbox(['After entering the value of the yield stress, pressing calculate will ',...
		  'generate the individual loads and moments for first yield. For the ',...
		  'bending calculation you may allow bending about the principal axes, ',...
          'or you may select restrained bending, which will only bend about the ',...
		  'x and z axes.',...
          'It also calculates the Bimoment developed in a simply supported ',...
          '(warping free ends) beam of length, L, with a concentrated torque, T, at position, x.'],...
          'Help for calculation of loads and moments','help')
end
%help on generation of stress based on loads and moments
if number==8   
msgbox(['The finite strip analysis is performed by defining the stress at each ',...
		'node. Select the button to the left to submit the stress generated on this page',...
        'to the analysis pre-processor. The stress distribution may be determined by entering in actions ',...
		'P, B, M etc. that you want to consider as the reference applied loads. Or ',...
        'you may also use the yield loads (or moments) for reference (press the Py button for example), or you ',...
		'can determine the actions that are closest to the current applied   ',...
		'stress by using the Generate stress button (or even read forces from Mastan..)'],...
        'Help for stress generation','help')
end
%help on plotting modes
if number==9   
msgbox(['The shape to the right is the buckling mode for the member at the length selected. ',...
        'The length should be interpreted as the half-wavelength for signature curve (traditional) analysis and the physical length for general boundary condition analysis.',...
		 'For 2D plots the deformation is shown at the longitudinal cross section position, i.e., y/L specified. It can be changed for convenience. ',...
		 '3D modes can be shown by checking 3D and then selecteing Plot Mode. ',...
         'Different lengths may be selected by pushing the arrows and then Plot Mode. ',...
         'If available, buckled shapes in higher modes may be selected in a similar fashion. ',...
		 'The scale of the plot may be changed (including negative numbers) for convenience. ',...
         'A box around a node indicates some degree of fixity at this node, a star around a node ',...
		 'indicates the presence of a spring at this node.'],...
         'Help for plotting modes','help')
end
%help on length
if number==10   
msgbox(['For a signature curve (traditional) finite strip analysis assumes a single half sine wave from end ',...
        'to end (longitudinal term m=1), the length of this half sine wave is known as the half-wavelength. ',...
        ' A large number of half-wavelngths are usually analyzed in order to understand the different possible buckling mode shapes.',...
	    'For general boundary condition solution, the length in finite strip analysis is the physical length. A large number of modes at a given physical lengths is usually analyzed, similar to shell finite element analysis. ',...
        'Press the arrow to select higher or lower half-wavelgnths and then press Plot Mode to ',...
		  'update the plot.       '],...
        'Help for length','help')
end
%help on mode selection
if number==11   
msgbox(['For signature curve (traditional) analysis in most situations the first buckling mode is of primary interest, and thus these controls',...
		  'are not needed. However, in some cases, the higher modes that exist at a given half-wavelength ',...
		  'are of interest - in this case press the arrows to select the mode of interest and ',...
        'then press Plot Mode to update the plot. An example where higher modes are of interest is ',...
		  'when the first mode is restricted in some way that is not included in the analysis. ',...
          'For general boundary condition solution, higher modes are always necessary to find all the three ',...
          'buckling modes (local, distortional, global).'],...
        'Help for mode selection','help')
end
%help on stress plot
if number==12   
msgbox(['Selection of the stress distribtuon button will provide a plot of the current stress ',...
		  'distribution on the member. To modify the stress distribution, select Input to modify the ',...
		  'values node by node, or select Sect. Prop. to place a new load or moment on the ',...
        'member. ',...
		  ' '],...
        'Help for stress plot','help')
end
%help on plotting buckling curve
if number==13   
msgbox(['The buckling curve summarizes the results of the analysis. The horizontal axis is the ',...
		  'length at which the analyses are performed at. Half-wavelength for singature curve and physical ',...
          'length for general boundary condition solution. The vertical axis is the load factor, or ',...
		  'the eigenvalue. The load factor is a multiplier times the applied stress distribution that ',...
         'indicates when the buckling occurs. For signature curve, the half-wavelength is the length of the single half ',...
		 'sine wave used to perform the analysis. Minimums of this plot are of particular interest. ',...
         'For general boundary condition solution, the curve of load factor vs mode number reveals more buckling info.',...
         'check the load factor vs mode # and then Plot Curve to examine the buckling behavior of a physical length for higher modes ',...
         'and its corresponding longitudinal term participation.',...
         ],...
        'Help for buckling curve','help')
end
%help on modes
if number==14   
msgbox(['Enter the mode number you want to plot on the curve of load factor vs length. For signature curve, in most situations ',...
        'the first buckling mode is of primary interest, and thus these controls ',...
		  'are not needed. However, in some cases, the higher modes that exist at a given half-wave',...
		  'length are of interest - in this case press the arrows to select the highest mode of ',...
        'interest, then press Plot Curve. An example where higher modes are of interest is when the ',...
		  'first mode is restricted in some way that is not included in the analysis. ',...
          'For general boundary condition solution, higher modes are always necessary to find all the three ',...
          'buckling modes (local, distortional, global).'],...
        'Help for mode selection','help')
end
%help on text output
if number==15   
msgbox(['Ascii Text output of the buckling curve of the lowest mode is available by selecting the ',...
		  'text output button.                                                                      '],...
        'Help for text output','help')
end
%help on inputting continuous springs
if number==16   
%version 4.2 and earlier
% msgbox(['Continuous elastic springs may be added to your model. Leave this entry blank if you do not ',...
% 		  'want to add any springs to the model. The required entry for each spring is the node number ',...
% 		  'where the spring will act, the global degree of freedom in which the spring will act where 1=x dir. (i.e., left-right) ',...
% 		  '2=z dir. (i.e., up-down) 3=y dir (i.e., axial) 4=q dir (twist), the stiffness of the spring. ',...
%           'and ''kflag'' to indicate if the entered value is the total stiffness (0) or a foundation stiffness (1). ',...
%           'For example 1 2 4.0 1 means node 1, dof 2, has a spring of 4 (force/length)/length. The total stiffness ',...
%           'added to the model will reflect the length of the current analysis step. ',...
%           'For example 2 3 60.0 0 means node 2, dof 3, has a spring of 60.0 force/length, the spring magnitude is a constant. ',...
%  		  'All springs are assumed grounded - so they cannot be used to connect two nodes. '],...
%          'Help for spring definition input','help')
msgbox(['**NOTE! Springs under active development, run your own verification before using.** Springs may be added to your model. The format of the input is [spring# nodei nodej ku kv kw kq local discrete y/L] ',...
		  'where spring# is an integer that labels the spring, nodei is n integer that indicates the starting node, ',...
		  'nodej is an integer that indicates the ending node, ku, kv, kw, kq are the u,v,w, and theta spring stiffness, ',...
		  'local is a logical referring to the coordinate system 1 means use the local coordinate system defined such that ',...
          'the x direction and ku are from nodei to nodej, a 0 for local means use the global coordinate system, discrete ',...
          'is a logical referring to the type of spring 1 means a discrete spring 0 means a foundation or continuous spring, ',...
          'y/L is the longitudinal locaiton of the discrete spring otherwise enter 0. Note to create a spring to ground instead ',...
          'of between nodes set nodej=0. Example [1 5 15 1.0 0.1 0.1 0.001 1 0 0] which is spring 1 connect node 5 to node 15 ',...
 		  'the ku stiffness is 1.0 the kv and kw stiffness are 0.1 and the kq stiffnes is 0.001 the local (1) coordinate system ',...
 		  'is used the spring are foundation springs (i.e. per length) and y/L location is thus irrelevant 0 used.'],...
         'Help for spring definition input','help')
end
%help on file selector in compareout
if number==17   
msgbox(['Press the arrows to select the number of the file that you wish to see the mode shape for. Then press Plot Mode.',...
		  ''],...
        'Help for selecting file','help')
end
%help on file selector in compareout
if number==18   
msgbox(['Enter the number of the files that you want displayed in the buckling curve, separate your entries with spaces. Then press Plot Curve.',...
		  ''],...
        'Help for selecting files','help')
end
%constraints modeled in the pre-processor
if number==19   
msgbox(['                                                                                                    ',...
        'Equation constraints may be added to your model. Leave this entry blank if you do not want to add ',...
		'any equation constraints to the model. Equation constraints are used to tie different degrees of ',...
        'freedom together. For instance say the vertical deflection at nodes 10 and 20 should be the same. ',...
        'this is achieved by (e)liminating one degree of freedom at one node and writing an equation that is ',...
        'in terms of the (k)ept degrees of freedom.         For the the example of w10 = 1.0*w20, the user should enter ',...
        '"10 2 1.0 20 2" which says node 10, degree of freedom 2 should be eliminated and set equal to 1.0 times ',...
        'node 20 degree of freedom 2. More complicated expression with multiple (k)ept degrees of freedom can be ',...
        'completed by adding additional columns. Eliminated degrees of freedom may not be used as kept degrees of ',...
        'freedom in any expression. Common use of this feature is to model an external strap where two portions ',...
        'of the model are forced to move together. Or, to partially restrict behavior to examine a particular mode. ',...
        '     MASTER_SLAVE: Use of the master-slave constraints button will allow you to slave any number of nodes to ',...
        'a master node. This allows rigid cross-section modes, and partial rigid cross-section modes to be modeled ',...
        'for example slave all nodes of a flange to the flange/web junction to enforce distortional buckling, or ',...
        'or slave all nodes to one single master node, to examine global torsional and flexural modes. All entries ',...
		'must have the same number of columns, trailing zeros may be added to the end of the entries.'],...
        'Help for constraints modeling','help')
end
%constraints modeled in the pre-processor
if number==20   
msgbox(['                                                                                                    ',...
        'This feature allows constraints to be added to any model in the spirit of GBT style calculations. ',...
        'View the constraints ',...
        'to see whether or not you desire to allow that particular deformation mode in your analysis. '...
        'One needs to choose the basis and activate the constrained Finite Strip Method (cFSM) by selecting ',...
        'the base vectors (1 means on). This allows users to experiment with modal decomposition and identification ',...
        'in an active way.',...
		  ''],...
        'Help for modal constraints modeling','help')
end
%warning on calculation of loads and moments using fy
if number==21   
msgbox(['Calculated loads and moments are based on the nodal coordinate locations, ',...
	    'NOT the extreme fiber. This approximation is generally valid for thin-walled ',...
	    'members, but the user should exercise care as the thickness increases.',...
        ''],...
          'Warning for calculation of loads and moments','help')
end
%help on modal classification
if number==22
msgbox(['This feature allows for determination of modal participationg factors using cFSM, ',...
	    'cFSM must be turned on in the pre-processor and the basis selected in the cFSM. ',...
	    'The modal participation also depends on the manner in which the cFSM base vectors are ',...
	    'normalized (in addition to the selected basis). Work norm, strain energy norm, and ',...
	    'a simple vector norm are available.',...
        ''],...
        'Help for classification by cFSM','help')
end
%




%help on solution type
if number==61
msgbox(['This feature allows selection of solution types between the signature curve (traditional solution) and general boundary condition solution. ',...
	    'The signature curve is a special case of the general boundary solution when longitudinal term employed is only m=1 ',...
	    'and end boundary conditions are set to be Simply-Simply supported (S-S). Usually the lengths in signature curve are seen ',...
	    'as a sweep of half-wavelengths. Thus, the solution of singature curve is in terms of load factor versus half-wavelength. ',...
	    ' While general boundary condition solution is load factor versus physical length similar to FEM solution, and usually a set of longitudinal terms should be included for each length. ',...
        'Carefully selecting the longitudinal terms (by Recommend m) is always advised to save computational time.',...
        ''],...
        'Help for solution type','help')
end
%
%help on boundary condition
if number==62
msgbox(['End boundary conditions can be specified as simple-simple (S-S), clamped-clamped (C-C), ',...
	    'simple-clamped (S-C), clamped-free (C-F), and clamped-guided (C-G). For signature curve, ',...
	    'end boundary conditions shall always be S-S.',...
        ''],...
        'Help for boundary conditions','help')
end
%
%help on eigenvalue number
if number==63
msgbox(['Specify the number of eigenvalues of the solution. For general boundary condition solution, ',...
	    'usually more than 10 eigen solutions are needed to possibly ensure the solution including all  ',...
	    'the three buckling modes (local, distortional, global) depending on the longitudinal terms ',...
        'employed in the solution.'],...
        'Help for number of eigenvalues','help')
end
%
%help on length and longitudinal term
if number==64
msgbox(['Enter the lengths and longitudinal terms for each length for finite strip analysis.   ',...
    '                                                            ',...
	    'For signature curve, finite strip analysis assumes a single half sine wave from end to end (longitudinal term m=1), the length of this ',...
        'half sine wave is known as the half-wavelength. A large number of half-wavelngths are ',...
        'usually analyzed in order to understand the different possible buckling mode shapes.',...
        '                                                                ',...
	    'For general boundary condition solution, the length in finite strip analysis is physical length. ',...
        'For each length, several longitudinal terms should be included to accurately capture the possible buckling behavior. ',...
        'Longitudinal shape function is plotted on the right to help understanding the longitudinal term. ',...
        '                                                               ',...
        'For signature curve, IN MATLAB ONLY, you may enter column vectors as well as singular numbers ',...
		'for example [10:10:100]'' is the same as entering [10; 20; 30; 40; 50; 60; 70; 80; 90; 100;] ',...
		'you may also space things evenly in logspace, enter logspace(0,1,20)'' puts 20 points ',...
        'between 10^0 and 10^1 (evenly spaced in a semilogx plot in column). ',...       
        'In the Standalone version in DOS enter the lengths separated by a semicolon.',...
         '                              ',...
        'For general boundary condition solution, each row is a length with associated longitudinal terms. ',...
        'You may enter the lengths similar to the way for signature curve without entering the longitudinal terms and longitudinal term 1''s are automatically assigned for each length.',...
         'Also one can either input/edit their own longitudinal terms for each length or by using the menu to make them all same ',...
        'or using the suggested longitudinal terms based on the signature curve. ',...
        ''],...
        'Help for length and longitudinal term','help')
end

%Help on the basis selection
if number==70   
msgbox(['Natural basis is defined by explicitly following the mechanical criteria similar to those in GBT,  which separates the deformations into the G, D, L, ST/O spaces . ',...
        'The modal basis (similar to GBT) by performing an auxiliary eigen problem within each space either for a unit axial stress, or for the actual applied stresses. ',...
        'For non-simply supported boundary conditions due to the loss of orthogonality of the stiffness matrices between longitudinal terms, ',...
        'whether the constrained eigenvalue problem is solved inside each longitudinal term or over all the longitudinal terms results in the uncoupled and coupled bases, respectively. ',...
        'Finally, ST space is the union of the shear and transverse extensions (generally preferred) and the O space is defined as the null of the GDL subspace, ',...
        'which also has the couple and uncoupled space for non-simply supported boundary conditions. ',...                                                                                                     ',...
		'In general, coupled basis is more accurate while uncoupled basis is more computational efficiency (preferred). ',...
        ''],...
        'Help for basis selection','help')
end
%Help on the base vector
if number==71   
msgbox(['                                                                                                    ',...
        'The lengths of base vector of each space is determined by the cross section of the member. ',...
        'For local and other/ST, they are also influenced by the mesh. Whenever the node and element have been changed, ',...
        'the base vectors should be accordingly updated as well. '...
        'For the base vectors, 1 means on. For modal classification these base vectors are not necessarily be active (save a little computational time). ',...
        'One can come back to activate them when performing modal classification. ',...
		''],...
        'Help for base vectors','help')
end
%Help on the m recommend--half-wave length
if number==72   
msgbox(['                                                                                                    ',...
        'Half-wave lengths of local and distortional buckling are necessary for recommedation of longitudinal terms. ',...
        'A signature curve with unique minima is defined as one that includes two distinct minima, corresponding to local and distortional buckling. ',...
        'In this case, the half-wave lengths of local and distortional buckling are distinct.',...
        'If either or both minima is ''indistinct'' the signature curve is characterized as having only non-unique minima. '...        
        'In this case, need user''s judgement to input the missing half-wave lengths. Or if the cross section of the member is straight-line model, ',...
		'one can turn on the cfsm for helping on the half-wave lengths. '],...
        'Help for half-wave lengths','help')
end




%
%About
if number==200   
msgbox(['                                                                                                    ',...
        '',...
        'CUFSM stands for Constrained and Unconstrained Finite Strip Method. CUFSM is authored by Ben ',...
        'Schafer with large contributions from Zhanjie Li, Sandor Adany, and Andrew Sarawit, and smaller ',...
        'contributions from many others. CUFSM is a cross-section analysis tool that provides a means to ',...
        'assess the stabiltiy of thin-walled members. CUFSM is available at www.ce.jhu.edu/bschafer/cufsm. ',...
        'CUFSM wants you to have fun with stability, so get at it.',...
        '',...
        '',...
        '',...
        '',...
        '',...
        '',...
        ''],...
        'About CUFSM')
end

%Extra tools
%Plastic Surface Builder
if number==300   
msgbox(['The plastic surface builder provides a means to numerically generate the plastic strength of ',...
        'a general thin-walled cross-section under any combination of longitudinal stress. The original ',...
        'code for the plastic surface builder was developed by Dr. Shahab Torabian. The method assumes ',...
        'elastic perfectly-plastic material behavior in its current incarnation. The method has its ',...
        'greatest use in the context of the generalized beam-column strength prediction.' ,...
        '',...
        '',...
        '',...
        ''],...
        'About Plastic Surface Builder')
end
%Hole Effect Analysis Tool
if number==400   
msgbox(['<Under Development.> ------------- The hole effect analysis tool allows you to use an approximate method to determine the influence ',...
        'of holes on the elastic buckling of a cross-section. The orginal code was developed by ',...
        'Junle Cai under the instruction of Dr. Cris Moen in 2016. The method uses different approximate techniques ',...
        'for local, distortional, and global buckling. The interface can handle multiple hole locations ',...
        'Note, in local buckling, buckling at the hole location and away from the hole must both be checked. ',...
        'In distortional buckling the impact of the hole must be assessed over the elastic distortional buckling ',...
        'halfwavelength without the holt (Lcrd). Lcrd is found from a separate finite strip analysis on the ',...
        'section without a hole. Global buckling uses a combination of net and smeared properties in the classical ',...
        'cubic equation to provide its solution. (i.e. finite strip analysis is not used for global buckling ',...
        'analysis, rather something like CUTWP is employed). See the report of Cai and Moen for futher details.',...
        '-----------------------------------------------------------------------------',...
        'Note, based on original coding for this tool the section must be constant thickness, open, no branches, ',...
        'and the numbering of nodes and elements needs to be sequential.',...
        ''],...
        'About Hole Effect Analysis Tool')
end
%Abaqusmaker
if number==500   
msgbox(['ABAQUS maker is legacy code created to provide a quick means to generate an 3D shell model in ABAQUS ',...
        'from a 2D CUFSM model. The code was standalone from CUFSM prior to version 5 of CUFSM, and provides a means ',...
        'for researchers to understand how to move from CUFSM to ABAQUS shell models. Modification of the base code ',...
        'is generally required for special situations.',...
        ''],...
        'About ABAQUS Maker')
end


