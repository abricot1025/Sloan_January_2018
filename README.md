# MOONs
Repository containing the code for MOONS (Multi-Object Optical and Near-infrared Spectrograph) collision avoidance.    
  
To better understand the state of the project, the chapter [Previous authors](#previous-authors) give an overview of what was done and the next point open for further research that the researcher has found or discussed at the end of their work.  
A good practice would be to keep updating this section at the end of each work ! 

## Introduction
The Multi Object Optical and Near-infrared Spectrograph for the VLT (MOONS) is a multi-scope instrument that is also aiming to probe galactic archaelogy, galaxy evolution and cosmology. It will obtain near infrared spectra for millions of stars and galaxies.    

MOONS will be installed on the VLT 8-meter telescope at Paranal Observatory in 2019. MOONS is a european project led by the Royal Observatory of Edinburgh.    

The MOONS fiber system will consistutes of 1000 fiber-positioning actuators following an hexagonal distribution on the focal plane with a ~3cm distance between actuators (center to center).     

This project involves groups from different countries (Switzerland, USA, Scotland, Italy...etc).   In europa, at least for the collision avoidance, different people are involved:
- EPFL: Responsible for collision avoidance, mechanical design, motor charaterization...etc
- Italy: responsible for creating the OPS file, a files defining the position in the focal space for each positioners and their target, and also assigning the parity and target to each positioner
- Scotland: Leader of this project, validation of code and research     

At the EPFL, three laboratories collaborate together with the Professor Jean-Paul Richard Kneib being the director of this project :  
- Prof. Jean-Paul Richard Kneib, Director of the LASTRO Lab. at the EPFL
- Dr. Mohamed Bouri, group leader of the LSRO Lab. at the EPFL
- Dr. Denis Gillet, Mer. Head of the REACT lab. at the EPFL  

While the others at the EPFL (LASTRO) are in charge of the design of the positioners (mechanical design, motors, ...), the REACT lab is specialized in multi-agent coordination.

## First step: 
### Administration:
When beginning the project, a first good step is to get introduced to all the people at the EPFL involved in this project.
- Ask Denis Gillet to get introduced to the member LASTRO (Astrobot) and also to Prof. Jean-Paul Richard Kneib and Dr. Mohamed Bouri
- Prof. Jean-Paul Richard Kneib or Dr. Mohamed Bouri (or a member of the LASTRO lab) to get introduced to Dr. Steven Beard, astronomer at the Royal Observatory of Edinburgh, Scotland.

### Technical:
##### Thesis, Conference paper and presentation 
A good introduction and explanation of the MOONs project are given by the following works concerning multi-agent coordination:
- Laleh's thesis : "decentralized multi-robot coordination in crowded workspaces" 
  - (DNF algorithm valid for DESI but not for MOONs)
- Paper : "Collision avoidance in next-generation fiber positioner robotic systems for large survey spectrographs"
- Paper : "Collision-Free Coordination of Fiber Positioners in Multi-object Spectrographs"
- Paper : Dominique Tao's soon to arrive
- presentation : [Link to presentation](https://taodominique.github.io/PrColAvoid_MOONsPresentation/#/) (or find it in presentation folder)   
Give a good overview of the Finite state machine logic

These paper can be found either in the "papers" and "presentation" folders from this github or online.

##### Datasheet:

- "RFE to Software ICD.docx" :  
  - gives a overview of the two differents referentials used for the focal and local plane of the positioners (the direction of the focal and local referentiels are not the same, ask astronomers or LASTRO guys why...)
- "VLT-TRE-MON-14620-3010 Fibre Positioner Software Test Cases.pdf" :  
  - Provided by Steven Beard from Scotland, a description of different test cased that are used to validate the robustness of the multi-agent coordination algorithm against collision and deadlock ! 

## Code: "MOONs":
A simulator of the positioners movement in the hexagonal focal plane is written in python (either 2.7 or 3.) (and soon in C++ as they want the project to be in C++). 
In order to have a good feeling about how it is organize, a brief description on how to launch the program, the different used files and results is given. A summary of the code architecture is also shown.
 
### start up ! 
With a good IDE in hand (I used pycharm and works like a "charm" ! ): 
1. in /PATH ANALYSIS Left Coordinate_Commited to SVN/, run "pa_interface.py" to open the following gui interface:  
![GitHub Logo](/imageMOONs_Readme/interface.png)  
 - [Configuration file](#configuration-file):
 	- .cfg file
 	- is located in PATH ANALYSIS Left Coordinate_Commited to SVN/from Stefano/"..."pos/"..."Positioners.cfg
 - [Target file](#target-file)
 	- .txt file
 	- is located in PATH ANALYSIS Left Coordinate_Commited to SVN/from Stefano/"..."pos/target/
 - two check box 
 	- Create a visual animation of positioners moving to their targets
 - "result folder" button
    - Create the result files indicating, with the last column, if the corresponding positioners managed to converge to its target (1) or was not able to (0)
    - Better create a result folder to dump all this result files...
 - "Go" button:
    - Launch the simulation
2. For collision avoidance, mainly work on "pa-dnf.py"

### Configuration file
The configuration file indicates how many positioners (with their targets) will be used for the python simulator! It also describe the polar focal coordinate of the positioner and its target in the focal plane.  
(Found after: "Specific parameters for each fibre positioner in the grid.")  
![GitHub Logo](/imageMOONs_Readme/configurationFile.png)  

### Target file
The target file is the one used to set the "positioner's movement parameters".  
It is organized as follows with polar coordinate (R,Theta) in the focal plane:  
![GitHub Logo](/imageMOONs_Readme/targetFile.png)  

As one could have noticed, we can find again the polar coordinate of the positioner and its target as it is necessary for the program to identify which positioner it corresponds in the configuration file    

The Parity is given by the pair (Arg 1, Arg2). Parity indicates the ability for a robotic arm to reach one position with different configuration depending on its degree of freedom (DOF). As we have only two DOF robots, it only has maximum two configurations where it can reach the target: either by coming from the right (Arg 1 = 0, Arg2 = 1) or from the left (Arg 1 = 1, Arg2 = 0)

The robot positioner priority indicates how important it is. The higher it is, the more chance it has to reach its target.

### Code architecture
![GitHub Logo](/imageMOONs_Readme/ArchitectureSystem.png)


## Test cases
Different test case exist:
1. in PATH ANALYSIS Left Coordinate_Commited to SVN/from Stefano
 - different folder with a "x" number of positioners.	
 	- in each of these folder, there are 9 differents files with the same number of positioners but with a different configuration
 	- NB: at the beginning, these 9 differents files were constructed from the ones with the 1006 positioner 
2. in VLT-TRE-MON-14620-3010 Fibre Positioner Software Test Cases/:
 - Different test cases to test the robustness of the code with scenario that can either end up with deadlocks or with the positioners never being able to reach the target
 - Provided by Steven Beard

## Result
In the "result" folder:
- in Video of before(using algorithm without Finite state machine and with)
- "AnalyisData.xlsx" 
	- show the result from running the code with and without the finite state machine on the different output configurations across different number of positioners ! 

## Conference to look for paper submission !!!
- SPIE. Astronomical Telescopes + Instrumentation 
- ICATOO
- ADASS Conference Series
- (Maybe others ?!?)

## Previous authors
- Dr. Laleh Makarem 
  - During her ph.d. thesis, she worked on “decentralized multi-robot coordination in crowded workspaces” (ph.d. thesis title). 
  - In order to efficiently coordinate the positioners of the spectograph in a non-collision manner, she studied, designed and used a “decentralized navigation function” (DNF), a family of potential functions that have only one minimum point. Her algorithm that ensures no collision among positioners is valid for the DESI project but not for the MOONs' one since the MOONs positioners have their second arm two times longer the the DESI's one. (bigger workspace and thus higher chance of collision/overlapping)  
  - Next tasks: 
  	- Make all the MOONs positoners converge (complete)
  	- Give priority to each positioner, as some distant object in the universe are more important to observe that others
  	- Solve deadlocks 

- Dominique Tao
  - Working as a scientist intern for 6 months, he extended the algorithm of Dr. Laleh Makarem for the MOONs project. He works on priority integration and collision avoidance.
  - In the study of the evolution of the universe, certain distant objects in the universe are more important to observe than others. Priorities (given beforehand by the OPS file)  are therefore assigned to each positioner. The higher the priority, the more important the object in the universe is to observe. 

  - He designed a local finite state machine imitating the principle of priority planning in order to integrate priority into the positioners movements such that higher prority positioner motion is favor to have a higher chance to reach its target. 
It also help  to locally solve deadlocks as they usually happen only in some them (the rest of the positioner being able to reach their target with the DNF) (despite the algorithm not being complete)
If  Laleh's work can be considered as a first layer making the positioners move, Dominique's one is the “decision” or “rules” layer coming on top of the first one, allowing the positioner to have a certain behaviours in specific situations . 
  - Mainly work on "pa_dnf.py" !
  - Results:
    - Improvement of results, convergence rate from 70 to 90 % overall!
    - a good step for priority integration, as priority algorithm works most of the time when positioners with higher priority encounter a lower one or two positioners with same priority meet each other but still not working all the time... 
  - Next tasks:
  	- Improve algorithms ! even if there is 20% of more convergence (from 65-70% to 85-90%) the algorithm is not complete
  	- Overparametrization... find a way to regulate the behavior of positiones without having to tune too much parameters ??? 
  	- Look for better parameterization value since sometime, positioners with higher priority still cannot bypass lower priority positioners...
  	- With the new project "SLOAN", one of the goal is to have three fiber in different EM spectrum (IR,visible,...)within the positioners in order to observe more objects in the universe. With this three optical fiber, investigate how the positioner should arrive on the target ? directly to its center and depending on if it recquires IR or visible optical fiber, move to the selected fiber? 
