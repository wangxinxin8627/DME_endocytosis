This is the project used in the paper 
"A Master Equation Approach to Actin Polymerization Applied to Endocytosis in Yeast", by Xinxin Wang and Anders Carlsson, submitted to PLoS Computational Biology

Compiling/running procedure:
(A) For Windows users
1. Install Microsoft Visual Studio 2013 (VS)
2. Download the zip package of the project
3. Decompress the package to get the folder "DME_endocytosis-master"
4. Double click the solution file ("calculation_RK.sln") in the folder
5. Now the project should be opened in VS
6. Click the "play" button to debug and run the program

(B) For Mac, Linux users
1. Download the zip package of the project
2. Decompress the package to get the folder "DME_endocytosis-master"
3. Open the terminal and go to the directory: "DME_endocytosis-master\calculation_RK"
4. Compile the file "for_linux_user.cpp" by typing in the terminal:
   g++ for_linux_user.cpp -o run
5. Run the program by typing
   ./run

after finishing, the output files will be in "DME_endocytosis-master\calculation_RK":
"data_timeCourse.txt" and so on.

membrane profiles are stored in '\calculation_RK\membrane', including:
1. nucleation, branching and forbidden regions
2. membrane related parameters
3. shape functions 

See more detailed description in the paper 
