# HMM-Target-Detection-Tracking-UAV

**Final year thesis project in QUT, working with Nova Systems Australia Pty Ltd**

**Testing Dataset provided by Nova Systems**

**Project Summary**

 The goal of the project was to implement a state of the art mechanism for detecting and tracking airborne targets from the perspective of a host UAV for a proposed Unmanned Traffic Management (UTM) system, where only SWIR/MIR/Electro-optical imagery source was considered as dataset.
- The aforementioned mechanism was implemented using Hidden Markov Model(HMM) filters. To gain more insight into what Hidden Markov Model is and how it works, please have a look here https://en.wikipedia.org/wiki/Hidden_Markov_model
- The coding aspect of this project was done using MATLAB
- Please have a read through the final thesis submission document in the current repository to get an overall insight into the procedures and outcomes of the project

**Executing the code**
1. To execute and evaluate the project's code, please download the Target Detection HMM.m file from the current repo
2. Download the folder named DataSet, as it contains some sample images for evaluating the performance of the HMM mechanism
3. Please copy the path of the DataSet folder in line **3 and 5** accordingly, code will be ready for execution without any bugs
4. Once the execution is finished, the output will be figures of the **transition matrix, observation probability matrix** and most importantly, figures containg the result- **estimated states for a given target** along with a output **CSV file** containing the numerical values of the **estimated states/pixel locations** of the target across the image frames.
