# TurboBC
TurboBC is a memory efficient and highly scalable GPU-based set of Betweenness Centrality (BC) algorithms in the language of linear algebra. These
algorithms are applicable to unweighted, directed and undirected graphs represented by sparse adjacency matrices in the Compressed Sparse Column (CSC) format, 
and in the transpose (COOC) of the Coordinate (COO) format. 
 
The TurboBC algorithms were designed to process sparse adjacency matrices selected from the SuiteSparse Matrix Collection in the Matrix Market format. 

More details of the design, implementation and experimental results obtained with the TurboBC algorithms are given in our paper cited below.

# Prerequisites
This software has been tested on the following dependences:
* CUDA 10.1.243 (GPU NVIDIA Titan Xp)
* gcc 8.4.0 
* Ubuntu 16.04.6

# Install
Instructions to install and compile the code:

1. Download the software: git clone https://github.com/pcdslab/TurboBC.git. The graphData folder contains some examples of Matrix Market files. You may add more files from https://suitesparse-collection-website.herokuapp.com.
2. Compile the code with the Makefile in the TurboBC folder. Please set the library and include directories paths on the Makefile available for your particular   machine.

# Run
1. Update the path of the graphData folder in the corresponding .sh file
2. Add permision for the .sh file, example: chmod +x run_bcugcsccooc.sh
3. run the experiments with: ./run_bcugcsccooc.sh

# Publications
If you use this software please cite our paper:

Oswaldo Artiles and Fahad Saeed,. 2021. TurboBC: A Memory Efficient and Scalable GPU Based Betweenness Centrality Algorithm in the Language of Linear Algebra. In 50th International Conference on Parallel Processing Workshop (ICPP Workshops ’21), August 9–12, 2021, Lemont, IL, USA. ACM, New York, NY, USA, 10 pages. https://doi.org/10.1145/3458744.3474047. https://camps.aptaracorp.com/ACM_PMS/PMS/ACM/ICPPWORKSHOPS21/35/d12b464b-debf-11eb-8d84-166a08e17233/OUT/icppworkshops21-35.html

# Acknowledgements
This research was supported by the National Science Foundations (NSF) under the Award Numbers CAREER OAC-1925960. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Science Foundation. We would also like to acknowledge the donation of a K-40c Tesla GPU and a TITAN Xp GPU from NVIDIA which was used for all the GPU-based experiments performed in our paper.

