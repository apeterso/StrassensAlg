StrassensAlg
============
The strassen application uses Strassen's algorithm to multiply two nXn matrices. It can also determine for what value of n Strassen's algorithm becomes more efficient than conventional "divide-and-conquer" matrix multiplication. It is also possible to multiply using a hybrid method of Strassen's algorithm and "divide-and-conquer" where the method begins using "divide-and-conquer" for all values of n lesser than a specified value.

The command line arguments are as follows: $ java -jar strassen.jar 0 dimension inputfile

0 is a placeholder for testing. However, if it is replaced by 1 the program will print out the value of n where larger values causes Strassen's algorithm to be more efficient than conventional "divide-and-conquer." Any value larger than 1 will be used as a "crossover point" where the "divide-and-conquer" method will begin to be used instead of Strassen's algorithm. Dimension specifies the size of the matrices to be multiplied. The input file must be standard ascii format containing enough integers to populate 2 matrices of the size specified in the dimension argument. Integer values in the input file must be listed on their own lines.
