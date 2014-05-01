package progassign2;

import java.io.*;
import java.util.Random;
import java.util.Scanner;

/**
 * @author Anders Peterson
 * Programming Project 2
 * Due 4/5/2014
 */
public class Main {
    private static long strassOperations = 0;
    private static long convOperations = 0;
    private static Random random = new Random();
    
    /**
     * Main prints the diagonal entries of a matrix that is the product of two matrices
     * whose values are dictated by a chosen text file. The matrices are multiplied using
     * Strassen's algorithm. If 1 is the first argument, main will print the optimal
     * crossover point between Strassen's algorithm and conventional divide-and-conquer.
     * @param args  expected in the form 0 <dimension> <inputfile>
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException{
        if(Integer.parseInt(args[0]) == 1){
            System.out.println("Crossover point is " + getCrossover());
        }
        else{
            int dimension = Integer.parseInt(args[1]);
            
            double[][] firstMat;
            double[][] secondMat;
            if(args.length < 3){
                firstMat = buildPowTwoMat(generateNNMatrix(dimension));
                secondMat = buildPowTwoMat(generateNNMatrix(dimension));
            }
            else{
                Scanner scanner = new Scanner(new File(args[2]));           
                firstMat = buildPowTwoMat(generateNNMatrix(scanner, dimension));
                secondMat = buildPowTwoMat(generateNNMatrix(scanner, dimension));
            }
            
            double[][] productMat;
            if(Integer.parseInt(args[0]) != 0){
                productMat = strassen(firstMat,secondMat,Integer.parseInt(args[0]));
            }
            else{
                productMat = strassen(firstMat,secondMat);
            }
            
            for(int i = 0; i < dimension; i++){
                System.out.println(productMat[i][i]);
            }
        }
    }
    
    /**
     * The getCrossover method tests the Strassen matrix multiplication method and the
     * conventional matrix multiplication method for an increasing n.
     * 
     * @return the first power of two that causes Strassen's method to be more efficient
     *         than the conventional method
     */
    private static int getCrossover(){
        int n = 1;
        double[][] firstMat;
        double[][] secondMat;
        while(strassOperations >= convOperations){
            n = n*2;
            System.out.println(n);
            strassOperations = 0;
            convOperations = 0;
            firstMat = generateNNMatrix(n);
            secondMat = generateNNMatrix(n);
            strassen(firstMat, secondMat);
            divideAndConquer(firstMat, secondMat);
            System.out.println("conventional took " + convOperations);
            System.out.println("strass took " + strassOperations);
            System.out.println("ratio conv/strass " + (double) convOperations/strassOperations);
            System.out.println("-----");
        }
        return n;
    }
    
    /**
     * The generateNNMatrix method creates an nXn matrix containing values from the input
     * text file.
     * @param scanner   the scanner from which the vales will be drawn 
     * @param n         the dimensions of the nXn matrix
     * @return          a populated nXn matrix
     */
    private static double[][] generateNNMatrix(Scanner scanner, int n){
        double[][] mat = new double[n][n];
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(scanner.hasNext()){
                    mat[i][j] = scanner.nextDouble();
                }
                else{
                    throw new RuntimeException(
                            "Not enough values in the input file for the specified"
                            + " dimension");
                }
            }
        }
        return mat;
    }
    
    /**
     * The overloaded generateNNMatrix method creates an nXn matrix with randomly
     * assigned real numbers in the range [0,1].
     * @param n     the dimensions of the nXn matrix
     * @return a populated nXn matrix
     */
    private static double[][] generateNNMatrix(int n){
        double[][] mat = new double[n][n];
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                mat[i][j] = Math.random();
            }
        }
        return mat;
    }
    
    /**
     * The overloaded generateNNMatrix method creates an nXn matrix with randomly
     * assigned values less than the specified maximum value
     * @param n     the dimensions of the nXn matrix
     * @param range the maximum value to be allowed in the matrix
     * @return a populated nXn matrix
     */
    private static double[][] generateNNMatrix(int n, int range){
        double[][] mat = new double[n][n];
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                mat[i][j] = random.nextInt(range);
            }
        }
        return mat;
    }
    
    /**
     * The matrixMultiply method multiplies two nXn matrices of the same size using the
     * conventional method.
     * @param firstMat  the first nXn matrix
     * @param secondMat the second nXn matrix
     * @return the product matrix
     */
    private static double[][] matrixMultiply(double[][] firstMat, double[][] secondMat){
        if(firstMat.length != secondMat.length){
            throw new RuntimeException("Matrices not equal dimensions");
        }
	double[][] productMat = new double[firstMat.length][firstMat.length];
	for(int i = 0; i < firstMat.length; i++){
            for(int j = 0; j < firstMat.length; j++){
                for(int k = 0; k < firstMat.length; k++){
                    productMat[i][j] += firstMat[i][k]*secondMat[k][j];
                    convOperations += 2; // one multiplication and one addition operation
                }
            }
	}
	return productMat;
    }
    
    /**
     * The strassen method multiplies two nXn matrices using strassen's method
     * @param firstMat  the first nXn matrix
     * @param secondMat the second nXn matrix
     * @return the product matrix
     */
    private static double[][] strassen(double[][] firstMat, double[][] secondMat){
        if(firstMat.length != secondMat.length){
            throw new RuntimeException("Matrices not equal dimensions");
        }
        if(firstMat.length == 1){
            firstMat[0][0] = firstMat[0][0]*secondMat[0][0];
            strassOperations++;
            return firstMat;
        }
        int n = firstMat.length;
        double[][] productMat = new double[n][n];
        
        // Strassen's method part 1
        double[][] a11 = getSubmatrix(firstMat,1);
        double[][] a12 = getSubmatrix(firstMat,2);
        double[][] a21 = getSubmatrix(firstMat,3);
        double[][] a22 = getSubmatrix(firstMat,4);
        
        double[][] b11 = getSubmatrix(secondMat,1);
        double[][] b12 = getSubmatrix(secondMat,2);
        double[][] b21 = getSubmatrix(secondMat,3);
        double[][] b22 = getSubmatrix(secondMat,4);
        
        /* the number of operations for each addition or subtraction of matrices can be
         * predetermined since each a and b matrix have the same dimensions. when adding
         * one matrix to another, each element is added to its counterpart, thus the cost
         * of adding one matrix to another is equal to the number of elements it has
         */
        int operationCost = a11.length*a11.length;
        
        // Strassen's algorithm part 2
        double[][] s1 = subtractMatrices(b12,b22);
        strassOperations += operationCost;
        double[][] s2 = addMatrices(a11,a12);
        strassOperations += operationCost;
        double[][] s3 = addMatrices(a21,a22);
        strassOperations += operationCost;
        double[][] s4 = subtractMatrices(b21,b11);
        strassOperations += operationCost;
        double[][] s5 = addMatrices(a11,a22);
        strassOperations += operationCost;
        double[][] s6 = addMatrices(b11,b22);
        strassOperations += operationCost;
        double[][] s7 = subtractMatrices(a12,a22);
        strassOperations += operationCost;
        double[][] s8 = addMatrices(b21,b22);
        strassOperations += operationCost;
        double[][] s9 = subtractMatrices(a11,a21);
        strassOperations += operationCost;
        double[][] s10 = addMatrices(b11,b12);
        strassOperations += operationCost;
        
        // Strassen's algorithm part 3 (7 multiplications)
        double[][] p1 = strassen(a11,s1);
        double[][] p2 = strassen(s2,b22);
        double[][] p3 = strassen(s3,b11);
        double[][] p4 = strassen(a22,s4);
        double[][] p5 = strassen(s5,s6);
        double[][] p6 = strassen(s7,s8);
        double[][] p7 = strassen(s9,s10);
        
        // the operation cost can again be precalculated
        operationCost = p1.length*p1.length;
        
        // Strassen's algorithm part 4
        double[][] c11 = addMatrices(subtractMatrices(addMatrices(p5,p4),p2),p6);
        strassOperations += operationCost + operationCost + operationCost;
        double[][] c12 = addMatrices(p1,p2);
        strassOperations += operationCost;
        double[][] c21 = addMatrices(p3,p4);
        strassOperations += operationCost;
        double[][] c22 = subtractMatrices(subtractMatrices(addMatrices(p5,p1),p3),p7);
        strassOperations += operationCost + operationCost + operationCost;
        
        // build the product matrix with the four submatrices
        for(int i = 0; i < c11.length; i++){
            System.arraycopy(c11[i], 0, productMat[i], 0, c11.length);
        }
        for(int i = 0; i < c12.length; i++){
            System.arraycopy(c12[i], 0, productMat[i], n/2, c12.length);
        }
        for(int i = 0; i < c21.length; i++){
            System.arraycopy(c21[i], 0, productMat[i+n/2], 0, c21.length);
        }
        for(int i = 0; i < c22.length; i++){
            System.arraycopy(c22[i], 0, productMat[i+n/2], n/2, c22.length);
        }
        return productMat;
    }
    
    /**
     * The overloaded strassen method multiplies two nXn matrices using strassen's method
     * but also specifies a crossover point where the method should multiply the matrices
     * conventionally.
     * @param firstMat  the first nXn matrix
     * @param secondMat the second nXn matrix
     * @param crossOver the n at which the method should just used conventional
     *                  multiplication
     * @return the product matrix
     */
    private static double[][] strassen(double[][] firstMat,
            double[][] secondMat, int crossOver){
        if(firstMat.length != secondMat.length){
            throw new RuntimeException("Matrices not equal dimensions");
        }
        if(firstMat.length <= crossOver){
            return divideAndConquer(firstMat, secondMat);
        }
        int n = firstMat.length;
        double[][] productMat = new double[n][n];
        
        // Strassen's method part 1
        double[][] a11 = getSubmatrix(firstMat,1);
        double[][] a12 = getSubmatrix(firstMat,2);
        double[][] a21 = getSubmatrix(firstMat,3);
        double[][] a22 = getSubmatrix(firstMat,4);
        
        double[][] b11 = getSubmatrix(secondMat,1);
        double[][] b12 = getSubmatrix(secondMat,2);
        double[][] b21 = getSubmatrix(secondMat,3);
        double[][] b22 = getSubmatrix(secondMat,4);
        
        /* the number of operations for each addition or subtraction of matrices can be
         * predetermined since each a and b matrix have the same dimensions. when adding
         * one matrix to another, each element is added to its counterpart, thus the cost
         * of adding one matrix to another is equal to the number of elements it has
         */
        int operationCost = a11.length*a11.length;
        
        // Strassen's algorithm part 2
        double[][] s1 = subtractMatrices(b12,b22);
        strassOperations += operationCost;
        double[][] s2 = addMatrices(a11,a12);
        strassOperations += operationCost;
        double[][] s3 = addMatrices(a21,a22);
        strassOperations += operationCost;
        double[][] s4 = subtractMatrices(b21,b11);
        strassOperations += operationCost;
        double[][] s5 = addMatrices(a11,a22);
        strassOperations += operationCost;
        double[][] s6 = addMatrices(b11,b22);
        strassOperations += operationCost;
        double[][] s7 = subtractMatrices(a12,a22);
        strassOperations += operationCost;
        double[][] s8 = addMatrices(b21,b22);
        strassOperations += operationCost;
        double[][] s9 = subtractMatrices(a11,a21);
        strassOperations += operationCost;
        double[][] s10 = addMatrices(b11,b12);
        strassOperations += operationCost;
        
        // Strassen's algorithm part 3
        double[][] p1 = strassen(a11,s1);
        double[][] p2 = strassen(s2,b22);
        double[][] p3 = strassen(s3,b11);
        double[][] p4 = strassen(a22,s4);
        double[][] p5 = strassen(s5,s6);
        double[][] p6 = strassen(s7,s8);
        double[][] p7 = strassen(s9,s10);
        
        // the operation cost can again be precalculated
        operationCost = p1.length*p1.length;
        
        // Strassen's algorithm part 4
        double[][] c11 = addMatrices(subtractMatrices(addMatrices(p5,p4),p2),p6);
        strassOperations += operationCost + operationCost + operationCost;
        double[][] c12 = addMatrices(p1,p2);
        strassOperations += operationCost;
        double[][] c21 = addMatrices(p3,p4);
        strassOperations += operationCost;
        double[][] c22 = subtractMatrices(subtractMatrices(addMatrices(p5,p1),p3),p7);
        strassOperations += operationCost + operationCost + operationCost;
        
        // build the product matrix with the four submatrices
        for(int i = 0; i < c11.length; i++){
            System.arraycopy(c11[i], 0, productMat[i], 0, c11.length);
        }
        for(int i = 0; i < c12.length; i++){
            System.arraycopy(c12[i], 0, productMat[i], n/2, c12.length);
        }
        for(int i = 0; i < c21.length; i++){
            System.arraycopy(c21[i], 0, productMat[i+n/2], 0, c21.length);
        }
        for(int i = 0; i < c22.length; i++){
            System.arraycopy(c22[i], 0, productMat[i+n/2], n/2, c22.length);
        }
        return productMat;
    }
    
    /**
     * The divide and conquer method multiplies two matrices using the divide and conquer
     * method.
     * @param firstMat  the first nXn matrix
     * @param secondMat the second nXn matrix
     * @return the product of the two matrices
     */
    private static double[][] divideAndConquer(double[][] firstMat, double[][] secondMat){
        if(firstMat.length != secondMat.length){
            throw new RuntimeException("Matrices not equal dimensions");
        }

        int n = firstMat.length;
        if(n == 1){
            firstMat[0][0] = firstMat[0][0]*secondMat[0][0];
            convOperations++;
            return firstMat;
        }
        else{
            double[][] productMat = new double[n][n];
            
            int fee = (firstMat.length/2)*(firstMat.length/2);
            
            convOperations += fee;
            double[][] topLeft = addMatrices(
                    divideAndConquer(getSubmatrix(firstMat, 1),
                    getSubmatrix(secondMat, 1)),
                    divideAndConquer(getSubmatrix(firstMat, 2),
                    getSubmatrix(secondMat, 3)));
            
            convOperations += fee;
            double[][] topRight = addMatrices(
                    divideAndConquer(getSubmatrix(firstMat, 1),
                    getSubmatrix(secondMat, 2)),
                    divideAndConquer(getSubmatrix(firstMat, 2),
                    getSubmatrix(secondMat, 4)));
            
            convOperations += fee;
            double[][] bottomLeft = addMatrices(
                    divideAndConquer(getSubmatrix(firstMat, 3),
                    getSubmatrix(secondMat, 1)),
                    divideAndConquer(getSubmatrix(firstMat, 4),
                    getSubmatrix(secondMat, 3)));
            
            convOperations += fee;
            double[][] bottomRight = addMatrices(
                    divideAndConquer(getSubmatrix(firstMat, 3),
                    getSubmatrix(secondMat, 2)),
                    divideAndConquer(getSubmatrix(firstMat, 4),
                    getSubmatrix(secondMat, 4)));
            
            // build the product matrix with the four submatrices
            for(int i = 0; i < topLeft.length; i++){
                System.arraycopy(topLeft[i],
                        0, productMat[i], 0, topLeft.length);
            }
            for(int i = 0; i < topRight.length; i++){
                System.arraycopy(topRight[i],
                        0, productMat[i], n/2, topRight.length);
            }
            for(int i = 0; i < bottomLeft.length; i++){
                System.arraycopy(bottomLeft[i],
                        0, productMat[i+n/2], 0, bottomLeft.length);
            }
            for(int i = 0; i < bottomRight.length; i++){
                System.arraycopy(bottomRight[i],
                        0,productMat[i+n/2], n/2, bottomRight.length);
            }
            
            return productMat;
        }   
    }
    
    /**
     * The getSubmatrix method returns a submatrix of another matrix of size (n/2)X(n/2).
     * @param oldMat    the matrix from which the submatrix will come
     * @param quadrant  dictates which submatrix to extract from the original matrix,
     *                  where quadrant 1 is the top left quarter of the matrix, 2 is
     *                  the top right, 3 is the bottom left, and 4 is the bottom right.
     * @return an (n/2)X(n/2) submatrix
     */
    private static double[][] getSubmatrix(double[][] oldMat, int quadrant){
        double[][] newMat = new double[oldMat.length/2][oldMat.length/2];
        
        int colBegin = quadrant == 1||quadrant == 3 ? 0 : oldMat.length/2;
        int colEnd = colBegin + oldMat.length/2;
        int rowBegin = quadrant == 1||quadrant == 2 ? 0 : oldMat.length/2;
        int rowEnd = rowBegin + oldMat.length/2;
        
        int colIndex = 0;
        int rowIndex = 0;

        for(int i = rowBegin; i < rowEnd; i++){
            for(int j = colBegin; j < colEnd; j++){
                newMat[rowIndex][colIndex] = oldMat[i][j];
                colIndex++;
            }
            colIndex = 0;
            rowIndex++;
        }
        return newMat;
    }
    
    /**
     * The addMatrices method adds two matrices together by conventional means.
     * @param firstMat  the first nXn matrix
     * @param secondMat the second nXn matrix
     * @return the sum matrix
     */
    private static double[][] addMatrices(double[][] firstMat, double[][] secondMat){
        double[][] newMat = new double[firstMat.length][firstMat.length];
        for(int i = 0; i < firstMat.length; i++){
            for(int j = 0; j < firstMat.length; j++){
                newMat[i][j] = firstMat[i][j] + secondMat[i][j];
            }
        }
        return newMat;
    }
    
    /**
     * The subtractMatrices method subtracts one matrix from another by conventional
     * means.
     * @param firstMat  the nXn matrix to be subtracted from
     * @param secondMat the nXn matrix that will be subtracted from firstMat
     * @return the difference matrix
     */
    private static double[][] subtractMatrices(double[][] firstMat, double[][] secondMat){
        double[][] newMat = new double[firstMat.length][firstMat.length];
        for(int i = 0; i < firstMat.length; i++){
            for(int j = 0; j < firstMat.length; j++){
                newMat[i][j] = firstMat[i][j] - secondMat[i][j];
            }
        }
        return newMat;
    }
    
    /**
     * Prints a given matrix.
     * @param mat the matrix to be printed
     */
    private static void printMat(double[][] mat){
        System.out.print("[");
        for(int i = 0; i < mat.length; i ++){
            System.out.print("[");
            for(int j = 0; j < mat[i].length; j++){
                System.out.print(mat[i][j] + ", ");
            }
            System.out.print("]" + "\n");
        }
    }
    
    /**
     * The equals method determines if two matrices are identical.
     * @param firstMat  the first nXn matrix
     * @param secondMat the second nXn matrix
     * @return boolean indicating whether or not the matrices are equal
     */
    private static boolean equals(double[][] firstMat, double[][] secondMat){
        for(int i = 0; i < firstMat.length; i++){
            for(int j = 0; j < firstMat.length; j++){
                /* occassionally the different multiplication methods produce negligibly
                 * different results when using any real number
                 */
                if(Math.abs(firstMat[i][j] - secondMat[i][j]) > 0.00000000001){
                    System.out.println(firstMat[i][j] + " and " + secondMat[i][j]);
                    return false;
                }
            }
        }
        return true;
    }
    
    /**
     * The overloaded equals method determines if three matrices are identical.
     * @param firstMat  the first nXn matrix
     * @param secondMat the second nXn matrix
     * @param thirdMat  the third nXn matrix
     * @return boolean indicating whether or not the matrices are equal
     */
    private static boolean equals(double[][] firstMat, double[][] secondMat,
            double[][] thirdMat){
        for(int i = 0; i < firstMat.length; i++){
            for(int j = 0; j < firstMat.length; j++){
                /* occassionally the different multiplication methods produce negligibly
                 * different results when using any real number
                 */
                if(Math.abs(firstMat[i][j] - secondMat[i][j]) > 0.000000000001 ||
                        Math.abs(firstMat[i][j] - thirdMat[i][j]) > 0.000000000001 ||
                        Math.abs(secondMat[i][j] - thirdMat[i][j]) > 0.000000000001){
                    return false;
                }
            }
        }
        return true;
    }
    
    /**
     * The buildPowTwoMat method pads a given matrix with zeroes so that its dimensions
     * are powers of two.
     * @param mat   the matrix to be padded with zeroes
     * @return an nXn matrix where n is a power of two
     */
    private static double[][] buildPowTwoMat(double[][] mat){
        int powTwoLength = closestPowTwo(mat);
        if(powTwoLength == 1){
            return mat;
        }
        double[][] powTwo = new double[powTwoLength][powTwoLength];
        for(int i = 0; i < mat.length; i++){
            System.arraycopy(mat[i], 0, powTwo[i], 0, mat[0].length);
        }
        return powTwo;
    }
    
    /**
     * The closestPowTwo method determines the nearest power of two to the dimensions of
     * a given nXn matrix.
     * @param firstMat  the matrix to be evaluated
     * @return the nearest power of two that is greater than n
     */
    private static int closestPowTwo(double[][] firstMat){
        int minValue = 1;
        for(int i = 0; i < firstMat.length; i++){
            int value = (int) Math.pow(2,i);
            if(Math.abs(firstMat.length - value) < firstMat.length - minValue
                    && value > firstMat.length){
                minValue = value;
            }
        }
        return minValue;
    }
}