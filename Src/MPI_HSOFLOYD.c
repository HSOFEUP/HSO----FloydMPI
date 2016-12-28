/*
 ============================================================================
 Name        : MPI_HSOFLOYD_OK.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "/usr/include/openmpi-x86_64/mpi.h"
//#include "mpi.h"

#define TAG 0
#define ERROR -1
#define ROOT 0
#define INF -1

/*
 * Macro: to compute displacement in  the vector
 */
#define GET_MATRIX_VAL(COL,ROW, WIDTH)  (COL) + (WIDTH) * (ROW)

/*
 * Prototypes
 */
int * allocate1D(int nRows, int nCols);
int processArguments(int nArgs, char *argv[], char * inName, char *outName);
void free2DMatrix (int ** matrix, int nRows);
int ** allocate2DMatrix (int nRows, int nCols);
int saveMatrixToFile(char *fileName, int ** matrix, int nRows);
int readFileToMatrixTest(char* fileName, int *nRows, int ***matAdd);
int replaceZerosByInf(int **matrix, int nRows);
void printMatrix(int **matrix, int nRows);
int checkConditions(int P, int Q, int nRows);
int replaceInfByZeros(int **matrix, int nRows);
void floydWahshal(int nSize, int *aMatrix, int *bMatrix, int *cMatrix);


/*
 * Function: allocate1D()
 * Description: allocates an N*M array vector
 * Inputs:
 * 	 nRows - number of rows
 * 	 nCols - number of columns
 * 	 inName - pointer of char filename with matrix to be filled
 * 	 inName - pointer of char filename of the result to be filled
 *
 * Returns:
 * 	 vector pointer
 *
 */
int * allocate1D(int nRows, int nCols){
	int * matrix = malloc((nRows * nCols) * sizeof(int));
	return matrix;
}


/*
 * Function: processArguments()
 * Description: parses the input arguments
 * Inputs:
 * 	 nArgs - number of args passed
 * 	 argv[] - vector containg the arguments strings
 * 	 inName - pointer of char filename with matrix to be filled
 * 	 inName - pointer of char filename of the result to be filled
 *
 * Returns:
 * 	 OK(>0) valid arguments, (<0) invalid parameters
 *
 */
int processArguments(int nArgs, char *argv[], char * inName, char *outName){// all process are synchronized at this point

	if (nArgs <3){
		printf(" Number o arguments insuficient");
		return -1;
	}else{
		// lets parse
		printf("Input File name: %s\n", argv[1]);
		printf("Output File name: %s\n", argv[2]);
		sprintf(inName,"%s", argv[1]);
		sprintf(outName,"%s", argv[2]);
		return 1;
	}
}


/*
 * Function: free2DMatrix()
 * Description: frees from heap an 2D matrix
 * Inputs:
 * 	 matrix - matrix pointer
 * 	 nRows - number rows
 *
 *
 * Returns:
 * 	 void
 *
 */
void free2DMatrix (int ** matrix, int nRows) {

	for (int i = 0; i <nRows; i++)
		free ((int *)matrix[i]);
	free ((int **)matrix);
}


/*
 * Function: sallocate2DMatrix()
 * Description: allocates an 2D matrix
 * Inputs:
 * 	 nRows - number rows
 * 	 nCols - number columns
 *
 *
 * Returns:
 * 	 return matrix pointer
 *
 */
int ** allocate2DMatrix (int nRows, int nCols) {
	int ** matrix;

	matrix = (int **) malloc (nRows*sizeof(int *));
	if (!matrix) return (NULL);

	for (int i=0; i < nRows; i++) {
		matrix[i] = (int *) malloc (nCols*sizeof(int));
		if (!matrix[i]) return (NULL);
	}
	return matrix;
}


/*
 * Function: saveMatrixToFile()
 * Description: saves matrix to file in the speciefied format
 * Inputs:
 * 	 filename - file name to be created
 * 	 matrix - pointer to matrix to be saved
 * 	 nRows - number rows fo the matrix
 *
 * Returns:
 * 	 return OK (>0) condition or error condition(<0)
 *
 */
int saveMatrixToFile(char *fileName, int ** matrix, int nRows){
	int i, j;
	FILE *fp;

	if((fp=fopen(fileName, "wb"))==NULL){
		printf("Error creating file.\n");
		return (-1);
	}

	for(i=0; i< nRows; i++){// all process are synchronized at this point
		for(j =0; j< nRows; j++){
			fprintf(fp, "%d\t", matrix[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return 1; // all ok
}



/*
 * Function: readFileToMatrixTest()
 * Description: read files, determines number Rows, and returns the matrix
 * Inputs:
 * 	 filename - filename to be opened
 * 	 nRows - pinter to be filled with the rows in file
 * 	 matAdd - address of the matrix pointer to be filled
 *
 * Returns:
 * 	 return OK (>0) condition or error condition(<0)
 *
 */
int readFileToMatrixTest(char* fileName, int *nRows, int ***matAdd){

	int i, j, n, value, ret;
	char lineString[1000];

	FILE *fp = fopen(fileName,"r");
	if(fp ==NULL){
		printf("Cannot open file\n\r");
		//exit(1);
		return -1;
	}
	// Get the number of lines to allocate the matrix
	fscanf(fp,"%d", &n);
	printf("Number Lines to read: %d\n", n);
	// allocate
	int ** mat = allocate2DMatrix(n,n); // allocate the matrix

	// Parse the files
	char *token;
	const char s[2]=" ";
	int nline=0; // current reading line

	// go to each line, fist line is clean
	while(fgets(lineString,sizeof(lineString),fp)){
		if(nline!=0)
			printf("Line %d is: %s", nline, lineString);
		/* Establish string and get the first token: */
		token = strtok( lineString,s);
		i=0;
		while( token != NULL )
		{
			/* While there are tokens in "string" */
			//printf("%s\n", token);
			//printf("Current line %d\n", nline);
			if(nline <=n && nline>0)// skip first line
				mat[nline-1][i]=atoi(token); // line -1
			//printf("%d\n", atoi(token));

			i++;
			/* Get next token: */
			token = strtok( NULL, s );
		}
		nline++;
	}

	fclose(fp);
	// read input file
	*nRows =n; // Fill n
	// return matrix
	printf("All read\n");
	*matAdd = mat;
	return 1;
}



/*
 * Function: replaceZerosByInf()
 * Description: replaces INF values of non-diagonal by zeros
 * Inputs:
 * 	 P - total number of processing units
 * 	 Q - sqrt of the number of process
 * 	 nRows - number rows of the original matrix
 *
 * Returns:
 * 	 return OK (>0) condition or error condition(<0)
 *
 */
int replaceZerosByInf(int **matrix, int nRows){
	int i, j;
	for(i=0; i < nRows; i++) {
		for(j=0; j < nRows; j++) {
			if(matrix[i][j] == 0  && i!= j) // for the own state we assume zero,
				matrix[i][j] = INF;
		}
	}
	return 1; // all ok
}


/*
 * Function: printMatrix()
 * Description: print the contents of matrix to console
 * Inputs:
 * 	 matrix - pointer to matrix
 * 	 nRows - number rows of the original matrix
 *
 * Returns:
 * 	 void
 *
 */
void printMatrix(int **matrix, int nRows) {
	int i, j;

	for(i = 0; i < nRows; i++) {
		for(j = 0; j < nRows; j++) {
			if(matrix[i][j] <0)
				printf("%7s", "INF");
			else
				printf("%7d ", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}


/*
 * Function: checkConditions()
 * Description: check if fox algorithm can be run (square process)
 * Inputs:
 * 	 P - total number of processing units
 * 	 Q - sqrt of the number of process
 * 	 nRows - number rows of the original matrix
 *
 * Returns:
 * 	 return OK (>0) condition or error condition(<0)
 *
 */
int checkConditions(int P, int Q, int nRows){
	if(nRows % Q !=0 || P != Q*Q){ // casted to int if decimal cannot give save value
		return -1;
	}else{
		return 1; // all ok
	}
}


/*
 * Function: replaceInfByZeros()
 * Description: replaces zeros on non-diagonal by INF
 * Inputs:
 * 	 matrix- pointer to the matrix to be replaced
 * 	 nRows - number rows of the matrix
 *
 * Returns:
 * 	 OK condition all went well
 */
int replaceInfByZeros(int **matrix, int nRows){
	int i, j;

	for(i = 0; i < nRows; i++) {
		for(j=0; j < nRows; j++) {
			if(matrix[i][j] < 0)
				matrix[i][j] = 0;
		}
	}
	return 1; // all ok
}


/*
 * Function: floydWharshal()
 * Description: computer the minimum adejancy path in the matrix
 * Inputs:
 * 	 nSize - submatrix size
 * 	 aMatrix - pointer to A
 * 	 bMatrix - pointer to B
 * 	 cMatrix - pointer to result matrix C
 *
 * Returns:
 * 	 void
 */
void floydWahshal(int nSize, int *aMatrix, int *bMatrix, int *cMatrix) {
	int aVal, bVal, cVal; // hold values

	for (int i = 0; i < nSize; i++) {
		for (int j = 0; j < nSize; j++) {
			for (int k = 0; k < nSize; k++) {
				// we waste three CPU cycles on this copy
				aVal = aMatrix[k + i * nSize];
				bVal = bMatrix[j + k * nSize];
				cVal = cMatrix[j + i * nSize];

				if (aVal != INF && bVal != INF) {
					if (cVal == INF || cVal > aVal + bVal)
						cMatrix[j + i * nSize] = aVal + bVal;
				}
			}
		}
	}
}


// MAIN
int main (int argc, char * argv[]) {
	int P, Q, rank, nRows, nCols;
	int subMatrixSize;
	int i, j;

	// Array pointers
	int **matrix; // 2 D
	// array vectors
	int *localMatrix;
	int *rowMatrix;
	int *colMatrix;
	int *resMatrix;
	int *solMatrix;

	int srcCol, destCol;


	// MPI vars
	MPI_Status status;
	MPI_Comm rowsComm;
	MPI_Comm colsComm;
	MPI_Comm gridComm;

	// Aux vars
	int errCode; // var to handle error codes
	double startTime, endTime;

	char inName[10];
	char outName[10];

	// Lets Start things
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &P); // Grab number Processing units
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Grab the rank

	// Get the squares of value// all process know this
	Q = sqrt(P);


	// SETUP things
	// Rank zero reads the file and check conditions and devides chunc and anounces
	if (rank == ROOT) {
		// parse the arguments
		errCode = processArguments(argc, argv, inName, outName);
		if(errCode <0){
			subMatrixSize = ERROR;

			// broadcast to all process this error for then terminate also
			MPI_Bcast(&subMatrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Finalize(); // I also terminate
			printf ("Invalid arguments\n");
			return 0;
		}
		// lets read the file
		errCode= readFileToMatrixTest(inName, &nRows, &matrix);
		//matrix = readFileToMatrix(inName, &nRows); // nrows comes populated
		if(errCode<0){
			subMatrixSize = ERROR;

			// broadcast to all process this error for then terminate also
			MPI_Bcast(&subMatrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Finalize(); // I also terminate
			printf ("Cannot read file\n");
			return 0;

		}

		// lest test the conditions
		if(checkConditions(P,Q,nRows) <0){ // we cannot run this
			subMatrixSize = ERROR;
			// broadcast to all process this condition for then terminate also
			MPI_Bcast(&subMatrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Finalize(); // I also terminate
			printf ("The current setup does not allow to run the algorithm\n");
			return 0;
		}
		// Make things square
		nCols=nRows;

		// Allocate the solution to save time
		solMatrix=allocate1D(nRows, nRows);
		if(solMatrix==NULL){
			subMatrixSize = ERROR;
			// broadcast to all process this condition for then terminate also
			MPI_Bcast(&subMatrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Finalize(); // I also terminate
			printf("Error allocating Chunk\n");
			return 0;

		}


		// Lets handle the matrix correctly
		errCode = replaceZerosByInf(matrix, nRows);
		// Test if all went well
		if(errCode <0){
			printf(" Something went wrong\n");
			//subMatrixSize=ERROR;
			// announce all process of error
			MPI_Bcast(&subMatrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Finalize();
			return 0;
		}

		// Lets show the matrix
		printMatrix(matrix,nRows);

		// From here all went well, lets proceed

		// Start tick
		startTime = MPI_Wtime();

		// Start slicing matrix and sent the chunks
		subMatrixSize = nRows / Q;

		// Inform the process of their size chunks
		MPI_Bcast(&subMatrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);


		int *tempMatrix; // Local aux chunk to be used to send the chunks to remainder process
		for (i = 0; i < P; i++) {
			if (i == 0)// My local copy
				localMatrix = allocate1D(subMatrixSize,subMatrixSize);
			else // For the others
				tempMatrix = allocate1D(subMatrixSize,subMatrixSize);

			int subMatrixPosX = i / Q;
			int subMatrixPosY = i % Q;
			for (int rows = 0; rows < subMatrixSize; rows++) {
				for (int cols = 0; cols < subMatrixSize; cols++) {
					if (i == 0) // if im ROOT get my local copy
						localMatrix[GET_MATRIX_VAL(cols, rows, subMatrixSize)] = matrix[rows + (subMatrixPosX * subMatrixSize)][cols + (subMatrixPosY * subMatrixSize)];
					else // For the other process, get an temp copy and send to the process
						tempMatrix[GET_MATRIX_VAL(cols,rows,subMatrixSize)] = matrix[rows + (subMatrixPosX * subMatrixSize)][cols + (subMatrixPosY * subMatrixSize)];

				}
			}
			// Except me, I don't need to send to myself, i identifies the process to be sent
			if (i != 0) MPI_Send(tempMatrix, subMatrixSize*subMatrixSize, MPI_INT, i, TAG, MPI_COMM_WORLD);
		}
	}

	// Now for all remainder process except ROOT and lets receive the matrix chunks, ROOT already has
	// its own copy
	else {
		// Lets receive the announced chunk size
		MPI_Bcast(&subMatrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// if size not valid lets abort
		if (subMatrixSize == ERROR) {
			MPI_Finalize();
			return -1;
		}
		// If all OK each process allocates an local chunk and receive the announced matrix
		localMatrix = allocate1D(subMatrixSize, subMatrixSize);
		if(localMatrix==NULL){
				printf("Error allocating Chunk\n");
				MPI_Finalize();
				return -1;
			}
		// All remainder process receive their chunks to their local matrices rank > 0
		MPI_Recv(localMatrix, subMatrixSize * subMatrixSize, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status);
	}

	// Setup the cartesian matrix for process
	int dims[2], periods[2];
	dims[0] = Q;
	dims[1] = Q;
	periods[0] = 1;
	periods[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &gridComm);

	// Setup Col communicator
	dims[0] = 1;
	dims[1] = 0;
	MPI_Cart_sub(gridComm, dims, &colsComm);

	// Setup row communicator
	dims[0] = 0;
	dims[1] = 1;
	MPI_Cart_sub(gridComm, dims, &rowsComm);
	// Calculate the coordinated in the process grid
	MPI_Cart_coords(gridComm, rank, 2, dims);

	// All process allocate their vector chunks
	rowMatrix = allocate1D(subMatrixSize, subMatrixSize);
	if(rowMatrix==NULL){
		printf("Error allocating Chunk\n");
		MPI_Finalize();
		return 0;
	}

	colMatrix = allocate1D(subMatrixSize, subMatrixSize);
	if(colMatrix==NULL){
			printf("Error allocating Chunk\n");
			MPI_Finalize();
			return 0;
	}

	resMatrix = allocate1D(subMatrixSize, subMatrixSize);
	if(resMatrix==NULL){
			printf("Error allocating Chunk\n");
			MPI_Finalize();
			return 0;
	}
	// Diagonal process send their matrix to the lines, this
	if ((rank / Q == (rank % Q)))// Make local copy for diagonal process in rowMatrix
		memcpy(rowMatrix, localMatrix, (subMatrixSize * subMatrixSize) * sizeof(int));

	// Make initial local copies
	memcpy(colMatrix, localMatrix, (subMatrixSize * subMatrixSize) * sizeof(int));
	memcpy(resMatrix, localMatrix, (subMatrixSize * subMatrixSize) * sizeof(int));


	// Calculate addresses for circular shift ∗/
	srcCol = (dims[0] + 1) % Q;
	destCol = ((dims[0] - 1) + Q) % Q;

	// Stride 2 in two, this is to avoid sending nRows original
	int len = 2*(subMatrixSize*Q);
	for (int df = 2; df < len; df = df<<1) { // Number of time FOX is ivoqued = logN

		// Based in the Article for matrix Multiplication
		for (int stage = 0; stage < Q; stage++) { // Each stage
			int bcast = (dims[0] + stage) % Q; // when o main diagonal

			if (bcast == dims[1]) { // I'm one on the main diagonal
				// Broadcast the local matrix to process in the same line
				MPI_Bcast(localMatrix, subMatrixSize*subMatrixSize, MPI_INT, bcast, rowsComm);
				floydWahshal(subMatrixSize, localMatrix, colMatrix, resMatrix);
			} else {
				//  Receive the broadcasted matrix line
				MPI_Bcast(rowMatrix, subMatrixSize*subMatrixSize, MPI_INT, bcast, rowsComm);
				floydWahshal(subMatrixSize, rowMatrix, colMatrix, resMatrix);
			}
			// Send the colValues value and receives rotated colValues in place the colMatrix to be used in next round using the torus
			MPI_Sendrecv_replace(colMatrix, subMatrixSize*subMatrixSize, MPI_INT, destCol, TAG, srcCol, TAG, colsComm, MPI_STATUS_IGNORE);
		}

		// We need to ensure that all process have finished for the next round to keep the algnment
		MPI_Barrier(MPI_COMM_WORLD);

		// All process ended this stage, prepare next round by saving to localMatrix and colMatrix
		memcpy(localMatrix, resMatrix, (subMatrixSize * subMatrixSize) * sizeof(int));
		// update colMatrix with Local for next round
		memcpy(colMatrix, resMatrix, (subMatrixSize * subMatrixSize) * sizeof(int));
	}


	// All process sincronize at this point// not nedded
	MPI_Barrier(MPI_COMM_WORLD); // can be skipped

	// Now all process make an copy of their partial results to final
	// Gather all partial result in resMatrix to the final solution Matrix
	MPI_Gather(resMatrix, subMatrixSize*subMatrixSize, MPI_INT, solMatrix, subMatrixSize*subMatrixSize, MPI_INT, 0, MPI_COMM_WORLD);

	int fullLength = subMatrixSize*subMatrixSize;
	// Root now build the result matrix
	if(rank == ROOT){
		// The number of runs is the total amount of processing units involved
		for (i = 0; i < P; i++) {
			int subMatrixPosX = i / Q;
			int subMatrixPosY = i % Q;
			int pos = i * fullLength;

			for (int rows = 0; rows < subMatrixSize; rows++) {
				for (int cols = 0; cols < subMatrixSize; cols++) {
					// Convert back to the original 2D form
					matrix[rows + (subMatrixPosX * subMatrixSize)][cols + (subMatrixPosY * subMatrixSize)] = solMatrix[pos + (cols + rows*subMatrixSize)];

				}
			}
		}

		// All done, grab finishing time
		endTime = MPI_Wtime();
		printf("Total time: %f\n", (endTime - startTime));
		// Back to original format
		replaceInfByZeros(matrix,nRows);
		// Print the result matrix
		printMatrix(matrix, nRows);
		// lets save to file
		saveMatrixToFile(outName, matrix, nRows);
		// free matrix
		free2DMatrix(matrix,nRows);
		// Free array
		free(solMatrix);

	}



	// all process free their local matrix
	free(localMatrix);
	free(rowMatrix);
	free(colMatrix);
	free(resMatrix);


	// Lets end things
	MPI_Finalize();
	return 0;
}
