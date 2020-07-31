///////////////////////////////////////////////////////////
// James Schubach, 29743338, jsch0026@student.monash.edu //
///////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <mpi.h>

#define RESULT 0
#define ASK_FOR_JOB 1
#define JOB_DATA 2
#define STOP 3
#define iXmax 8000 // default
#define iYmax 8000 // default

static MPI_File fp;
FILE * fpp;

struct tuple {
	int a;
	int b;
};





void compute(int y, int header) {

	MPI_Status stat3;
	/* world ( double) coordinate = parameter plane*/
	double Cx, Cy;
	const double CxMin = -2.5;
	const double CxMax = 1.5;
	const double CyMin = -2.0;
	const double CyMax = 2.0;

	/* */
	double PixelWidth = (CxMax - CxMin)/iXmax;
	double PixelHeight = (CyMax - CyMin)/iYmax;

	/* color component ( R or G or B) is coded from 0 to 255 */
	/* it is 24 bit color RGB file */
	const int MaxColorComponentValue = 255; 
	// FILE * fp;
	// char *filename = "Mandelbrot.ppm";
	// char *comment = "# ";	/* comment should start with # */

	// RGB color array
	static unsigned char color[iXmax*3];

	/* Z = Zx + Zy*i;	Z0 = 0 */
	double Zx, Zy;
	double Zx2, Zy2; /* Zx2 = Zx*Zx;  Zy2 = Zy*Zy  */
	/*  */
	int Iteration;
	const int IterationMax = 1000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;

	int colorsize; 


	Cy = CyMin + (y * PixelHeight);

	if (fabs(Cy) < (PixelHeight / 2))
	{
		Cy = 0.0; /* Main antenna */
	}

	for(int x = 0; x < (iXmax); x++)
	{

		Cx = CxMin + (x * PixelWidth);
		/* initial value of orbit = critical point Z= 0 */
		Zx = 0.0;
		Zy = 0.0;
		Zx2 = Zx * Zx;
		Zy2 = Zy * Zy;

		/* */
		for(Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++)
		{
			Zy = (2 * Zx * Zy) + Cy;
			Zx = Zx2 - Zy2 + Cx;
			Zx2 = Zx * Zx;
			Zy2 = Zy * Zy;
		};

		/* compute  pixel color (24 bit = 3 bytes) */
		if (Iteration == IterationMax)
		{
			// Point within the set. Mark it as black
			color[x*3+0] = 0;
			color[x*3+1] = 0;
			color[x*3+2] = 0;
		}
		else 
		{
			// Point outside the set. Mark it as white
			double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0);
			if (c < 1)
			{
				color[x*3+0] = 0;
				color[x*3+1] = 0;
				color[x*3+2] = 255*c;
			}
			else if (c < 2)
			{
				color[x*3+0] = 0;
				color[x*3+1] = 255*(c-1);
				color[x*3+2] = 255;
			}
			else
			{
				color[x*3+0] = 255*(c-2);
				color[x*3+1] = 255;
				color[x*3+2] = 255;
			}
		}	
	}

	MPI_File_seek(fp, header+3*y*iYmax, MPI_SEEK_SET);
	MPI_File_write(fp, &color, iXmax*3, MPI_UNSIGNED_CHAR, &stat3);

}



/* Master */
void master_io(int header) {
	MPI_Status stat, stat2;
	MPI_Request request;
	int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int z = 0;
	int *job = NULL;
	job = (int*)malloc(1 * sizeof(int));
	int buff;
	int *result_data = NULL;
	result_data = (int*)malloc(1 * sizeof(int));
	int one = 1;
	int slaves = size-1;
	int *jobs = NULL;
	jobs = malloc((iYmax * sizeof(int)));
	int count = 0;


	printf("Initializing\n");
	for (int y = 0; y < iYmax; y++) {
			jobs[count] = y;
			count++;
	}
	printf("Initializing complete\n");
	printf("Computing Mandelbrot Set. Please wait...\n");
	float start = MPI_Wtime();


	while (slaves > 1 || z < count) {
		if (slaves == size-1 && size == 2) {
			compute(jobs[z], header);
			z++;
		}
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
		// Wait for any incomming message
		int slave_r = stat.MPI_SOURCE;
		if (stat.MPI_TAG == ASK_FOR_JOB) {
			slaves--;
			MPI_Irecv(&buff, 1, MPI_INT, slave_r, ASK_FOR_JOB, MPI_COMM_WORLD, &request);
			if (z <= count) {
				/* PAck data of job into msg_buff */
				job[0] = jobs[z];
				MPI_Isend(job, 1, MPI_INT, slave_r, JOB_DATA, MPI_COMM_WORLD, &request);
				//job = (int*)malloc(3 * sizeof(int));
				z++;
				slaves++;
			}
			else {
				for (int q = 1; q<size; q++) {
					MPI_Send(&one, 1, MPI_INT, q, STOP, MPI_COMM_WORLD);
				}
			}
		}
	}
	float end = MPI_Wtime();
	printf("Time: %f\n", end-start);
	for (int q = 1; q<size; q++) {
		MPI_Isend(&one, 1, MPI_INT, q, STOP, MPI_COMM_WORLD, &request);
	}
	
}


void slave_io(int header)
{

	MPI_Status stat, stat2;
	MPI_Request request;
	int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int stopp = 0;
	int *job = NULL;
	job = (int*)malloc(1 * sizeof(int));
	int *buff = NULL;
	buff = (int*)malloc(1*sizeof(int));
	int one = 1;

	MPI_Isend(&one, 1, MPI_INT, 0, ASK_FOR_JOB, MPI_COMM_WORLD, &request);
	do {
		
		// Here we send a message to the master asking for a job
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
		if (stat.MPI_TAG == JOB_DATA) {
            // Retrieve job data from master into msg_buffer
            MPI_Recv(job, 1, MPI_INT, 0, JOB_DATA, MPI_COMM_WORLD, &stat2);
			/* compute and write image data bytes to the file */
			compute(job[0], header);
            // send result to master
            MPI_Isend(&one, 1, MPI_INT, 0, ASK_FOR_JOB, MPI_COMM_WORLD, &request);
			
        } 
		else {
            // We got a stop message we have to retrieve it by using MPI_Recv
            // But we can ignore the data from the MPI_Recv call
            MPI_Recv(&buff, 1, MPI_INT , 0, STOP, MPI_COMM_WORLD, &stat2);
            stopp = 1;
        }

	} while (stopp == 0);
}

int main(int argc, char **argv) {

    int rank, size;

	MPI_Status stat;
	char *filename = "Mandelbrot.ppm";
	char *comment = "# ";	/* comment should start with # */
	const int MaxColorComponentValue = 255;
	
	int bufferlen = (sizeof(" ")*4) + sizeof("P6") + strlen(comment)*sizeof(char) + sizeof(iXmax)*4 + sizeof(iYmax)*4 + sizeof(MaxColorComponentValue)*3;

    MPI_Init(&argc , & argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp); 

	if (size >= 2) {
		if (rank == 0) {
			fpp = fopen(filename, "wb");
			fprintf(fpp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);
			fclose(fpp);
			master_io(bufferlen);
		} 
		else {
			// printf("Calling slave");
			slave_io(bufferlen);
		}
	}
	else {
		fpp = fopen(filename, "wb");
		fprintf(fpp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);
		fclose(fpp);
		printf("Initializing\n");
		printf("Initializing complete\n");
		printf("Computing Mandelbrot Set. Please wait...\n");
		float start = MPI_Wtime();
		for (int y = 0; y < iYmax; y++) {
			compute(y, bufferlen);
		}
		float end = MPI_Wtime();
		printf("Time: %f\n", end-start);
	}
	MPI_File_close(&fp); 
	printf("Finished %d\n", rank);
    MPI_Finalize();
	return(0);
}