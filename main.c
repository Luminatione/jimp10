#include "points.h"
#include "splines.h"
#include "makespl.h"

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

char* usage =
"Usage: %s -s spline-file [-p points-file] [ -g gnuplot-file [-f from_x -t to_x -n n_points ] ]\n"
"            if points-file is given then\n"
"               reads discrete 2D points from points-file\n"
"               writes spline approximation to spline-file\n"
"               - number of points should be >= 4\n"
"            else (points-file not given)\n"
"               reads spline from spline-file\n"
"            endfi\n"
"            if gnuplot-file is given then\n"
"               makes table of n_points within <from_x,to_x> range\n"
"               - from_x defaults to x-coordinate of the first point in points-file,\n"
"               - to_x defaults to x-coordinate of the last point\n"
"               - n_points defaults to 100\n"
"               - n_points must be > 1\n"
"            endif\n";

int
main(int argc, char** argv)
{
	int opt;
	char* inputFileName = NULL;
	char* splineFileName = NULL;
	char* plotFileName = NULL;
	double fromX = 0;
	double toX = 0;
	int n = 100;
	char* executableName = argv[0];

	points_t points;
	spline_t line;

	points.n = 0;
	line.n = 0;

	/* process options, save user choices */
	while ((opt = getopt(argc, argv, "p:s:g:f:t:n:")) != -1) {
		switch (opt) {
		case 'p':
			inputFileName = optarg;
			break;
		case 's':
			splineFileName = optarg;
			break;
		case 'g':
			plotFileName = optarg;
			break;
		case 'f':
			fromX = atof(optarg);
			break;
		case 't':
			toX = atof(optarg);
			break;
		case 'n':
			n = atoi(optarg);
			break;
		default:                   /* '?' */
			fprintf(stderr, usage, executableName);
			exit(EXIT_FAILURE);
		}
	}
	if (optind < argc) {
		fprintf(stderr, "\nBad parameters!\n");
		for (; optind < argc; optind++)
			fprintf(stderr, "\t\"%s\"\n", argv[optind]);
		fprintf(stderr, "\n");
		fprintf(stderr, usage, executableName);
		exit(EXIT_FAILURE);
	}

	/* if points-file was given, then read points, generate spline, save it to file */
	if (inputFileName != NULL) {
		FILE* inputFile = fopen(inputFileName, "r");
		if (inputFile == NULL) {
			fprintf(stderr, "%s: can not read points file: %s\n\n", argv[0], inputFileName);
			exit(EXIT_FAILURE);
		}

		if (read_pts_failed(inputFile, &points)) {
			fprintf(stderr, "%s: bad contents of points file: %s\n\n", argv[0],
				inputFileName);
			exit(EXIT_FAILURE);
		}
		fclose(inputFile);

		FILE* outputFile = fopen(splineFileName, "w");
		if (outputFile == NULL) {
			fprintf(stderr, "%s: can not write spline file: %s\n\n", argv[0], splineFileName);
			exit(EXIT_FAILURE);
		}

		make_spl(&points, &line);

		if (line.n > 0)
			write_spl(&line, outputFile);

		fclose(outputFile);
	}
	else if (splineFileName != NULL) {  /* if point-file was NOT given, try to read splines from a file */
		FILE* splinesFile = fopen(splineFileName, "r");
		if (splinesFile == NULL) {
			fprintf(stderr, "%s: can not read spline file: %s\n\n", argv[0], inputFileName);
			exit(EXIT_FAILURE);
		}
		if (read_spl(splinesFile, &line)) {
			fprintf(stderr, "%s: bad contents of spline file: %s\n\n", argv[0],
				inputFileName);
			exit(EXIT_FAILURE);
		}
	}
	else { /* ponts were not given nor spline was given -> it is an error */
		fprintf(stderr, usage, argv[0]);
		exit(EXIT_FAILURE);
	}

	if (line.n < 1) { /* check if there is a valid spline */
		fprintf(stderr, "%s: bad spline: n=%d\n\n", argv[0], line.n);
		exit(EXIT_FAILURE);
	}

	/* check if plot was requested and generate it if yes */
	if (plotFileName != NULL && n > 1) {
		FILE* plotFile = fopen(plotFileName, "w");

		if (fromX == 0.0 && toX == 0.0) { /* calculate plot range if it was not specified */
			if (points.n > 1) {
				fromX = points.x[0];
				toX = points.x[points.n - 1];
			}
			else if (line.n > 1) {
				fromX = line.x[0];
				toX = line.x[line.n - 1];
			}
			else {
				fromX = 0;
				toX = 1;
			}
		}
		double dx = (toX - fromX) / (n - 1);

		if (plotFile == NULL) {
			fprintf(stderr, "%s: can not write gnuplot file: %s\n\n", argv[0],
				plotFileName);
			exit(EXIT_FAILURE);
		}
		int i;
		for (i = 0; i < n; i++)
			fprintf(plotFile, "%g %g\n", fromX + i * dx,
				value_spl(&line, fromX + i * dx));

		fclose(plotFile);
	}

	return 0;
}
