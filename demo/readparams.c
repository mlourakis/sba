/* Loading of camera, 3D point & image projection parameters from disk files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "readparams.h"


#define MAXSTRLEN  1024 /* 1K */

/* get rid of the rest of a line upto \n or EOF */
#define SKIP_LINE(f){                                                       \
char buf[MAXSTRLEN];                                                        \
  while(!feof(f))                                                           \
    if(!fgets(buf, MAXSTRLEN-1, f) || buf[strlen(buf)-1]=='\n') break;      \
}


/* reads (from "fp") "nvals" doubles into "vals".
 * Returns number of doubles actually read, EOF on end of file, EOF-1 on error
 */
static int readNDoubles(FILE *fp, double *vals, int nvals)
{
register int i;
int n, j;

  for(i=n=0; i<nvals; ++i){
    j=fscanf(fp, "%lf", vals+i);
    if(j==EOF) return EOF;

    if(j!=1 || ferror(fp)) return EOF-1;

    n+=j;
  }

  return n;
}

/* reads (from "fp") "nvals" doubles without storing them.
 * Returns EOF on end of file, EOF-1 on error
 */
static int skipNDoubles(FILE *fp, int nvals)
{
register int i;
int j;

  for(i=0; i<nvals; ++i){
    j=fscanf(fp, "%*f");
    if(j==EOF) return EOF;

    if(ferror(fp)) return EOF-1;
  }

  return nvals;
}

/* reads (from "fp") "nvals" ints into "vals".
 * Returns number of ints actually read, EOF on end of file, EOF-1 on error
 */
static int readNInts(FILE *fp, int *vals, int nvals)
{
register int i;
int n, j;

  for(i=n=0; i<nvals; ++i){
    j=fscanf(fp, "%d", vals+i);
    if(j==EOF) return EOF;
    
    if(j!=1 || ferror(fp)) return EOF-1;

    n+=j;
  }

  return n;
}

/* returns the number of cameras defined in a camera parameters file.
 * Each line of the file corresponds to the parameters of a single camera
 */
static int findNcameras(FILE *fp)
{
int lineno, ncams, ch;

  lineno=ncams=0;
  while(!feof(fp)){
    if((ch=fgetc(fp))=='#'){ /* skip comments */
      SKIP_LINE(fp);
      ++lineno;
      continue;
    }

    if(feof(fp)) break;

    SKIP_LINE(fp);
    ++lineno;
    if(ferror(fp)){
      fprintf(stderr, "findNcameras(): error reading input file, line %d\n", lineno);
      exit(1);
    }
    ++ncams;
  }

  return ncams;
}


/* returns the number of doubles found in the first line of the supplied text file. rewinds file.
 */
static int countNDoubles(FILE *fp)
{
int lineno, ch, np, i;
char buf[MAXSTRLEN], *s;
double dummy;

  lineno=0;
  while(!feof(fp)){
    if((ch=fgetc(fp))=='#'){ /* skip comments */
      SKIP_LINE(fp);
      ++lineno;
      continue;
    }

    if(feof(fp)) return 0;

    ++lineno;
    if(!fgets(buf, MAXSTRLEN-1, fp)){ /* read the line found... */
      fprintf(stderr, "countNDoubles(): error reading input file, line %d\n", lineno);
      exit(1);
    }
    /* ...and count the number of doubles it has */
    for(np=i=0, s=buf; 1 ; ++np, s+=i){
      ch=sscanf(s, "%lf%n", &dummy, &i);
      if(ch==0 || ch==EOF) break;
    }

    rewind(fp);
    return np;
  }

  return 0; // should not reach this point
}

/* reads into "params" the camera parameters defined in a camera parameters file.
 * "params" is assumed to point to sufficiently large memory.
 * Each line contains the parameters of a single camera, "cnp" parameters per camera
 *
 * infilter points to a function that converts a motion parameters vector stored
 * in a file to the format expected by eucsbademo. For instance, it can convert
 * a four element quaternion to the vector part of its unit norm equivalent.
 */
static void readCameraParams(FILE *fp, int cnp,
                             void (*infilter)(double *pin, int nin, double *pout, int nout), int filecnp,
                             double *params)
{
int lineno, n, ch;
double *tofilter;

  if(filecnp>0 && infilter){
    if((tofilter=(double *)malloc(filecnp*sizeof(double)))==NULL){;
      fprintf(stderr, "memory allocation failed in readInitialSBAEstimate()\n");
      exit(1);
    }
  }
  else{ // camera params will be used as read
    infilter=NULL;
    tofilter=NULL;
    filecnp=cnp;
  }
  /* make sure that the parameters file contains the expected number of parameters per line */
  if((n=countNDoubles(fp))!=filecnp){
    fprintf(stderr, "readCameraParams(): expected %d camera parameters, first line contains %d!\n", filecnp, n);
    exit(1);
  }

  lineno=0;
  while(!feof(fp)){
    if((ch=fgetc(fp))=='#'){ /* skip comments */
      SKIP_LINE(fp);
      ++lineno;
      continue;
    }

    if(feof(fp)) break;

    ungetc(ch, fp);
    ++lineno;
    if(infilter){
      n=readNDoubles(fp, tofilter, filecnp);
      (*infilter)(tofilter, filecnp, params, cnp);
    }
    else
      n=readNDoubles(fp, params, cnp);
    if(n==EOF) break;
    if(n!=filecnp){
      fprintf(stderr, "readCameraParams(): line %d contains %d parameters, expected %d!\n", lineno, n, filecnp);
      exit(1);
    }
    if(ferror(fp)){
      fprintf(stderr, "findNcameras(): error reading input file, line %d\n", lineno);
      exit(1);
    }

    params+=cnp;
  }
  free(tofilter);
}


/* determines the number of 3D points contained in a points parameter file as well as the
 * total number of their 2D image projections across all images. Each 3D point is assumed
 * to be described by "pnp" parameters. The file format is
 * X_0...X_{pnp-1}  nframes  frame0 x0 y0  frame1 x1 y1 ...
 */
static void readNpointsAndNprojections(FILE *fp, int *n3Dpts, int pnp, int *nprojs)
{
int lineno, npts, nframes, ch, n;

  *n3Dpts=*nprojs=lineno=npts=0;
  while(!feof(fp)){
    if((ch=fgetc(fp))=='#'){ /* skip comments */
      SKIP_LINE(fp);
      ++lineno;
      continue;
    }

    if(feof(fp)) break;

    ungetc(ch, fp);
    ++lineno;
    skipNDoubles(fp, pnp);
    n=readNInts(fp, &nframes, 1);
    if(n!=1){
      fprintf(stderr, "readNpointsAndNprojections(): error reading input file, line %d: "
                      "expecting number of frames for 3D point\n", lineno);
      exit(1);
    }
    SKIP_LINE(fp);
    *nprojs+=nframes;
    ++npts;
  }

  *n3Dpts=npts;
}


/* reads a points parameter file.
 * "params", "projs" & "vmask" are assumed preallocated, pointing to
 * memory blocks large enough to hold the parameters of 3D points, 
 * their projections in all images and the point visibility mask, respectively.
 * Each 3D point is assumed to be defined by pnp parameters.
 * File format is X_{0}...X_{pnp-1}  nframes  frame0 x0 y0  frame1 x1 y1 ...
 */
static void readPointParamsAndProjections(FILE *fp, double *params, int pnp, double *projs, char *vmask, int ncams)
{
int nframes, ch, lineno, ptno, frameno, n;
const int mnp=2;
register int i;

  lineno=ptno=0;
  while(!feof(fp)){
    if((ch=fgetc(fp))=='#'){ /* skip comments */
      SKIP_LINE(fp);
      lineno++;

      continue;
    }

    if(feof(fp)) break;

    ungetc(ch, fp);

    n=readNDoubles(fp, params, pnp); /* read in point parameters */
    if(n==EOF) break;
    if(n!=pnp){
      fprintf(stderr, "readPointParamsAndProjections(): error reading input file, line %d: expecting %d parameters for 3D point\n",
                       lineno, pnp);
      exit(1);
    }
    params+=pnp;

    n=readNInts(fp, &nframes, 1);  /* read in number of image projections */
    if(n!=1){
      fprintf(stderr, "readPointParamsAndProjections(): error reading input file, line %d: expecting number of frames for 3D point\n",
                       lineno);
      exit(1);
    }

    for(i=0; i<nframes; ++i){
      n=readNInts(fp, &frameno, 1); /* read in frame number... */
      n+=readNDoubles(fp, projs, mnp); /* ...and image projection */
      if(n!=mnp+1){
        fprintf(stderr, "readPointParamsAndProjections(): error reading image projections from line %d [n=%d].\n"
                        "Line contains fewer than %d projections?\n", lineno+1, n, nframes);
        exit(1);
      }

      if(frameno>=ncams){
        fprintf(stderr, "readPointParamsAndProjections(): line %d contains an image projection for frame %d "
                        "but only %d cameras have been specified!\n", lineno+1, frameno, ncams);
        exit(1);
      }

      projs+=mnp;
      vmask[ptno*ncams+frameno]=1;
    }

    fscanf(fp, "\n"); // consume trailing newline

    lineno++;
    ptno++;
  }
}


/* combines the above routines to read the initial estimates of the motion + structure parameters from text files.
 * Also, it loads the projections of 3D points across images. The routine dynamically allocates the required amount
 * of memory (last 3 arguments).
 */
void readInitialSBAEstimate(char *camsfname, char *ptsfname, int cnp, int pnp, int mnp,
                            void (*infilter)(double *pin, int nin, double *pout, int nout), int filecnp,
                            int *ncams, int *n3Dpts, int *n2Dprojs,
                            double **motstruct, double **imgpts, char **vmask)
{
FILE *fpc, *fpp;

  if((fpc=fopen(camsfname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s, exiting\n", camsfname);
    exit(1);
  }

  if((fpp=fopen(ptsfname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s, exiting\n", ptsfname);
    exit(1);
  }

  *ncams=findNcameras(fpc);
  readNpointsAndNprojections(fpp, n3Dpts, pnp, n2Dprojs);

  *motstruct=(double *)malloc((*ncams*cnp + *n3Dpts*pnp)*sizeof(double));
  if(*motstruct==NULL){
    fprintf(stderr, "memory allocation failed in readInitialSBAEstimate()\n");
    exit(1);
  }
  *imgpts=(double *)malloc(*n2Dprojs*mnp*sizeof(double));
  if(*imgpts==NULL){
    fprintf(stderr, "memory allocation failed in readInitialSBAEstimate()\n");
    exit(1);
  }
  *vmask=(char *)malloc(*n3Dpts * *ncams * sizeof(char));
  if(*vmask==NULL){
    fprintf(stderr, "memory allocation failed in readInitialSBAEstimate()\n");
    exit(1);
  }
  memset(*vmask, 0, *n3Dpts * *ncams * sizeof(char)); /* clear vmask */


  /* prepare for re-reading files */
  rewind(fpc);
  rewind(fpp);

  readCameraParams(fpc, cnp, infilter, filecnp, *motstruct);
  readPointParamsAndProjections(fpp, *motstruct+*ncams*cnp, pnp, *imgpts, *vmask, *ncams);

  fclose(fpc);
  fclose(fpp);
}

/* reads the 3x3 intrinsic calibration matrix contained in a file */
void readCalibParams(char *fname, double ical[9])
{
  FILE *fp;
  int i, ch=EOF;

  if((fp=fopen(fname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s, exiting\n", fname);
    exit(1);
  }

  while(!feof(fp) && (ch=fgetc(fp))=='#') /* skip comments */
    SKIP_LINE(fp);

  if(feof(fp)){
    ical[0]=ical[1]=ical[2]=ical[3]=ical[4]=ical[5]=ical[6]=ical[7]=ical[8]=0.0;
    return;
  }

  ungetc(ch, fp);

  for(i=0; i<3; i++){
    fscanf(fp, "%lf%lf%lf\n", ical, ical+1, ical+2);
    ical+=3;
  }

  fclose(fp);
}

/* routines for printing the motion and structure parameters, plus the projections
 * of 3D points across images. Mainly for debugging purposes.
 *
 * outfilter points to a function that converts a motion parameters vector from 
 * the interrnal representation used by eucsbademo to a format suitable for display.
 * For instance, it can expand a quaternion to 4 elements from its vector part.
 */

/* motion parameters only */
void printSBAMotionData(double *motstruct, int ncams, int cnp,
                        void (*outfilter)(double *pin, int nin, double *pout, int nout), int outcnp)
{
register int i;

  printf("Motion parameters:\n");
  if(!outfilter || outcnp<=0){ // no output filtering
    for(i=0; i<ncams*cnp; ++i){
      printf("%lf ", motstruct[i]);
      if((i+1)%cnp==0) putchar('\n');
    }
  }
  else{
    register int j;
    double *filtered;

    if((filtered=(double *)malloc(outcnp*sizeof(double)))==NULL){;
      fprintf(stderr, "memory allocation failed in printSBAMotionData()\n");
      exit(1);
    }
    for(i=0; i<ncams*cnp; i+=cnp){
      (*outfilter)(motstruct+i, cnp, filtered, outcnp);
      for(j=0; j<outcnp; ++j)
        printf("%lf ", filtered[j]);
      putchar('\n');
    }
    free(filtered);
  }
}

/* structure parameters only */
void printSBAStructureData(double *motstruct, int ncams, int n3Dpts, int cnp, int pnp)
{
register int i;

  motstruct+=ncams*cnp;
  printf("\n\nStructure parameters:\n");
  for(i=0; i<n3Dpts*pnp; ++i){
    printf("%lf ", motstruct[i]);
    if((i+1)%pnp==0) putchar('\n');
  }
}

/* prints the estimates of the motion + structure parameters. It also prints the projections
 * of 3D points across images.
 */
void printSBAData(double *motstruct, int cnp, int pnp, int mnp, 
                  void (*outfilter)(double *pin, int nin, double *pout, int nout), int outcnp,
                  int ncams, int n3Dpts, double *imgpts, int n2Dprojs, char *vmask)
{
register int i, j, k, l;
int nframes;

  printSBAMotionData(motstruct, ncams, cnp, outfilter, outcnp);
  motstruct+=ncams*cnp;

  printf("\n\nStructure parameters and image projections:\n");
  for(i=k=0; i<n3Dpts; ++i){
    for(j=0; j<pnp; ++j)
      printf("%lf ", motstruct[i*pnp+j]);
    
    for(j=nframes=0; j<ncams; ++j)
      if(vmask[i*ncams+j]) ++nframes;
    printf("%d", nframes);
    for(j=0; j<ncams; ++j)
      if(vmask[i*ncams+j]){
        printf(" %d", j);
        for(l=0; l<mnp; ++l, ++k)
          printf(" %lf", imgpts[k]);
      }
    putchar('\n');
  }

#if 0
  /* alternative output format */
  printf("\n\nStructure parameters:\n");
  for(i=0; i<n3Dpts*pnp; ++i){
    printf("%lf ", motstruct[i]);
    if((i+1)%pnp==0) putchar('\n');
  }

  printf("\n\nImage projections:\n");
  for(i=0; i<n2Dprojs*mnp; ++i){
    printf("%lf ", imgpts[i]);
  }

  printf("\n\nVisibility mask\n");
  for(i=0; i<ncams*n3Dpts; ++i)
    printf("%d%s", (int)vmask[i], ((i+1)%ncams)? " " : "\n");
  putchar('\n');
#endif

  fflush(stdout);
}
