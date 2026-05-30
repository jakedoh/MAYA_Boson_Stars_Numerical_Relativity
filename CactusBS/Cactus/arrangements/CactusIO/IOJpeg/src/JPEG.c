#include "cctk.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 

#include "ioJpegGH.h"

#include "jconfig.h"
#include "jpeglib.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusIO_IOJpeg_JPEG_c)

typedef struct jpeg_compress_struct JpgComp;

#ifndef JpgErr
typedef struct jpeg_error_mgr JpgErr;
#endif

/* prototypes of routines defined in this source file */
GLOBAL(void)
jpeg_memory_dest (j_compress_ptr cinfo, JOCTET *buffer,int bufsize);


/*
  Image data is an array of unsigned character array of
  RGB data.  The data is stored in interleaved order.
  IE, the first three elements are a byte of Red followed
  by a byte of Green and then a byte of Blue for the first
  pixel.  Data is stored in fortran order (ie. x is fastest
  moving dimension).
 */
int WriteJPEGToFileRGB(int nx, /* width of image in pixels */
                       int ny, /* height of the image in pixels */
                       void *data, /* buffer containing image data */
                       int Quality, /* Integer from 0 to 100 */
                       FILE* outfile){  /* name of file to store in */
  JpgComp cinfo;
  JpgErr jerr;  
  /*  FILE * outfile;*/
  unsigned char *dataRGB = (unsigned char *)data;
  JSAMPROW row_pointer=(JSAMPROW)dataRGB;
  
  memset (&cinfo,0,sizeof(cinfo));
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
        
  /* Setup JPEG */
  cinfo.image_width = nx ;      /* image width and height, in pixels */
  cinfo.image_height = ny;
  cinfo.input_components = 3;   /* # of color components per pixel=3 RGB */
  cinfo.in_color_space = JCS_RGB;
  /*  if ((outfile = fopen(FileName, "wb")) == NULL) {
    printf("Cannot open file [%s]\n",FileName);
    return 0; 
  } */
  jpeg_stdio_dest(&cinfo, outfile);
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality (&cinfo,Quality,TRUE);
  /* Starting compress */
  jpeg_start_compress(&cinfo, TRUE);
  /* Now compress everything one scanline at-a-time */
  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer = (JSAMPROW)(dataRGB+(cinfo.next_scanline*3*nx)); /* in bytes or words? */
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);
  /* All done! */
  /*  fclose(outfile);*/
  return 1;
}

/*--------------
  A hack to hijack JPEG's innards to write into a memory buffer
----------------
/  this defines a new destination manager to store images in memory
/  derived by jdatadst.c */
typedef struct {
  struct jpeg_destination_mgr pub;      /* public fields */
  JOCTET *buffer;                                       /* start of buffer */
  int bufsize;                                          /* buffer size */
  int datacount;                                        /* finale data size */
} memory_destination_mgr;

typedef memory_destination_mgr *mem_dest_ptr;

/*----------------------------------------------------------------------------
  /  Initialize destination --- called by jpeg_start_compress before any data is actually written. */

METHODDEF(void)
init_destination (j_compress_ptr cinfo)
{
  mem_dest_ptr dest = (mem_dest_ptr) cinfo->dest;

  dest->pub.next_output_byte = dest->buffer;
  dest->pub.free_in_buffer = dest->bufsize;
  dest->datacount=0;
}



/*----------------------------------------------------------------------------
  /  Empty the output buffer --- called whenever buffer fills up. */
METHODDEF(boolean)
empty_output_buffer (j_compress_ptr cinfo)
{
  mem_dest_ptr dest = (mem_dest_ptr) cinfo->dest;

  dest->pub.next_output_byte = dest->buffer;
  dest->pub.free_in_buffer = dest->bufsize;

  return TRUE;
}


/*----------------------------------------------------------------------------

  /  Terminate destination --- called by jpeg_finish_compress
  /  after all data has been written.  Usually needs to flush buffer. */
METHODDEF(void)
term_destination (j_compress_ptr cinfo)
{
  /* expose the finale compressed image size */
  
  mem_dest_ptr dest = (mem_dest_ptr) cinfo->dest;
  dest->datacount = dest->bufsize - dest->pub.free_in_buffer;
  
}

/*----------------------------------------------------------------------------
/ Prepare for output to a memory buffer. The caller must have allocate memory
/ to store the compressed image, and supply its size */
GLOBAL(void)
jpeg_memory_dest (j_compress_ptr cinfo, JOCTET *buffer,int bufsize)
{
  mem_dest_ptr dest;
  if (cinfo->dest == NULL) {    /* first time for this JPEG object? */
    cinfo->dest = (struct jpeg_destination_mgr *)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,
                                  sizeof(memory_destination_mgr));
  }

  dest = (mem_dest_ptr) cinfo->dest;
  dest->bufsize=bufsize;
  dest->buffer=buffer;
  dest->pub.init_destination = init_destination;
  dest->pub.empty_output_buffer = empty_output_buffer;
  dest->pub.term_destination = term_destination;
}

/********************************************
Identical in nearly every way to WriteJPEGToFileRGB(), but
it writes into the memory buffer specified.  To be safe, its
good to make the memorybuffer the same size as the input image
+ 1024 bytes.  It is guaranteed that the image will be less
than this size.  In fact, if you use a typical "quality" level
of 75, you can get away with an image which is one quarter that
size.

 ******************************************** */
int WriteJPEGToMemoryRGB(int nx,int ny, void *data, int Quality, char *memorybuffer,int bufsize){       
  JpgComp cinfo;
  JpgErr jerr;  
  unsigned char *dataRGB = (unsigned char *)data;
  JSAMPROW row_pointer=(JSAMPROW)dataRGB;
  JOCTET *jpgbuff;
  mem_dest_ptr dest;
  int csize=0;

  /* zero out the compresion info structures and
     allocate a new compressor handle */
  memset (&cinfo,0,sizeof(cinfo));
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
 
  /* Setup JPEG datastructures */
  cinfo.image_width = nx ;      /* image width and height, in pixels */
  cinfo.image_height = ny;
  cinfo.input_components = 3;   /* # of color components per pixel=3 RGB */
  cinfo.in_color_space = JCS_RGB;               
  jpgbuff = (JOCTET*)memorybuffer;

  /* Setup compression and do it */
  jpeg_memory_dest(&cinfo,jpgbuff,bufsize);
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality (&cinfo,Quality,TRUE);
  jpeg_start_compress(&cinfo, TRUE);
  /* compress each scanline one-at-a-time */
  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer = (JSAMPROW)(dataRGB+(cinfo.next_scanline*3*nx));
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }
  jpeg_finish_compress(&cinfo);
  /* Now extract the size of the compressed buffer */
  dest=(mem_dest_ptr)cinfo.dest;
  csize=dest->datacount; /* the actual compressed datasize */
  /* destroy the compressor handle */
  jpeg_destroy_compress(&cinfo);
  return csize;
}

#define ClampToByte(f) (unsigned char)(f=(f>255)?255:f)
/* Werner's nifty non-linear red-to-blue auto-colormap.
   I faked a CCTK_REAL, you can make it the right thing
   when integrated with cactus.

   Anyways, the "rdfac" concentrates the colormap when you
   use higher values.  Typical default value is 32 based
   on the current JPEG thorn.

   Next step is to embed some IDL colormaps into the
   program.  I can do that a bit later.
*/

void AutoColorDataSlice(int nx,int ny, /* size of the image x & y */
                        const CCTK_REAL *datain, /* 2D slice of data input */
                        unsigned char *dataout, /* RGB image data output */
                        CCTK_REAL min,CCTK_REAL max, /* range of the entire 3D dataset
                                                  This could be ranged based 
                                                  on the values of the slice,
                                                  but then that would be a
                                                  bit untrustworthy.  Its
                                                  best to pass in the 
                                                  range for the entire 
                                                  dataset or a pre-defined
                                                  fixed range.  It does
                                                  handle clamping of the
                                                  range. */
                        CCTK_REAL bias,
                        int rdfac){
  /* Bias allows you to move the spectrum from the Red to the
     Blue.  A value of 0.5 is the median.  So 0.4 biases towards
     the Red and 0.6 would bias a little towards the blue.
     The range for this parameter would be -1.0 to +1.0. */
  /*  CCTK_REAL bias=0.4;*/
  int i,last;
  CCTK_REAL F=(CCTK_REAL)rdfac; /* cast to CCTK_REAL... don't know how the original worked at all without a cast */
  for(i=0,last=nx*ny;i<last;i++,dataout+=3){
    CCTK_REAL f;

    /* check for division by 0 */
    if (min == max)
    {
      f = bias;
    }
    else
    {
      f = bias-(*datain++ - min)/(max-min);
    }
    /* f-=(max-min); zero-center it */
    /* well it can't be less than 0 */
    if(f>0){
      f*=F;
      /* Color components biased towards blue */
      dataout[0]=ClampToByte(f);
      f*=F;
      dataout[1]=ClampToByte(f);
      f*=F;
      dataout[2]=ClampToByte(f);
    }
    else { /* f<0 */
      f=-f;
      f*=F;
      /* reverse color components to bias towards red */
      dataout[2]=ClampToByte(f);
      f*=F;
      dataout[1]=ClampToByte(f);
      f*=F;
      dataout[0]=ClampToByte(f);
    }
  }
}

#if 0
int main(int argc,char *argv[]){
  /* OK, lets create a bogus image */
  int nx=512,ny=512;
  unsigned char *datargb=(unsigned char *)malloc(3*nx*ny);
  CCTK_REAL *data= (CCTK_REAL*)malloc(nx*ny*sizeof(CCTK_REAL));
  int i,j,idx=0;
  CCTK_REAL radius=((CCTK_REAL)(nx/4));
  CCTK_REAL min=(CCTK_REAL)(nx*ny),max=-1.0;
  int bufsize = nx*ny*3+1024;
  char *memorybuffer=(char *)malloc(bufsize); /* safe size for membuf */
  int compressed_size=0;
  FILE *outfile;
  /* compute a circle by the most inefficient means possible */
  for(j=0,idx=0;j<ny;j++){
    CCTK_REAL y=((CCTK_REAL)(j-ny/2));
    y*=y;
    /* fprintf(stderr,"\n"); */
    for(i=0;i<nx;i++,idx++){
      CCTK_REAL x = ((CCTK_REAL)(i-nx/2));
      CCTK_REAL val;
      x*=x;
      data[idx]= val=sqrt(x+y);
     
      if(val>max) max=val;
      if(val<min) min=val;
    }
  }
  printf("Slice Data Min=%g Max=%g\n",min,max);
  /* Autocolor extracted from Werner's JPEG thorn */
  AutoColorDataSlice(nx,ny,data,datargb,min,max,32);
  WriteJPEGToFileRGB(nx,ny, datargb,75, "write2file.jpg");
  /* write to mem: It returns the size of the compressed image 
     cdntained in the memorybuffer */
  compressed_size = WriteJPEGToMemoryRGB(nx,ny, datargb,75, memorybuffer,bufsize);
  if ((outfile = fopen("write2mem.jpg", "wb")) == NULL) {
    printf("Cannot open file write2mem.jpg\n");
    return 0; /* failure */
  } 
  fwrite(memorybuffer,1,compressed_size,outfile);
  fclose(outfile);
  return 1;
}
#endif
