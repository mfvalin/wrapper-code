
#include <stdio.h>
#include <stdlib.h>
#include "jasper/jasper.h"

int TO_jpeg2000(unsigned char *cin,int width,int height,int nbits,
                 int ltype, int ratio, int retry, char *outjpc, 
                 int jpclen)
/*$$$  SUBPROGRAM DOCUMENTATION BLOCK
*                .      .    .                                       .
* SUBPROGRAM:    TO_jpeg2000      Encodes JPEG2000 code stream
*   PRGMMR: Gilbert          ORG: W/NP11     DATE: 2002-12-02
*
* ABSTRACT: This Function encodes a grayscale image into a JPEG2000 code stream
*   specified in the JPEG2000 Part-1 standard (i.e., ISO/IEC 15444-1) 
*   using JasPer Software version 1.700.2  or better written by the 
*   University of British Columbia, Image Power Inc, and others.
*   JasPer is available at http://www.ece.uvic.ca/~mdadams/jasper/.
*
* PROGRAM HISTORY LOG:
* 2002-12-02  Gilbert
* 2004-12-16  Gilbert - Added retry argument/option to allow option of
*                       increasing the maximum number of guard bits to the
*                       JPEG2000 algorithm.
* 2015-10-26  M.Valin - removed GRIB library connections, light refactoring,
*                       changed function name
*
* USAGE:    int TO_jpeg2000(unsigned char *cin,int width,int height,
*                            int nbits, int ltype, int ratio, 
*                            int retry, char *outjpc, int jpclen)
*
*   INPUT ARGUMENTS:
*      cin   - Packed matrix of Grayscale image values to encode.
*              (byte stream, MSB first, LSB last for each multibyte item)
*     width  - width of image
*     height - height of image
*     nbits  - depth (in bits) of image.  i.e number of bits
*              used to hold each data value
*    ltype   - indicator of lossless or lossy compression
*              = 1, for lossy compression
*              != 1, for lossless compression
*    ratio   - target compression ratio.  (ratio:1)
*              Used only when ltype == 1.
*    retry   - Pointer to option type.
*              1 = try increasing number of guard bits
*              otherwise, no additional options
*    jpclen  - Number of bytes allocated for new JPEG2000 code stream in
*              outjpc.
*
*   INPUT ARGUMENTS:
*     outjpc - Output encoded JPEG2000 code stream
*
*   RETURN VALUES :
*        > 0 = Length in bytes of encoded JPEG2000 code stream
*         -3 = Error decode jpeg2000 code stream.
*         -5 = decoded image had multiple color components.
*              Only grayscale is allowed.
*
* REMARKS:
*
*      Requires JasPer Software version 1.700.2 or better
*
* ATTRIBUTES:
*   LANGUAGE: C
*   OS:  Linux
*
*$$$*/
{
    int ier,rwcnt;
    jas_image_t image;
    jas_stream_t *jpcstream,*istream;
    jas_image_cmpt_t cmpt,*pcmpt;
#define MAXOPTSSIZE 1024
    char opts[MAXOPTSSIZE];

/*
    printf(" TO_jpeg2000:width %ld\n",width);
    printf(" TO_jpeg2000:height %ld\n",height);
    printf(" TO_jpeg2000:nbits %ld\n",nbits);
    printf(" TO_jpeg2000:jpclen %ld\n",jpclen);
*/
//    jas_init();

//
//    Set lossy compression options, if requested.
//
    if ( ltype != 1 ) {
       opts[0]=(char)0;
    }
    else {
       snprintf(opts,MAXOPTSSIZE,"mode=real\nrate=%f",1.0/(float)ratio);
    }
    if ( retry == 1 ) {             // option to increase number of guard bits
       strcat(opts,"\nnumgbits=4");
    }
    //printf("SAGopts: %s\n",opts);
    
//
//     Initialize the JasPer image structure describing the grayscale
//     image to encode into the JPEG2000 code stream.
//
    image.tlx_=0;
    image.tly_=0;
    image.brx_=(jas_image_coord_t)width;
    image.bry_=(jas_image_coord_t)height;
    image.numcmpts_=1;
    image.maxcmpts_=1;
    image.clrspc_=JAS_CLRSPC_SGRAY;         /* grayscale Image */
    image.cmprof_=0; 
    image.inmem_=1;

    cmpt.tlx_=0;
    cmpt.tly_=0;
    cmpt.hstep_=1;
    cmpt.vstep_=1;
    cmpt.width_=(jas_image_coord_t)width;
    cmpt.height_=(jas_image_coord_t)height;
    cmpt.type_=JAS_IMAGE_CT_COLOR(JAS_CLRSPC_CHANIND_GRAY_Y);
    cmpt.prec_=nbits;
    cmpt.sgnd_=0;
    cmpt.cps_=(nbits+7)/8;

    pcmpt=&cmpt;
    image.cmpts_=&pcmpt;

#if defined(DEBUG)
    printf(" SAGOUT ENCODE:\n");
    printf(" tlx %d \n",image.tlx_);
    printf(" tly %d \n",image.tly_);
    printf(" brx %d \n",image.brx_);
    printf(" bry %d \n",image.bry_);
    printf(" numcmpts %d \n",image.numcmpts_);
    printf(" maxcmpts %d \n",image.maxcmpts_);
    printf(" colorspace %d \n",image.clrspc_);
    printf(" inmem %d \n",image.inmem_);
    printf(" COMPONENT:\n");
    printf(" tlx %d \n",pcmpt->tlx_);
    printf(" tly %d \n",pcmpt->tly_);
    printf(" hstep %d \n",pcmpt->hstep_);
    printf(" vstep %d \n",pcmpt->vstep_);
    printf(" width %d \n",pcmpt->width_);
    printf(" height %d \n",pcmpt->height_);
    printf(" prec %d \n",pcmpt->prec_);
    printf(" sgnd %d \n",pcmpt->sgnd_);
    printf(" cps %d \n",pcmpt->cps_);
    printf(" type %d \n",pcmpt->type_);
#endif
//
//    Open a JasPer stream containing the input grayscale values
//
    istream=jas_stream_memopen((char *)cin,height*width*cmpt.cps_);
    cmpt.stream_=istream;

//
//    Open an output stream that will contain the encoded jpeg2000
//    code stream.
//
    jpcstream=jas_stream_memopen(outjpc,(int)jpclen);

//
//     Encode image.
//
    ier=jpc_encode(&image,jpcstream,opts);
    if ( ier != 0 ) {
       printf(" jpc_encode return = %d \n",ier);
       return -3;
    }
//
//     Clean up JasPer work structures.
//    
    rwcnt=jpcstream->rwcnt_;
    ier=jas_stream_close(istream);
    ier=jas_stream_close(jpcstream);
//
//      Return size of jpeg2000 code stream
//
    return (rwcnt);

}

   int FROM_jpeg2000(char *injpc,int bufsize,unsigned int *outfld)
/*$$$  SUBPROGRAM DOCUMENTATION BLOCK
*                .      .    .                                       .
* SUBPROGRAM:    FROM_jpeg2000      Decodes JPEG2000 code stream
*   PRGMMR: Gilbert          ORG: W/NP11     DATE: 2002-12-02
*
* ABSTRACT: This Function decodes a JPEG2000 code stream specified in the
*   JPEG2000 Part-1 standard (i.e., ISO/IEC 15444-1) using JasPer 
*   Software version 1.700.2 or better written by the University of British
*   Columbia and Image Power Inc, and others.
*   JasPer is available at http://www.ece.uvic.ca/~mdadams/jasper/.
*
* PROGRAM HISTORY LOG:
* 2002-12-02  Gilbert
* 2015-10-26  M.Valin - removed GRIB library connections, light refactoring,
*                       changed function name
*
* USAGE:     int FROM_jpeg2000(char *injpc,int bufsize,int *outfld)
*
*   INPUT ARGUMENTS:
*      injpc - Input JPEG2000 code stream.
*    bufsize - Length (in bytes) of the input JPEG2000 code stream.
*
*   OUTPUT ARGUMENTS:
*     outfld - Output matrix of grayscale image values.
*
*   RETURN VALUES :
*          0 = Successful decode
*         -3 = Error decode jpeg2000 code stream.
*         -5 = decoded image had multiple color components.
*              Only grayscale is allowed.
*
* REMARKS:
*
*      Requires JasPer Software version 1.700.2 or better
*
* ATTRIBUTES:
*   LANGUAGE: C
*   OS:  Linux
*
*$$$*/

{
    int ier;
    int i,j,k;
    jas_image_t *image=0;
    jas_stream_t *jpcstream;
    jas_image_cmpt_t *pcmpt;
    char *opts=0;
    jas_matrix_t *data;
#define ROW unsigned long
//    ROW *my_row;
    unsigned int value;

//    jas_init();

    ier=0;
//   
//     Create jas_stream_t containing input JPEG200 codestream in memory.
//       

    jpcstream=jas_stream_memopen(injpc,bufsize);

//   
//     Decode JPEG200 codestream into jas_image_t structure.
//       
    image=jpc_decode(jpcstream,opts);
    if ( image == 0 ) {
       printf(" jpc_decode return\n");
       return -3;
    }
    
    pcmpt=image->cmpts_[0];
#if defined(DEBUG)
    printf(" SAGOUT DECODE:\n");
    printf(" tlx %d \n",image->tlx_);
    printf(" tly %d \n",image->tly_);
    printf(" brx %d \n",image->brx_);
    printf(" bry %d \n",image->bry_);
    printf(" numcmpts %d \n",image->numcmpts_);
    printf(" maxcmpts %d \n",image->maxcmpts_);
    printf(" colorspace %d \n",image->clrspc_);
    printf(" inmem %d \n",image->inmem_);
    printf(" COMPONENT:\n");
    printf(" tlx %d \n",pcmpt->tlx_);
    printf(" tly %d \n",pcmpt->tly_);
    printf(" hstep %d \n",pcmpt->hstep_);
    printf(" vstep %d \n",pcmpt->vstep_);
    printf(" width %d \n",pcmpt->width_);
    printf(" height %d \n",pcmpt->height_);
    printf(" prec %d \n",pcmpt->prec_);
    printf(" sgnd %d \n",pcmpt->sgnd_);
    printf(" cps %d \n",pcmpt->cps_);
    printf(" type %d \n",pcmpt->type_);
#endif
//   Expecting jpeg2000 image to be grayscale only.
//   No color components.
//
    if (image->numcmpts_ != 1 ) {
       printf("FROM_jpeg2000: Found color image.  Grayscale expected.\n");
       return (-5);
    }

// 
//    Create a data matrix of grayscale image values decoded from
//    the jpeg2000 codestream.
//
    data=jas_matrix_create(jas_image_height(image), jas_image_width(image));
    jas_image_readcmpt(image,0,0,0,jas_image_width(image),
                       jas_image_height(image),data);
//
//    Copy data matrix to output integer array.
//
    k=0;
//    printf ("deltap = %d\n",data->rows_[1]-data->rows_[0]);
    for (i=0;i<pcmpt->height_;i++) {
//      my_row = (ROW *)data->rows_[i];
//      printf ("p = %16p\n",my_row);
      for (j=0;j<pcmpt->width_;j++) {
//        value = my_row[j];
//        outfld[k++] = my_row[j];
//        printf ("%4ld ",my_row[j]);
//        printf ("%16.16x ",my_row[j]);
        outfld[k++]=data->rows_[i][j];
      }
//      printf ("%d \n",k);
    }
//
//     Clean up JasPer work structures.
//
    jas_matrix_destroy(data);
    ier=jas_stream_close(jpcstream);
    jas_image_destroy(image);

    return 0;

}
#if defined(SELF_TEST)
// int TO_jpeg2000(unsigned char *cin,int width,int height,int nbits,
//                  int ltype, int ratio, int retry, char *outjpc, 
//                  int jpclen)
//    int FROM_jpeg2000(char *injpc,int bufsize,int *outfld)
main()
{
  unsigned int matrix[256];
  unsigned int restored[256];
  int packed[1024];
  int nbytes, status;
  int i, j;
  for(i=0 ; i<256 ; i++) { matrix[i] = i<<24 ; }
  for(j=0 ; j<16 ; j++){
    for(i=0 ; i<16 ; i++) {
      fprintf(stderr,"%4d",matrix[i+16*j]>>24);
    }
    fprintf(stderr,"\n");
  }
  nbytes = TO_jpeg2000((unsigned char *) matrix, 16,16,25,0,1,1,(char *)packed,sizeof(packed));
  fprintf(stderr,"nbytes=%d\n",nbytes);

  status = FROM_jpeg2000((char *)packed,nbytes,restored);
  fprintf(stderr,"status=%d\n",status);

  for(j=0 ; j<16 ; j++){
    for(i=0 ; i<16 ; i++) {
      fprintf(stderr,"%4d",restored[i+16*j]);
    }
    fprintf(stderr,"\n");
  }
}
#endif