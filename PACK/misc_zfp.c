//  useful routines for C and FORTRAN programming
//  Copyright (C) 2022  Environnement Canada
//
//  This is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation,
//  version 2.1 of the License.
//
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
// functions using the zfp (de)compressor
// this code relies on zfp version 1.0.0
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <zfp.h>
#include <misc_analyze.h>

#if defined(DEBUG)
static int debug = DEBUG;
#else
static int debug = 0;      // debug messages are OFF by default
#endif
static int diag = 1 ;      // error messages are ON by default

int ZfpCompress_set_debug(int flag){ 
  int oldval = debug ;
  debug = !!flag ; 
  return oldval;
}

int ZfpCompress_set_diag(int flag){ 
  int oldval = diag ;
  diag = !!flag ; 
  return oldval;
}

int32_t zfp_codec_consistent(){
  return zfp_codec_version == ZFP_CODEC ;
}

uint32_t get_zfp_codec_version(){ return zfp_codec_version ; }

uint32_t get_ZFP_CODEC(){ return ZFP_CODEC ; }

// use zfp to compress a 1/2/3 Dimensional float array (Fortran ordering)
// array     [IN] : float array to compress
// nx        [IN] : first storage dimension
// ny        [IN] : second array dimension
// nz        [IN] : third array dimension
// maxerror  [IN] : largest acceptable absolute error
// precision [IN] : number of bit planes to keep (related to relative error)
// Ssize    [OUT] : size in bytes of compressed stream
// return         : pointer to compressed stream
//
// N.B.  it is up to the caller to free the space allocated for the compressed stream
//       see zfp documentation for more precise explanations about precision (Fixed-Precision mode)
//       https://zfp.readthedocs.io/en/release1.0.0/modes.html#fixed-precision-mode
void *ZfpCompress(void* array, int nx, int ny, int nz, float maxerror, int precision , int *Ssize) {
  zfp_type type ;                       // array scalar type
  zfp_field* field = NULL ;             // array meta data
  zfp_stream* zfp = NULL ;              // compressed stream
  void* buffer = NULL ;                 // storage for compressed stream
  size_t bufsize ;                      // byte size of compressed buffer
  bitstream* stream = NULL ;            // bit stream to write to or read from
  size_t zfpsize ;                      // byte size of compressed stream
  double ratio ;
  double tolerance = maxerror ;
  int nd = 0 ;
  size_t nheader ;

  if(nx < 1 || ny < 1 || nz < 1)     return NULL ;  // invalid dimensions

  // allocate meta data for a compressed stream
  zfp = zfp_stream_open(NULL) ;
  if(zfp == NULL) goto error ;

  //  allocate meta data for a float array
  type = zfp_type_float ;                            // set field type and dimensions
  if(nz > 1)             { field = zfp_field_3d(array, type, nx, ny, nz); nd = 3 ; } // 3D array
  if(nz == 1 && ny > 1)  { field = zfp_field_2d(array, type, nx, ny)    ; nd = 2 ; } // 2D array
  if(ny == 1 && nz == 1) { field = zfp_field_1d(array, type, nx)        ; nd = 1 ; } // 1D array
  if(field == NULL) goto error ;
  if(debug) fprintf(stderr, "dimensions = %d  ", nd) ;

  if(precision > 0) {              // precision mode
    ratio = zfp_stream_set_precision(zfp, precision) ;    // bit precision mode (uncompressed bits per value)
    if(debug) fprintf(stderr, "precision in/out = %d, %g ", precision, ratio) ;
  }
  else if(tolerance > 0) {         // absolute error bound mode
    ratio = zfp_stream_set_accuracy(zfp, tolerance) ;
    if(debug) fprintf(stderr, "tolerance in/out = %g, %g  ",tolerance,ratio) ;
  }
  else {
    if(diag) fprintf(stderr, "ERROR: one of maxerror or precision MUST be specified\n");
    goto error ;
  }
  bufsize = zfp_stream_maximum_size(zfp, field) ; // get worst case size of compression buffer
  buffer = malloc(bufsize) ;                      // allocate buffer for compressed stream 
  if(buffer == NULL) goto error ;                 // buffer allocation failed

  stream = stream_open(buffer, bufsize) ;         // open stream
  zfp_stream_set_bit_stream(zfp, stream) ;        // associate stream with buffer
  zfp_stream_rewind(zfp) ;                        // rewind to beginning of stream

  if(  (nheader = zfp_write_header(zfp, field, ZFP_HEADER_FULL)) != 0) {   // store compression metadata into stream
    if(debug) fprintf(stderr, "write header size = %ld  ", nheader);
    zfpsize = zfp_compress(zfp, field) ;
    if (zfpsize == 0) {
      if(diag) fprintf(stderr, "ERROR: compression failed\n") ;
      goto error ;
    }
    *Ssize = zfpsize ;
    ratio = nx*ny*nz*4 ; ratio /= zfpsize ;
    if(debug) fprintf(stderr, "compression ratio = %4.2f ", ratio) ;
  }else{
    if(diag) fprintf(stderr, "ERROR: write_header failed\n") ;
    goto error ;
  }

end:
  /* clean up */
  if(field)  zfp_field_free(field);
  if(zfp)    zfp_stream_close(zfp);
  if(stream) stream_close(stream);

  if(debug) fprintf(stderr, "\n") ;
  return buffer ;

error:
  if(buffer) free(buffer) ;
  buffer = NULL ;
  goto end ;
}

// get the dimensions and strides of a compressed zfp float array
// Stream    [IN] : compressed stream
// d        [OUT] : dimensions (up to 4)
// s        [OUT] : strides (up to 4)
// Ssize     [IN] : size in bytes of the compressed stream
// return   : number of dimensions
// N.B.  elements in d and s arrays corresponding to non existent dimensions are set to 0
int32_t ZfpArrayDims(void *Stream, int *d, int *s, int Ssize){
  size_t zfpsize = Ssize ;
  bitstream* stream = stream_open(Stream, zfpsize) ;  // create bit stream structure
  zfp_stream* zfp   = zfp_stream_open(NULL) ;         // create zfp structure
  zfp_field* field  = zfp_field_alloc() ;             // create field structure
  size_t dims[4] ;
  size_t dimsize ;
  ptrdiff_t strides[4] ;
  zfp_bool is_strided ;
  uint32_t ndims = 0 ;
  size_t nheader ;

  if(stream == NULL) goto error ;
  if(zfp == NULL)    goto error ;
  if(field == NULL)  goto error ;

  zfp_stream_set_bit_stream(zfp, stream) ;            // associate stream to zfp stream
  zfp_stream_rewind(zfp) ;                            // set stream pointer to beginning of stream
  if( (nheader = zfp_read_header(zfp, field, ZFP_HEADER_FULL)) != 0){
    if(debug)fprintf(stderr, "read header size = %ld  ", nheader);
    d[0] = d[1] = d[2] = d[3] = 0 ;
    s[0] = s[1] = s[2] = s[3] = 0 ;
    ndims = zfp_field_dimensionality(field) ;         // get number of dimensions
    dimsize = zfp_field_size(field, dims) ;           // get dimensions
    is_strided = zfp_field_stride(field, strides) ;   // get stride information
    if(ndims>0) { d[0] = dims[0] ; s[0] = strides[0] ; }
    if(ndims>1) { d[1] = dims[1] ; s[1] = strides[1] ; }
    if(ndims>2) { d[2] = dims[2] ; s[2] = strides[2] ; }
    if(ndims>3) { d[3] = dims[3] ; s[3] = strides[3] ; }
  }else{
    d[0] = d[1] = d[2] = d[3] = 0 ;
  }

end:
  // cleanup
  if(field)  zfp_field_free(field);   // release field structure
  if(zfp)    zfp_stream_close(zfp);   // release zfp structure
  if(stream) stream_close(stream);    // release bit stream structure
  if(debug)
    fprintf(stderr, " %dD (%d,%d,%d,%d), %scontiguous (%d,%d,%d,%d)\n", 
            ndims, d[0], d[1], d[2], d[3], is_strided ? "NOT":" ", s[0], s[1], s[2], s[3]) ;
  return ndims ;

error:
  ndims = 0 ;
  goto end ;
}

int32_t ZfpExpand(void* array, int nx, int ny, int nz, void *Stream, int Ssize) {
  size_t zfpsize = Ssize ;
  bitstream* stream = stream_open(Stream, zfpsize) ;  // create bit stream structure
  zfp_stream* zfp   = zfp_stream_open(NULL) ;         // create zfp structure
  zfp_field* field  = zfp_field_alloc() ;             // create field structure
  size_t dims[4] ;
  size_t dimsize ;
  uint32_t ndims = 0 ;
  int status = 0 ;
  float ratio ;
  size_t nheader ;

  if(stream == NULL) goto error ;
  if(zfp == NULL)    goto error ;
  if(field == NULL)  goto error ;

  zfp_stream_set_bit_stream(zfp, stream) ;            // associate stream to zfp
  zfp_stream_rewind(zfp) ;                            // set stream pointer to beginning of stream
  if( (nheader = zfp_read_header(zfp, field, ZFP_HEADER_FULL)) != 0){
    if(debug) fprintf(stderr, "read header size = %ld  ", nheader);
    ndims   = zfp_field_dimensionality(field) ;       // get number of dimensions
    dimsize = zfp_field_size(field, dims) ;           // get dimensions
    if(ndims>0 && nx != dims[0]) goto error ;         // check against rquested dimensions
    if(ndims>1 && ny != dims[1]) goto error ;
    if(ndims>2 && nz != dims[2]) goto error ;
    if(ndims>3                 ) goto error ;
  }else{
    goto error ;
  }
  zfp_field_set_pointer(field, array);
  if (!zfp_decompress(zfp, field)) {
    if(diag) fprintf(stderr, "ERROR: decompression failed\n");
    goto error ;
  }
  ratio = (dimsize*4.0f) / zfpsize ;
  if(debug) fprintf(stderr, "expansion ratio = %4.2f\n", ratio) ;

end :
  // cleanup
  if(field)  zfp_field_free(field);   // release field structure
  if(zfp)    zfp_stream_close(zfp);   // release zfp structure
  if(stream) stream_close(stream);    // release bit stream structure
  return status ;

error:
  status = 1 ;
  goto end ;
}
