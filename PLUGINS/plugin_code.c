//  RMNLIB - useful routines for C and FORTRAN programming
//  Copyright (C) 2020  Environnement Canada
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dlfcn.h>

#include <plugins.h>

typedef int (*fnptr)();             // pointer to function

typedef struct {
  void *handle;       // blind pointer to plugin(library) handle (from dlopen)
  char *name;         // pointer to plugin(library) name
  const char* *symbol;    // pointer to list of symbol names (from plugin) (NULL if no such list)
  fnptr *addr;        // list of functions addresses (function names in symbol) (from dlsym) (NULL if no such list)
  int nentries;       // number of functions advertised in plugin (size of symbol table/address list)
  int ordinal;        // index in plugin table
} plugin;

#define MAX_PLUGINS 256
static plugin plugin_table[MAX_PLUGINS];
static int last_plugin=0;
static int verbose = 0;

// interface   !InTf!
//----------------------------------------------------------------------------------------
//   function Unload_plugin(handle) result(status) BIND(C,name='Unload_plugin') !InTf!
//     import :: C_INT, C_PTR                                                   !InTf!
//     integer(C_INT) :: status                                                 !InTf!
//     type(C_PTR), value :: handle                                             !InTf!
//   end function Unload_plugin                                                 !InTf!
//
// unload a plugin library
//
//  h  : handle obtained from Load_plugin
//  function return : 0 if no error, -1 if there was an error
//
int Unload_plugin(void *h)  //InTc
{
  plugin *p = (plugin *) h;
  if(p == NULL) return(-1) ; // bad handle
  if( (p - plugin_table) >= last_plugin) return(-1);  // out of table

  if(p->addr) free(p->addr);           // free address table if allocated
  p->addr = NULL;
  p->symbol = NULL;
  dlclose(p->handle);
  p->handle = NULL;
  p->ordinal = -1;
  p->nentries = 0;
  if(verbose) printf("INFO: plugin %s has been closed (slot %ld)\n",p->name,p-plugin_table);
  if(p->name) free(p->name);
  p->name = NULL;
  return (0);   // entry has been freed
}

//----------------------------------------------------------------------------------------
//   subroutine Set_plugin_diag(level) BIND(C,name='Set_plugin_diag')  !InTf!
//     import C_INT                                                    !InTf!
//     integer(C_INT), value :: level                                  !InTf!
//   end subroutine Set_plugin_diag                                    !InTf!
// set diagnostic verbosity level
// 0 = silent
// 1 = verbose
// there is no return
void Set_plugin_diag(int diag)  //InTc
{
  verbose = diag;
}
//----------------------------------------------------------------------------------------
//   function Load_plugin(plugin_name) result(handle) BIND(C,name='Load_plugin')  !InTf!
//     import :: C_CHAR, C_INT, C_PTR                                             !InTf!
//     character(C_CHAR), dimension(*), intent(IN) :: plugin_name                 !InTf!
//     type(C_PTR) :: handle                                                      !InTf!
//   end function Load_plugin                                                     !InTf!
//
//  load a plugin (shared object)
//
//  lib     : character string, name of shared object
//  verbose : int, if non zero some diagnostics are printed
//
//  function return :  handle (pointer to a plugin structure)
//
//  the rules for finding "lib" are the rules followed by dlopen (see: man dlopen)
//
void *Load_plugin(const char *lib, int diag)  //InTc
{
  const char **temp;
  int nsym;
  int i;
  plugin *p;
  fnptr func_nb;
  int slot ;

  verbose = diag;
  for (slot=0 ; slot < MAX_PLUGINS ; slot++){
    if(plugin_table[slot].name == NULL) break;  // free slot
  }
  if(slot >= MAX_PLUGINS) {
    fprintf(stderr,"ERROR: plugin table is full\n");
    return(NULL);   // table is full
  }

  p = &plugin_table[slot];
  p->ordinal = slot;

  p->name = (char *)malloc(strlen(lib)+1);
  strncpy(p->name, lib, strlen(lib)+1);
  if(verbose) fprintf(stderr,"\nINFO: attempting to load plugin %s (slot %d)\n",p->name,slot);
  p->handle = dlopen(p->name,RTLD_LAZY);
  if( p->handle == NULL ) {
    fprintf(stderr,"ERROR: load plugin failed for %s\n",p->name);
    free(p->name);
    p->name = NULL;
    return(NULL);
  }

  func_nb = dlsym(p->handle,"get_symbol_number");         // provide get_symbol_number (may initialize entry_list) (optional)
  temp = (const char **)dlsym(p->handle,"entry_list");    // provide **entry_list  (mandatory for plugin, absent in regular shared object)
  nsym = 0;
  if(func_nb) nsym = (*func_nb)();                        // optional function (mainly for Fortran usage)

  if( (temp == NULL) && (func_nb == NULL)) {
    if(verbose) fprintf(stderr,"WARNING: no function table found in plugin %s\n",p->name);
    dlclose(p->handle);     // nothing useful, close
    free(p->name);
    p->name = NULL;
    return(NULL);
  }

  if(nsym==0){
    nsym = 0 ; while(temp[nsym]) nsym++;
  }

  if (nsym == 0) {
    if(verbose) fprintf(stderr,"WARNING: no functions advertised in plugin %s\n",p->name);
    dlclose(p->handle);     // nothing useful, close
    free(p->name);
    p->name = NULL;
    return(NULL);     // no useful entry points found
  }

  p->symbol = temp;                                       // symbol table is at address of symbol entry_list
  p->nentries = nsym;
  if(verbose) fprintf(stderr,"INFO: %d functions advertised in plugin %s\n",p->nentries,p->name);
  p->addr = (void *)malloc(nsym*sizeof(void *));   // allocate address table
  if(slot >= last_plugin) last_plugin = slot + 1;

  for(i=0 ; i<nsym ; i++){
    p->addr[i] =  dlsym(p->handle,p->symbol[i]);   // fill address table
    if(verbose) fprintf(stderr,"INFO:   %p %s\n",p->addr[i],p->symbol[i]);
  }
  return(p);
}
//----------------------------------------------------------------------------------------
//   function Plugin_n_functions(handle) result(number) BIND(C,name='Plugin_n_functions')  !InTf!
//     import :: C_PTR, C_INT                                                              !InTf!
//     type(C_PTR), value :: handle                                                        !InTf!
//     integer(C_INT) :: number                                                            !InTf!
//   end function Plugin_n_functions                                                       !InTf!
//
//  find number of functions in a given plugin
//
//  h  : handle obtained from Load_plugin
//
//  function return :  number of symbols advertised in the shared object, 0 if error
//
int Plugin_n_functions(const void *h)      //InTc how many functions are defined in this plugin
{
  plugin *p = (plugin *) h;
  if(p == NULL) return(0) ; // bad handle
  if( (p - plugin_table) >= last_plugin) return(0);  // out of table

  return(p->nentries);           // address of list of names
}
//----------------------------------------------------------------------------------------
//   function Plugin_function_name(handle,ordinal) result(string) BIND(C,name='Plugin_function_name')  !InTf!
//     import :: C_PTR, C_INT                                                                          !InTf!
//     type(C_PTR), value :: handle                                                                    !InTf!
//     integer(C_INT), value :: ordinal                                                                !InTf!
//     type(C_PTR) :: string                                                                           !InTf!
//   end function Plugin_function_name                                                                 !InTf!
//
// get name of entry number ordinal from plugin entry name table
//
//  h       : handle obtained from Load_plugin
//  ordinal : int, ordinal in entry list of desired name (first entry has ordinal 1)
//
//  function return : pointer to the name of requested entry (char*)
//
const char* Plugin_function_name(const void *h, int ordinal)  //InTc
{
  plugin *p = (plugin *) h;
  if(p == NULL) return(NULL) ; // bad handle
  if( (p - plugin_table) >= last_plugin) return(NULL);  // out of plugin table

  if(ordinal<1 || ordinal>p->nentries) return(NULL);  // out of range for this plugin

  return((const char*)p->symbol[ordinal-1]);           // address of list of names
}
//----------------------------------------------------------------------------------------
//
//  get list of advertised entry names in a given plugin
//
//  p : handle obtained from Load_plugin
//
//  function return : pointer to the list of names (char**)
//
const char* *Plugin_function_names(const void *h)  //InTc
{
  plugin *p = (plugin *) h;
  if(p == NULL) return(NULL) ; // bad handle
  if( (p - plugin_table) >= last_plugin) return(NULL);  // out of table

  return((const char* *)p->symbol);           // address of list of names
}
//----------------------------------------------------------------------------------------
//   function Plugin_function(handle,fname) result(faddress) BIND(C,name='Plugin_function')  !InTf!
//     import :: C_PTR, C_FUNPTR, C_CHAR                                                     !InTf!
//     type(C_PTR), value :: handle                                                          !InTf!
//     character(C_CHAR), dimension(*), intent(IN) :: fname                                  !InTf!
//     type(C_FUNPTR) :: faddress                                                            !InTf!
//   end function Plugin_function                                                            !InTf!
//
//  get address of plugin entry name
//  if pointer to plugin is NULL, scan all known plugins
//
//  h    : handle obtained from Load_plugin
//  name : char *, null terminated string, name of entry 
//
//  function return : address of requested entry
//
void *Plugin_function(const void *h, const char *name)  //InTc
{
  plugin *p = (plugin *) h;
  int i, j;
  void *faddr;

  if(p == NULL) {           // scan all plugins for name
    for(j=0 ; j<last_plugin ; j++){
      p = &plugin_table[j];
      for(i=0 ; i<p->nentries ; i++){
        if(strcmp(name, p->symbol[i]) == 0) return(p->addr[i]) ;
      }
      // no entries because there was no entry_list symbol or unadvertized symbol, try dlsym
      faddr = dlsym(p->handle,name) ;
      if(faddr != NULL) return (faddr) ;
    }
    return(NULL);   // nothing found
  }

  if( (p - plugin_table) >= last_plugin) return(NULL);  // out of table
  // if no entries because there was no entry_list symbol, dlsym will be called directly
  for(i=0 ; i<p->nentries ; i++){
    if(strcmp(name, p->symbol[i]) == 0) return(p->addr[i]) ;
  }
  return ( dlsym(p->handle,name) ); // try unadvertised symbol. if not found, NULL will be returned
}
//----------------------------------------------------------------------------------------
// end interface   !InTf!
