//  libplugin - functions for C and Fortran programming
//  Copyright (C) 2020-2022  Environnement Canada
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
#if 0
//****p* libplugin/plugins
// Synopsis
// Fortran anc C callable package to use dynamic libraries (shared objects) at runtime
//
// a "plugin" is a dynamic library (shared object) that contains some special elements
//
// http://github.com/mfvalin/wrapper-code/tree/master/PLUGINS
// 
// step 0: declarations
//         Fortran :
//           use ISO_C_BINDING
//           include 'plugins.inc'
//         Fortran with module :
//           use ISO_C_BINDING
//           use fortran_plugins
//         C:
//           #include <plugins.h>
// 
// step 1: Load a dynamic library (shared object)
//         Fortran :
//           type(C_PTR) :: handle
//           handle = load_plugin("my_library_name.so")  ! NULL pointer if not successful
//         Fortran with module :
//           type(plugin) :: sharedf
//           logical :: status     ! .true. if successful
//           status = sharedf % load('libsharedf.so')
//         C :
//           void *handle;
//           handle = load_plugin("my_library_name.so");  ! NULL pointer if not successful
// 
// step 2: Get number of advertised entry points in dynamic library
//         Fortran:
//           integer(C_INT) :: nsym
//           nsym = plugin_n_functions(handle)  ! zero if not successful
//         Fortran with module :
//           type(plugin) :: sharedf
//           nsym = sharedf % symbols()         ! zero if not successful
//         C:
//           int nsym = plugin_n_functions(handle);   ! zero if not successful
// 
// step 3: Get the name of advertised entry point n
//         Fortran:
//           type(C_PTR) :: string
//           integer(C_INT) :: n
//           string = plugin_function_name(handle,n)  ! NULL pointer if not successful
//         Fortran with module :
//           type(plugin) :: sharedf
//           logical :: status               ! .true. if successful
//           character(len=128) :: longstr   ! "" if not successful
//           status = sharedf % fname(n, longstr)
//         C:
//           int n;
//           char *string = plugin_function_name(handle,n);  ! NULL pointer if not successful
// 
// step 4: Get address of entry point by name and call it
//         Fortran:
//           type(C_FUNPTR) :: faddress
//           procedure(xxx), pointer :: fptr
//           character(C_CHAR), dimension(*) :: name
//           faddress = plugin_function(handle,name)  ! NULL pointer if not successful
//           call c_f_procpointer(faddress,fptr)
//           call fptr(...arguments...)
//         Fortran with module :
//           type(plugin) :: sharedf
//           character(len=128) :: longstr
//           faddress = sharedf % fnptr(trim(longstr))  ! NULL pointer if not successful
//           call c_f_procpointer(faddress,fptr)
//           call fptr(...arguments...)
//         C:
//           void *faddress;
//           char *name;
//           faddress = plugin_function(handle,name);  ! NULL pointer if not successful
//           .. = (*faddress)(...arguments...);
// 
// step n: unload plugin
//         Fortran:
//           integer(C_INT) :: status
//           status = unload_plugin(handle)        ! 0 if successful
//         Fortran with module :
//           type(plugin) :: sharedf
//           logical :: status                     ! .true. if successful
//           status = sharedf % unload()
//         C:
//           int status = unload_plugin(handle);   ! 0 if successful
// 
// other : set diagnostics verbosity
//         Fortran:
//           integer(C_INT) :: verbose
//           call set_plugin_diag(verbose)
//         Fortran with module :
//           type(plugin) :: sharedf
//           call sharedf % diag(VERBOSE)
//         C:
//           int verbose;
//           set_plugin_diag(verbose);
// 
// NOTES:
//   the presence of entry point "EntryList_" is mandatory in a "plugin", it is
//   a NULL pointer terminated list of pointers to NULL terminated strings
//   providing the names of the advertised entry points in the plugin
//   (see examples below)
// 
//   function "get_symbol_number" is optional and may be used to get the
//   number of values in "EntryList_" (see Fortran)
//
//   if neither "EntryList_" nor "get_symbol_number" is present, plugin_n_functions
//   will return 0, and the user must then known in advance the names available 
//   before calling "plugin_function"
// 
// EXAMPLES
// ----------------------- Example of C plugin -----------------------
// c_compiler -shared -fpic -o libxxx.so xxx.c
// 
// #include <stdio.h>
// #include <string.h>
// 
// char *EntryList_[4] = { "name1","name2",NULL} ;
// 
// int name1(int arg){
// printf("name1: %d\n",arg);
// return(arg);
// }
// 
// int name2(int arg){
// printf("name2: %d\n",arg);
// return(arg);
// }
// 
// int get_symbol_number(){  // function to get number of symbols, optional
//   return(2);
// }
// void __attribute__ ((constructor)) PluginConstructor(void) { // executed before load
//    printf("plugin constructor for plugin\n");
// }
// 
// void __attribute__ ((destructor)) PluginDestructor(void) {   // executed before unload
//    printf("plugin destructor for plugin\n");
// }
// 
// ----------------------- Example of Fortran plugin -----------------------
//               ( needs a little more extra code than C )
// fortran_compiler -shared -fpic -o libxxx.so xxx.F90 -lbuildplugin
// libbuildplugin.a is a special library that contains constructor ans destructor functions
//                  that call fortran_constructor and fortran_destructor
// 
// integer function  fn1(arg) BIND(C,name='name1f')
// integer, intent(IN) :: arg
// print *,'Fortran name1 =',arg
// fn1 = arg
// return
// end
// 
// integer function  fn2(arg) BIND(C,name='name2f')
// integer, intent(IN) :: arg
// print *,'Fortran name2 =',arg
// fn2 = arg
// return
// end
// !
// ! what follows is boiler plate code
// ! to be adjusted by user : MAX_NAMES, MAX_NAME_LENGTH, 
// !                          calls to insert_in_name_table in subroutine user_symbols
// !                          initialization (ONCE ONLY) code
// ! fortran_constructor will be called by the plugin library constructor
// ! fortran_destructor will be called by the plugin library destructor
// #define MAX_NAMES 2
// #define MAX_NAME_LENGTH 8
// #include <library_plugin_mod.hf>
// subroutine fortran_constructor() bind(C,name='fortran_constructor')
//   use library_plugin_mod
//   implicit none
// ! START of user adjustable code
//   call insert_in_name_table('name1f')
//   call insert_in_name_table('name2f')
// ! perform any tasks needed for library initialization here
// ! END of user adjustable code
//   return
// end
// subroutine fortran_destructor() bind(C,name='fortran_destructor')
// ! START of user adjustable code
// ! perform any tasks needed before closing library
// ! END of user adjustable code
// end 
// 
//****
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dlfcn.h>

#include <plugins.h>

typedef int (*fnptr)();             // pointer to function

typedef struct {
  void *handle;        // blind pointer to plugin(library) handle (from dlopen)
  char *name;          // pointer to plugin(library) name
  const char* *symbol; // pointer to list of symbol names (from plugin) (NULL if no such list)
  fnptr *addr;         // list of functions addresses (function names in symbol) (from dlsym) (NULL if no such list)
  int nentries;        // number of functions advertised in plugin (size of symbol table/address list)
  int ordinal;         // index in plugin table
} plugin;

#define MAX_PLUGINS 256
static plugin plugin_table[MAX_PLUGINS];
static int last_plugin=0;
static int verbose = 0;

// interface   !InTf!
//----------------------------------------------------------------------------------------
//****f* libplugin/Unload_plugin
// Synopsis
// unload a plugin library
//
// Fortran interface
//   function Unload_plugin(handle) result(status) BIND(C,name='Unload_plugin') !InTf!
//     import :: C_INT, C_PTR                                                   !InTf!
//     integer(C_INT) :: status                                                 !InTf!
//     type(C_PTR), intent(IN), value :: handle                                 !InTf!
//   end function Unload_plugin                                                 !InTf!
//
//  handle          : handle obtained from Load_plugin
//  function return : 0 if no error, -1 if there was an error
//
// ARGUMENTS
int Unload_plugin(void *handle)  //InTc
//****
{
  plugin *p = (plugin *) handle;
  if(p == NULL) return(-1) ; // bad handle
  if( (p - plugin_table) >= last_plugin) return(-1);  // beyond table end

  if(p->addr) free(p->addr);           // free address table if allocated
  p->addr = NULL;
  p->symbol = NULL;
  dlclose(p->handle);
  p->handle = NULL;
  p->ordinal = -1;
  p->nentries = 0;
  if(verbose) printf("INFO: plugin %s has been closed (slot %ld/%d)\n",p->name,(p-plugin_table)+1,MAX_PLUGINS);
  if(p->name) free(p->name);
  p->name = NULL;
  return (0);   // entry has been freed
}

//----------------------------------------------------------------------------------------
//****f* libplugin/Set_plugin_diag
// Synopsis
// set diagnostic verbosity level
//
// Fortran interface
//   subroutine Set_plugin_diag(level) BIND(C,name='Set_plugin_diag')  !InTf!
//     import C_INT                                                    !InTf!
//     integer(C_INT), value :: level                                  !InTf!
//   end subroutine Set_plugin_diag                                    !InTf!
//
// diag  : diagnotic verbosity level
// 0 = silent
// 1 = verbose
// there is no return value
//
// ARGUMENTS
void Set_plugin_diag(int diag)  //InTc
//****
{
  verbose = diag;
}
//----------------------------------------------------------------------------------------
//****f* libplugin/Load_plugin
// Synopsis
//  load a plugin (shared object)
//
// Fortran interface
//   function Load_plugin(lib) result(handle) BIND(C,name='Load_plugin')  !InTf!
//     import :: C_CHAR, C_INT, C_PTR                                     !InTf!
//     character(C_CHAR), dimension(*), intent(IN) :: lib                 !InTf!
//     type(C_PTR) :: handle                                              !InTf!
//   end function Load_plugin                                             !InTf!
//
//  lib     : character string, name of shared object
//
//  function return :  handle (blind pointer to a plugin structure)
//
//  the rules for finding "lib" are the rules followed by dlopen (see: man dlopen)
//
// ARGUMENTS
void *Load_plugin(const char *lib)  //InTc
//****
{
  const char **temp;
  int nsym;
  int i;
  plugin *p;
  fnptr func_nb;
  int slot ;

  for (slot=0 ; slot < MAX_PLUGINS ; slot++){
    if(plugin_table[slot].name == NULL) break;  // free slot
  }
  if(slot >= MAX_PLUGINS) {
    fprintf(stderr,"ERROR: plugin table is full\n");
    return(NULL);   // table is full
  }
// ===============================================================================================
//                                             open shared library
// ===============================================================================================
  p = &plugin_table[slot];                      // use free slot just found
  if(verbose) fprintf(stderr,"\nINFO: attempting to load shared object %s (slot %d/%d)\n",lib,slot+1,MAX_PLUGINS);
  p->ordinal = slot;
  p->name    = NULL;
  p->symbol  = NULL;
  p->nentries= 0;
  p->addr    = NULL;
  p->handle = dlopen(lib,RTLD_LAZY);

  if( p->handle == NULL ) {
    fprintf(stderr,"ERROR: load failed for %s\n",lib);
    return(NULL);
  }
  p->name = (char *)malloc(strlen(lib)+1);      // only allocate p->name if load successful
  strncpy(p->name, lib, strlen(lib)+1);         // copy name

// ===============================================================================================
//     get number of available symbols with get_symbol_number if available
//     (may be used as init call by shared library)
//     then look for EntryList_ and use it (if present) to count advertized entries
// ===============================================================================================
  nsym = 0;
  func_nb = dlsym(p->handle,"get_symbol_number");     // get_symbol_number found ?
  if(func_nb == NULL){                                // optional function (mainly for Fortran plugin usage)
    if(verbose) fprintf(stderr,"INFO: get_symbol_number not found\n");
  }else{
    nsym = (*func_nb)(); 
    if(verbose) fprintf(stderr,"INFO: %d symbols found by get_symbol_number \n", nsym);
  }

  temp = (const char **)dlsym(p->handle,"EntryList_");    // EntryList_  (mandatory for plugin, absent in regular shared object)
  if(temp == NULL){
    if(verbose) fprintf(stderr,"INFO: EntryList_ not found\n");
  }
  if(nsym == 0 && temp != NULL){                            // count names in list if get_symbol_number not found
    nsym = 0 ; while(temp[nsym]) nsym++;
    if(verbose) fprintf(stderr,"INFO: %d symbols found in EntryList_ \n", nsym);
  }

  p->symbol   = temp;                                     // symbol table is at address of symbol EntryList_
  p->nentries = nsym;                                     // number of entries advertised
  if(verbose) fprintf(stderr,"INFO: %d functions advertised in %s %s\n",p->nentries,nsym ? "plugin" : "shared object",p->name);
  if(nsym > 0) p->addr = (void *)malloc(nsym*sizeof(void *));   // allocate address table if necessary
  if(slot >= last_plugin) last_plugin = slot + 1;

  for(i=0 ; i<nsym ; i++){                                // fill address table if advertised symbol list is present
    p->addr[i] =  dlsym(p->handle,p->symbol[i]);
    if(verbose) fprintf(stderr,"INFO:   %p %s\n",p->addr[i],p->symbol[i]);
  }
  return(p);  // pointer to "plugin" structure
}
//----------------------------------------------------------------------------------------
//****f* libplugin/Plugin_n_functions
// Synopsis
//  find number of functions in a given plugin
//
// Fortran interface
//   function Plugin_n_functions(handle) result(number) BIND(C,name='Plugin_n_functions')  !InTf!
//     import :: C_PTR, C_INT                                                              !InTf!
//     type(C_PTR), intent(IN), value :: handle                                            !InTf!
//     integer(C_INT) :: number                                                            !InTf!
//   end function Plugin_n_functions                                                       !InTf!
//
//  handle          : handle obtained from Load_plugin
//
//  function return :  number of symbols advertised in the shared object, 0 if error
//
// ARGUMENTS
int Plugin_n_functions(const void *handle)      //InTc how many functions are advertised in this plugin
//****
{
  plugin *p = (plugin *) handle;

  if(p == NULL) return(0) ;                // bad handle

  if( (p - plugin_table) >= last_plugin) return(0);  // beyond table end

  return(p->nentries);                     // number of advertised names
}
//----------------------------------------------------------------------------------------
//****f* libplugin/Plugin_function_name
// Synopsis
// get name of entry number ordinal from plugin entry name table
//
// Fortran interface
//   function Plugin_function_name(handle,ordinal) result(string) BIND(C,name='Plugin_function_name')  !InTf!
//     import :: C_PTR, C_INT                                                                          !InTf!
//     type(C_PTR), intent(IN), value :: handle                                                        !InTf!
//     integer(C_INT), value :: ordinal                                                                !InTf!
//     type(C_PTR) :: string                                                                           !InTf!
//   end function Plugin_function_name                                                                 !InTf!
//
//  handle  : handle obtained from Load_plugin
//  ordinal : int, ordinal in entry list of desired name (first entry has ordinal 1)
//
//  function return : pointer to the name of requested entry ( char * / type(C_PTR)  )
//
// ARGUMENTS
const char* Plugin_function_name(const void *handle, int ordinal)  //InTc
//****
{
  plugin *p = (plugin *) handle;

  if(p == NULL) return(NULL) ; // bad handle

  if( (p - plugin_table) >= last_plugin) return(NULL);    // out of plugin table

  if(ordinal < 1 || ordinal > p->nentries || p->symbol == NULL) return(NULL);  // out of range or no symbols

  return((const char*)p->symbol[ordinal-1]);              // address of list of names
}
//----------------------------------------------------------------------------------------
//****f* libplugin/Plugin_function_names
// Synopsis
//  get list of advertised entry names in a given plugin
//
// there is no Fortran interface
//
//  handle : handle obtained from Load_plugin
//
//  function return : pointer to the list of names (char**) or NULL if none available
//
// ARGUMENTS
const char* *Plugin_function_names(const void *handle)  //InTc
//****
{
  plugin *p = (plugin *) handle;

  if(p == NULL) return(NULL) ;                          // bad handle

  if( (p - plugin_table) >= last_plugin) return(NULL);  // beyond table end

  return((const char* *)p->symbol);                     // address of list of names (can be NULL)
}
//----------------------------------------------------------------------------------------
//****f* libplugin/Plugin_function
// Synopsis
//  get address of plugin entry name
//  if pointer to plugin is NULL, scan all known plugins
//
// Fortran interface
//   function Plugin_function(handle,fname) result(faddress) BIND(C,name='Plugin_function')  !InTf!
//     import :: C_PTR, C_FUNPTR, C_CHAR                                                     !InTf!
//     type(C_PTR), intent(IN), value :: handle                                              !InTf!
//     character(C_CHAR), dimension(*), intent(IN) :: fname                                  !InTf!
//     type(C_FUNPTR) :: faddress                                                            !InTf!
//   end function Plugin_function                                                            !InTf!
//
//  handle : handle obtained from Load_plugin (if NULL, scan all known plugins and return first match)
//  name   : null terminated string, name of entry 
//
//  function return : address of requested entry ( void * / type(C_PTR)  )
//
// ARGUMENTS
void *Plugin_function(const void *handle, const char *name)  //InTc
//****
{
  plugin *p = (plugin *) handle;
  int i, j;
  void *faddr;

  if(p == NULL) {           // scan all plugins for name
    for(j=0 ; j<last_plugin ; j++){
      p = &plugin_table[j];
      for(i=0 ; i<p->nentries ; i++){
        if(strcmp(name, p->symbol[i]) == 0) return(p->addr[i]) ;
      }
      // no entries because there was no EntryList_ symbol or unadvertised symbol, try dlsym
      faddr = dlsym(p->handle,name) ;
      if(faddr != NULL) return (faddr) ;
    }
    return(NULL);   // nothing found
  }

  if( (p - plugin_table) >= last_plugin) return(NULL);  // beyond table end
  // if no entries because there was no EntryList_ symbol, dlsym will be called directly
  for(i=0 ; i<p->nentries ; i++){
    if(strcmp(name, p->symbol[i]) == 0) return(p->addr[i]) ;
  }
  return ( dlsym(p->handle,name) ); // try unadvertised symbol. if not found, NULL will be returned
}
//----------------------------------------------------------------------------------------
// end interface   !InTf!
