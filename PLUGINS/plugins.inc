 interface
   function Unload_plugin(handle) result(status) BIND(C,name='Unload_plugin')
     import :: C_INT, C_PTR
     integer(C_INT) :: status
     type(C_PTR), intent(IN), value :: handle
   end function Unload_plugin
   subroutine Set_plugin_diag(level) BIND(C,name='Set_plugin_diag')
     import C_INT
     integer(C_INT), value :: level
   end subroutine Set_plugin_diag
   function Load_plugin(lib) result(handle) BIND(C,name='Load_plugin')
     import :: C_CHAR, C_INT, C_PTR
     character(C_CHAR), dimension(*), intent(IN) :: lib
     type(C_PTR) :: handle
   end function Load_plugin
   function Plugin_n_functions(handle) result(number) BIND(C,name='Plugin_n_functions')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: handle
     integer(C_INT) :: number
   end function Plugin_n_functions
   function Plugin_function_name(handle,ordinal) result(string) BIND(C,name='Plugin_function_name')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: handle
     integer(C_INT), value :: ordinal
     type(C_PTR) :: string
   end function Plugin_function_name
   function Plugin_function(handle,fname) result(faddress) BIND(C,name='Plugin_function')
     import :: C_PTR, C_FUNPTR, C_CHAR
     type(C_PTR), intent(IN), value :: handle
     character(C_CHAR), dimension(*), intent(IN) :: fname
     type(C_FUNPTR) :: faddress
   end function Plugin_function
 end interface
