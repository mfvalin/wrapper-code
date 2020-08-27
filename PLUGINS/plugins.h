int Unload_plugin(void *handle);
void Set_plugin_diag(int diag);
void *Load_plugin(const char *lib);
int Plugin_n_functions(const void *handle);
const char* Plugin_function_name(const void *handle, int ordinal);
const char* *Plugin_function_names(const void *handle);
void *Plugin_function(const void *handle, const char *name);
