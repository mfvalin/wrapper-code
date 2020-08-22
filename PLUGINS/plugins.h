int Unload_plugin(void *h);
void Set_plugin_diag(int diag);
void *Load_plugin(const char *lib, int diag);
int Plugin_n_functions(const void *h);
const char* Plugin_function_name(const void *h, int ordinal);
const char* *Plugin_function_names(const void *h);
void *Plugin_function(const void *h, const char *name);
