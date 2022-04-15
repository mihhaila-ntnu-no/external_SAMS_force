#if defined WIN
#define CAL_CONV __stdcall
#else
#define CAL_CONV
#endif

#if defined WIN
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif