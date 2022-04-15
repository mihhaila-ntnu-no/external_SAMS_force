call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" --force x64
cl /D "WIN" /LD print_state.c /link /def:print_state.def /out:print_state.dll
copy /Y print_state.dll "C:\Users\mihha\SIMA Workspaces\Workspace_print_state\StorageTask"
cmd /k