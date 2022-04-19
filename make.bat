call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" --force x64
cl /D "WIN" /LD simo_extfunc.c csvwriter.c getline.c /link /def:external_SAMS_force.def /out:external_SAMS_force.dll
copy /Y external_SAMS_force.dll "C:\Users\mihha\SIMA Workspaces\Workspace_unitcube\StorageTask"
cmd /k


