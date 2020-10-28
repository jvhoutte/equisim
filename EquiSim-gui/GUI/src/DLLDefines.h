#ifndef _DLLDEFINES_GUI_H_
#define _DLLDEFINES_GUI_H_

/* Cmake will define MyLibrary_EXPORTS on Windows when it
configures to build a shared library. If you are going to use
another build system on windows or create the visual studio
projects by hand you need to define MyLibrary_EXPORTS when
building a DLL on windows.
*/
// We are using the Visual Studio Compiler and building Shared libraries


#if defined (_WIN32)
  #if defined(VIS_EXPORTS)
    #define  GUI_EXPORT __declspec(dllexport)
  #else
    #define  GUI_EXPORT __declspec(dllimport)
  #endif /* VIS_EXPORTS */
#else /* defined (_WIN32) */
 #define GUI_EXPORT
#endif

#endif /* _DLLDEFINES_H_ */
