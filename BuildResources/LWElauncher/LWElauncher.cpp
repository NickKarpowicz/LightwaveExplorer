// This just launches the LightwaveExplorerGTK executable in /bin - exists just
// because GTK wants a non-windows-like directory structure to live inside
// This lets me have a cleaner structure to the program directory.
//

#include "framework.h"
#include "LWElauncher.h"
void setSYCLvars() {
    wchar_t loadBuffer[1024];
    DWORD envcount = GetEnvironmentVariableW(L"SYCL_CACHE_PERSISTENT", loadBuffer, 16);
    if (envcount == 0) {
        STARTUPINFO si;
        PROCESS_INFORMATION pi;

        ZeroMemory(&si, sizeof(si));
        si.cb = sizeof(si);
        ZeroMemory(&pi, sizeof(pi));
        wchar_t setSYCLpersistent[] = L"setx SYCL_CACHE_PERSISTENT 1";
        // Start the child process. 
        CreateProcess(NULL,   // No module name (use command line)
            setSYCLpersistent,        // Command line
            NULL,           // Process handle not inheritable
            NULL,           // Thread handle not inheritable
            false,          // Set handle inheritance to false
            CREATE_NO_WINDOW,              // No creation flags
            NULL,           // Use parent's environment block
            NULL,           // Use parent's starting directory 
            &si,            // Pointer to STARTUPINFO structure
            &pi);           // Pointer to PROCESS_INFORMATION structure

        // Wait until child process exits.
        WaitForSingleObject(pi.hProcess, INFINITE);

        // Close process and thread handles. 
        CloseHandle(pi.hProcess);
        CloseHandle(pi.hThread);
    }
}

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);
    setSYCLvars();

    STARTUPINFO si;
    PROCESS_INFORMATION pi;

    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));
    wchar_t launchLWE[] = L"bin\\LightwaveExplorerGTK.exe";
    // Start the child process. 
    CreateProcess(NULL,   // No module name (use command line)
        launchLWE,        // Command line
        NULL,           // Process handle not inheritable
        NULL,           // Thread handle not inheritable
        false,          // Set handle inheritance to false
        CREATE_NO_WINDOW,              // No creation flags
        NULL,           // Use parent's environment block
        NULL,           // Use parent's starting directory 
        &si,            // Pointer to STARTUPINFO structure
        &pi);           // Pointer to PROCESS_INFORMATION structure

    // Wait until child process exits.
    WaitForSingleObject(pi.hProcess, INFINITE);
    ExitProcess(0);

    return 0;
}



