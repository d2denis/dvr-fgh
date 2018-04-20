#include <Python.h>

int main(int argc, char *argv[])
{
PyObject *presult = NULL;
PyObject *pName = NULL;
PyObject *pModule = NULL;
PyObject *pFunc = NULL;

    if (argc < 3)
    {
        printf("Usage: exe_name python_source function_name\n");
        return 1;
    }

    // Initialize the Python Interpreter
    Py_Initialize();
PySys_SetArgv(argc,argv);
if ((pName = PyString_FromString(argv[1]))) {
    if ((pModule = PyImport_Import(pName))) {
        if ((pFunc = PyObject_GetAttrString(pModule, argv[2]))) {
            if ((presult = PyObject_CallObject(pFunc, NULL))) {
                printf("Result is %d\n", (int) PyInt_AsLong(presult));
            }
        }
    }
}
    // Clean up
    Py_DECREF(pModule);
    Py_DECREF(pName);

    // Finish the Python Interpreter
    Py_Finalize();

    return 0;
}

