# global PyObject constants that get initialized at runtime.  We
# initialize them here (rather than via "global foo = ..." in __init__)
# so that their type is known at compile-time.

const PY_THERMO = PythonCall.pynew()
const PY_CHEM = PythonCall.pynew()
const PY_COPY = PythonCall.pynew()

function __init__()
    ccall(:jl_generating_output, Cint, ()) == 1 && return nothing
    PythonCall.pycopy!(PY_THERMO, pyimport("thermo"))
    PythonCall.pycopy!(PY_CHEM, PY_THERMO.chemical)
    PythonCall.pycopy!(PY_COPY, pyimport("copy"))
end
