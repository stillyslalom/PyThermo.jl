# global PyObject constants that get initialized at runtime.  We
# initialize them here (rather than via "global foo = ..." in __init__)
# so that their type is known at compile-time.

const PY_THERMO = PyNULL()
const PY_CHEM = PyNULL()
const PY_COPY = PyNULL()

function __init__()
    copy!(PY_THERMO, pyimport_conda("thermo", "thermo", "conda-forge"))
    copy!(PY_CHEM, pyimport("thermo.chemical"))
    copy!(PY_COPY, pyimport("copy"))
end