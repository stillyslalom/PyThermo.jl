# global PyObject constants that get initialized at runtime.  We
# initialize them here (rather than via "global foo = ..." in __init__)
# so that their type is known at compile-time.

const PY_THERMO = PythonCall.pynew()
const PY_CHEM = PythonCall.pynew()
const PY_COPY = PythonCall.pynew()

# If the user has set JULIA_CONDAPKG_OFFLINE, use that value.  Otherwise,
# if .CondaPkg exists in the current environment, use "true" to skip
# Conda pkg resolution.  Otherwise, use "false" to allow Conda to
# resolve the pkg environment.
get!(ENV, "JULIA_CONDAPKG_OFFLINE") do
    if ispath(joinpath(dirname(Pkg.project().path), ".CondaPkg"))
        "true"
    else
        "false"
    end
end

function __init__()
    PythonCall.pycopy!(PY_THERMO, pyimport("thermo"))
    PythonCall.pycopy!(PY_CHEM, PY_THERMO.chemical)
    PythonCall.pycopy!(PY_COPY, pyimport("copy"))
end
