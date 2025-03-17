var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = PyThermo","category":"page"},{"location":"#PyThermo","page":"Home","title":"PyThermo","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for PyThermo.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Species\nMixture","category":"page"},{"location":"#PyThermo.Species","page":"Home","title":"PyThermo.Species","text":"Species(ID::String; T=298.15, P=101325)\n\nCreates a Species object which contains basic information such as molecular weight and the structure of the species, as well as thermodynamic and transport properties as a function of temperature and pressure. If T and P are given as non-Unitful numbers, they must have units of K and Pa.\n\nParameters\n\nID : One of the following [-]:\n\nName, in IUPAC form or common form or a synonym registered in PubChem \nInChI name, prefixed by \"InChI=1S/\" or \"InChI=1/\" \nInChI key, prefixed by \"InChIKey=\" \nPubChem CID, prefixed by \"PubChem=\" \nSMILES (prefix with \"SMILES=\" to ensure smiles parsing) \nCAS number\n\nT : temperature of the chemical (default 298.15 K) \n\nP : pressure of the chemical (default 101325 Pa)\n\nExamples\n\njulia> He = Species(\"He\")\nSpecies(He, 298.1 K, 1.013e+05 Pa)\n\njulia> density(He)\n0.16360253235815483 kg m^-3\n\njulia> using Unitful\n\njulia> He.T = 30u\"K\"\n30 K\n\njulia> density(He)\n1.623503007493497 kg m^-3\n\nA wide variety of unexported properties can be accessed from the underlying Python object:\n\njulia> SF6 = Species(\"SF6\", P=30u\"psi\")\nSpecies(SF6, 298.1 K, 2.068e+05 Pa)\n\njulia> SF6.MW\n146.055419\n\nhelp?> SF6.<tab>\nA                                 SGs                                __init_subclass__\nAPI                               STEL                               __le__\nAm                                STEL_source                        __lt__\nBond                              STEL_sources                       __module__\nBvirial                           S_dep_Tb_P_ref_g                   __ne__\nCAS                               S_dep_Tb_Pb_g                      __new__\nCapillary                         S_dep_Tb_Pb_l                      __reduce__\nCarcinogen                        S_dep_ref_g                        __reduce_ex__\nCarcinogen_source                 S_int_Tb_to_T_ref_g                __repr__\nCarcinogen_sources                S_int_l_Tm_to_Tb                   __setattr__\n[...]\n\n\n\n\n\n","category":"type"},{"location":"#PyThermo.Mixture","page":"Home","title":"PyThermo.Mixture","text":"Mixture(chemnames::Vector{String}; kwargs...)\nMixture(chemnames::Vector{Pair{String, Float64}}, kwargs...)\n\nCreates a Mixture object which contains basic information such as molecular weight and the structure of the species, as well as thermodynamic and transport properties as a function of temperature and pressure.\n\nThe components of the mixture are specified by the names of the chemicals; the composition can be specified by providing any one of the following parameters as a keyword argument:\n\nMass fractions ws\nMole fractions zs\nLiquid volume fractions (based on pure component densities) Vfls\nGas volume fractions (based on pure component densities) Vfgs\n\nThe composition can also be specified by providing a vector of \"ID\" => molefrac pairs.\n\nExamples\n\njulia> air = Mixture([\"N2\" => 0.78, \"O2\" => 0.21, \"Ar\" => 0.01])\nMixture(78% N2, 21% O2, 1% Ar, 298.1 K, 1.013e+05 Pa)\n\njulia> soundspeed(air)\n346.1345044609487 m s^-1\n\n\n\n\n\n","category":"type"},{"location":"#Interaction-with-Conda","page":"Home","title":"Interaction with Conda","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PyThermo's Python dependencies are managed by CondaPkg.jl, which registers a project's dependencies in CondaPkg.toml (similar to Julia's Project.toml). These dependencies are installed automatically in a shared Conda environment located at ~/.julia/conda_environments/Thermo when PyThermo is first loaded. If you'd like to use a different Conda environment, you can set the corresponding preference as described in the CondaPkg.jl documentation.","category":"page"}]
}
