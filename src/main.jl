
periodic_table = """
H                                                                                            He
Li Be                                                                         B  C  N  O  F  Ne
Na Mg                                                                         Al Si P  S  Cl Ar
K  Ca                                           Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
Rb Sr                                           Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
Fr Ra Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
"""

sym_to_int::Dict{String,Int} =
    Dict([s => i for (i, s) in enumerate(eachsplit(periodic_table))])

function parse_xyz(filename)
    open(filename) do io
        lines = Iterators.Stateful(eachline(io))

        n = parse(Int, popfirst!(lines))
        popfirst!(lines)

        atoms = zeros(Int, n)
        coords = zeros(3, n)

        for (i, l) in zip(1:n, lines)
            atom, x, y, z = eachsplit(l)

            atoms[i] = if haskey(sym_to_int, atom)
                sym_to_int[atom]
            else
                parse(Int, atom)
            end

            coords[1, i] = parse(Float64, x)
            coords[2, i] = parse(Float64, y)
            coords[3, i] = parse(Float64, z)
        end

        atoms, coords
    end
end
