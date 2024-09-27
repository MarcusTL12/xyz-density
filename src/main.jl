using Printf

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

atom_masses::Vector{Float64} = [1.00797, 4.00260, 6.941, 9.01218, 10.81, 12.011,
    14.0067, 15.9994, 18.998403, 20.179, 22.98977, 24.305, 26.98154, 28.0855,
    30.97376, 32.06, 35.453, 39.0983, 39.948, 40.08, 44.9559, 47.90, 50.9415,
    51.996, 54.9380, 55.847, 58.70, 58.9332, 63.546, 65.38, 69.72, 72.59,
    74.9216, 78.96, 79.904, 83.80, 85.4678, 87.62, 88.9059, 91.22, 92.9064,
    95.94, 98.0, 101.07, 102.9055, 106.4, 107.868, 112.41, 114.82, 118.69,
    121.75, 126.9045, 127.60, 131.30, 132.9054, 137.33, 138.9055, 140.12,
    140.9077, 144.24, 145.0, 150.4, 151.96, 157.25, 158.9254, 162.50, 164.9304,
    167.26, 168.9342, 173.04, 174.967, 178.49, 180.9479, 183.85, 186.207, 190.2,
    192.22, 195.09, 196.9665, 200.59, 204.37, 207.2, 208.9804, 209.0, 210.0,
    222.0, 223.0, 226.0254, 227.0278, 231.0359, 232.0381, 237.0482, 238.029,
    242.0, 243.0, 247.0, 247.0, 250.0, 251.0, 252.0, 255.0, 256.0, 257.0, 258.0,
    260.0, 261.0, 262.0, 262.0, 263.0, 269.0, 272.0, 277.0]

const Å2B::Float64 = 1.8897259886

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

            coords[1, i] = parse(Float64, x) * Å2B
            coords[2, i] = parse(Float64, y) * Å2B
            coords[3, i] = parse(Float64, z) * Å2B
        end

        atoms, coords
    end
end

function map_density_naive(atoms, coords, masses, widths, box_min, box_max, Δ)
    x_range = box_min[1]:Δ:box_max[1]
    y_range = box_min[2]:Δ:box_max[2]
    z_range = box_min[3]:Δ:box_max[3]

    cube_data = zeros(
        length(z_range),
        length(y_range),
        length(x_range))

    for (atom, (x0, y0, z0)) in zip(atoms, eachcol(coords))
        a = widths[atom]
        c = masses[atom] / (a * √π)^3
        am2 = 1 / a^2

        for (ix, x) in enumerate(x_range),
            (iy, y) in enumerate(y_range),
            (iz, z) in enumerate(z_range)

            r2 = (x - x0)^2 + (y - y0)^2 + (z - z0)^2

            cube_data[iz, iy, ix] += c * exp(-am2 * r2)
        end
    end

    (x_range, y_range, z_range), cube_data
end

function write_cube_file(filename, atoms, coords, ranges, cube_data,
    print_geometry=false)
    open(filename, "w") do io
        println(io, "Atomic density\n")

        if print_geometry
            @printf io "%10d" length(atoms)
        else
            @printf io "%10d" 0
        end
        for r in ranges
            @printf io " %13.6f" first(r)
        end
        @printf io "\n"

        @printf(io, "%10d %13.6f %13.6f %13.6f\n",
            length(ranges[1]), step(ranges[1]), 0.0, 0.0)
        @printf(io, "%10d %13.6f %13.6f %13.6f\n",
            length(ranges[1]), 0.0, step(ranges[2]), 0.0)
        @printf(io, "%10d %13.6f %13.6f %13.6f\n",
            length(ranges[1]), 0.0, 0.0, step(ranges[3]))

        if print_geometry
            for (atom, r) in zip(atoms, eachcol(coords))
                @printf(io, "%10d %13.6f %13.6f %13.6f %13.6f\n",
                    atom, 0.0, r...)
            end
        end

        for ix in axes(cube_data, 3), iy in axes(cube_data, 2)
            extra_newline = false
            for iz in axes(cube_data, 1)
                @printf io " %13.6e" cube_data[iz, iy, ix]
                if iz % 6 == 0
                    @printf io "\n"
                    extra_newline = false
                else
                    extra_newline = true
                end
            end
            if extra_newline
                @printf io "\n"
            end
        end
    end
end
