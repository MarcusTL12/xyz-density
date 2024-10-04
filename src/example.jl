
include("main.jl")

function example1()
    atoms, coords = @time parse_xyz("testdata/tmp.xyz")

    box_ranges, atom_ranges, atoms_sorted, coords_sorted = @time partitioning(
        atoms, coords, (-144, -144, -144), (144, 144, 144), 2.0
    )

    masses = [18.0, 1.0]
    widths = [1.0, 1.0]

    ranges, cube_data = @time map_density_partitioned(
        box_ranges, atom_ranges, atoms_sorted, coords_sorted,
        masses, widths, 1.0
    )

    @time write_cube_file("testdata/test.cube", ranges, cube_data)
end
