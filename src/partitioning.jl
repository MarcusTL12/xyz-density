
# r[i] = r0 + (i - 1) * Δ
# r[i] <= x < r[i + 1]
# r0 + (i - 1) * Δ <= x < r0 + i * Δ
# (i - 1) * Δ <= x - r0 < i * Δ
# i - 1 <= (x - r0) / Δ < i
function locate_domain_index(r, x)
    Int(fld(x - first(r), step(r))) + 1
end

function partitioning_naive(atoms, coords, box_min, box_max, n_domains_xyz)
    box_ranges = ((range(x0, xend, length=(n + 1))[1:end-1]
                   for (x0, xend, n) in
                   zip(box_min, box_max, n_domains_xyz))...,)

    atoms_partitioned = [Int[] for
                         _ in box_ranges[1],
                         _ in box_ranges[2],
                         _ in box_ranges[3]]

    coords_partitioned = [NTuple{3,Float64}[] for
                          _ in box_ranges[1],
                          _ in box_ranges[2],
                          _ in box_ranges[3]]

    for (atom, coord) in zip(atoms, eachcol(coords))
        ixyz = ((locate_domain_index(r, x)
                 for (r, x) in zip(box_ranges, coord))...,)

        push!(atoms_partitioned[ixyz...], atom)
        push!(coords_partitioned[ixyz...], (coord...,))
    end

    box_ranges, atoms_partitioned, coords_partitioned
end

function partitioning(atoms, coords, box_min, box_max, Δ)
    ns = ((length(x0:Δ:x1) for (x0, x1) in zip(box_min, box_max))...,)

    box_ranges = ((range(x0, xend, length=(n + 1))[1:end-1]
                   for (x0, xend, n) in
                   zip(box_min, box_max, ns))...,)

    sortable_data = [(locate_domain_index.(box_ranges, coord), atom, coord)
                     for (atom, coord) in zip(atoms, coords)]

    sort!(sortable_data)

    atoms_sorted = [atom for (_, atom, _) in sortable_data]
    coords_sorted = [coord for (_, _, coord) in sortable_data]

    atom_ranges = fill(1:0, length.(reverse(box_ranges)))
    cur_start = 1
    cur_inds = (1, 1, 1)
    for (i, (inds, _, _)) in enumerate(sortable_data)
        if inds != cur_inds
            n_counted = i - cur_start
            atom_ranges[reverse(cur_inds)...] =
                cur_start:(cur_start+n_counted-1)
            cur_start = i
            cur_inds = inds
        end
    end

    atom_ranges[reverse(cur_inds)...] = cur_start:length(atoms)

    box_ranges, atom_ranges, atoms_sorted, coords_sorted
end

function testfunclength(start, stop, step)
    n = length(start:step:stop)
    range(start, stop, length=n)
end

function map_density_partitioned(box_ranges, atom_ranges, atoms, coords,
    masses, widths, Δ)

    ns = ((length(first(r):Δ:first(r)+step(r)) for r in box_ranges)...,)

    cube_data = zeros(((ns .- 1) .* length.(box_ranges) .+ 1)...,)

    ranges = ((range(first(r), last(r) + step(r); length=l)
               for (r, l) in zip(box_ranges, size(cube_data)))...,)

    box_data = [(x, y, z)
                for x in enumerate(box_ranges[1])
                for y in enumerate(box_ranges[2])
                for z in enumerate(box_ranges[3])]

    @fastmath @inbounds Threads.@threads for (
        (box_xi, x_start),
        (box_yi, y_start),
        (box_zi, z_start)) in box_data

        x_range = range(x_start, x_start + step(box_ranges[1]), length=ns[1])
        y_range = range(y_start, y_start + step(box_ranges[2]), length=ns[2])
        z_range = range(z_start, z_start + step(box_ranges[3]), length=ns[3])

        (box_xi == length(box_ranges[1])) || (x_range = x_range[1:end-1])
        (box_yi == length(box_ranges[2])) || (y_range = y_range[1:end-1])
        (box_zi == length(box_ranges[3])) || (z_range = z_range[1:end-1])

        xi_start = (box_xi - 1) * (ns[1] - 1) + 1
        yi_start = (box_yi - 1) * (ns[2] - 1) + 1
        zi_start = (box_zi - 1) * (ns[3] - 1) + 1

        xi_range = xi_start:xi_start+ns[1]-1
        yi_range = yi_start:yi_start+ns[2]-1
        zi_range = zi_start:zi_start+ns[3]-1

        (box_xi == length(box_ranges[1])) || (xi_range = xi_range[1:end-1])
        (box_yi == length(box_ranges[2])) || (yi_range = yi_range[1:end-1])
        (box_zi == length(box_ranges[3])) || (zi_range = zi_range[1:end-1])

        atom_box_x_start = box_xi == 1 ? box_xi : box_xi - 1
        atom_box_y_start = box_yi == 1 ? box_yi : box_yi - 1
        atom_box_z_start = box_zi == 1 ? box_zi : box_zi - 1

        atom_box_x_stop = box_xi == length(box_ranges[1]) ? box_xi : box_xi + 1
        atom_box_y_stop = box_yi == length(box_ranges[2]) ? box_yi : box_yi + 1
        atom_box_z_stop = box_zi == length(box_ranges[3]) ? box_zi : box_zi + 1

        atom_box_x_range = atom_box_x_start:atom_box_x_stop
        atom_box_y_range = atom_box_y_start:atom_box_y_stop
        atom_box_z_range = atom_box_z_start:atom_box_z_stop

        for atom_box_x in atom_box_x_range,
            atom_box_y in atom_box_y_range,
            atom_box_z in atom_box_z_range

            for atom_i in atom_ranges[atom_box_z, atom_box_y, atom_box_x]
                atom = atoms[atom_i]
                x0, y0, z0 = coords[atom_i]

                a = widths[atom]
                c = masses[atom] / (a * √π)^3
                am2 = 1 / a^2

                for (ix, x) in zip(xi_range, x_range),
                    (iy, y) in zip(yi_range, y_range),
                    (iz, z) in zip(zi_range, z_range)

                    r2 = (x - x0)^2 + (y - y0)^2 + (z - z0)^2

                    cube_data[iz, iy, ix] += c * exp(-am2 * r2)
                end
            end
        end
    end

    ranges, cube_data
end
