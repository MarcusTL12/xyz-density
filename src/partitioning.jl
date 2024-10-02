
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

function partitioning(atoms, coords, box_min, box_max, n_domains_xyz)
    box_ranges = ((range(x0, xend, length=(n + 1))[1:end-1]
                   for (x0, xend, n) in
                   zip(box_min, box_max, n_domains_xyz))...,)

    sortable_data = [(locate_domain_index.(box_ranges, coord), atom, coord)
                     for (atom, coord) in zip(atoms, coords)]

    sort!(sortable_data)

    atoms_sorted = [atom for (_, atom, _) in sortable_data]
    coords_sorted = [coord for (_, _, coord) in sortable_data]

    atom_ranges = fill(1:0, length.(box_ranges))
    cur_start = 1
    cur_inds = (1, 1, 1)
    for (i, (inds, _, _)) in enumerate(sortable_data)
        if inds != cur_inds
            n_counted = i - cur_start
            atom_ranges[cur_inds...] = cur_start:(cur_start+n_counted-1)
            cur_start = i
            cur_inds = inds
        end
    end

    atom_ranges[cur_inds...] = cur_start:length(atoms)

    atom_ranges, atoms_sorted, coords_sorted
end
