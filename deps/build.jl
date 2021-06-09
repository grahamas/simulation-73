@info "Active project name: $(Pkg.project().name)"

dev_path(dirs...) = joinpath(Pkg.devdir(), dirs...)

using GitCommand

Pkg.develop("AxisIndices")
@warn "Checking out custom dev copy of AxisIndices [FIXME: might be updated]"
git() do git
    axisindices_path = dev_path("AxisIndices")
    try
        run(`$git -C $axisindices_path remote add grahamas git@github.com:grahamas/AxisIndices.jl.git`)
    catch e
        @show e
    end
	run(`$git -C $axisindices_path fetch grahamas`)
	run(`$git -C $axisindices_path checkout logical_indexing`)
end


