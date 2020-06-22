
using AbstractPlotting, MakieLayout

function heatmap_slices_execution(exec::AbstractExecution, n_slices=5, resolution=(1600,1200)
                                 ; max_ts = 1000, max_xs = 1000)
    scene, layout = layoutscene(resolution=resolution)
    
    # adding timepoint slices
    soln = exec.solution
    t = soln.t
    xs = frame_xs(exec)
    pop_names = exec.simulation.model.pop_names
    pop_idxs = 1:length(pop_names)
    
    n_x, n_p, n_t = size(soln)
    step = (length(soln.t)) ÷ n_slices
    t_idxs = 2:step:length(soln.t)
    
    hm_axes = [LAxis(scene, title = "$pop_name activity", aspect=1.0) for pop_name in pop_names]
    t_idx_subs = 1:max(1,(length(t) ÷ max_ts)):length(t) 
    xs_idx_subs = 1:max(1,(length(xs) ÷ max_xs)):length(xs) 
    heatmaps = map(1:length(pop_names)) do idx_pop
        ax = hm_axes[idx_pop]
        pop_activity = cat(population.(soln.u, idx_pop)..., dims=2)
        heatmap!(ax, t[t_idx_subs], xs[xs_idx_subs], pop_activity[xs_idx_subs,t_idx_subs]')
    end
    tightlimits!.(hm_axes)
    linkaxes!(hm_axes...)
    hideydecorations!.(hm_axes[2:end])
    
    layout[2,pop_idxs] = map(pop_idxs) do pop_idx
        slices_layout = GridLayout(rowsizes=[Auto()], alignmode=Outside(10), tellheight=true)
        slice_axes = slices_layout[:h] = map(t_idxs) do t_idx
            ax = LAxis(scene, aspect=1.0, tellheight=true)
            lines!(ax, xs[xs_idx_subs], soln[xs_idx_subs,pop_idx,t_idx])
            tightlimits!(ax)
            hideydecorations!(ax)
            hidexdecorations!(ax)
            ax
        end
        linkaxes!(slice_axes...)
        slices_layout[2,1:length(t_idxs)] = [LText(scene, "t=$(round(time, digits=1))", textsize=14, tellwidth=false) for time in t[t_idxs]]
        trim!(slices_layout)
        slices_layout        
    end
    layout[1,pop_idxs] = hm_axes
    cbar = layout[1, end+1] = LColorbar(scene, heatmaps[1], label = "Activity Level")
    cbar.width = 25

    ylabel = layout[:,0] = LText(scene, "space (μm)", rotation=pi/2, tellheight=false)
    xlabel = layout[end+1,2:3] = LText(scene, "time (ms)")
    return (scene, layout)
end

function animate_execution(filename, execution::AbstractFullExecution{T,<:Simulation{T}}; fps=20, kwargs...) where T
    solution = execution.solution
    pop_names = execution.simulation.model.pop_names
    x = coordinate_axes(Simulation73.reduced_space(execution))[1]
    t = timepoints(execution)
    max_val = maximum(solution)
	min_val = minimum(solution)
    
    scene = Scene();
    time_idx_node = Node(1)
    single_pop = lift(idx -> population_timepoint(solution, 1, idx), time_idx_node)
    lines!(scene, x, single_pop)
    ylims!(scene, (min_val, max_val))
    
    record(scene, filename, 1:length(t); framerate=fps) do time_idx # TODO @views
        time_idx_node[] = time_idx
    end
end
