using PowerModels
using JSON

data = PowerModels.parse_file("data/MATPOWER/case9.m")
println(json(data, 4))
buses = [x.second for x in data["bus"] if (x.second[pm_component_status["bus"]] != pm_component_status_inactive["bus"])]
sort!(buses, by=x->x["index"])
#
idx_to_bus = [x["index"] for x in buses]
bus_to_idx = Dict(x["index"] => i for (i,x) in enumerate(buses))
for (i, branch) in data["branch"]
	f_bus = bus_to_idx[branch["f_bus"]]
	t_bus = bus_to_idx[branch["t_bus"]]
	y = inv(branch["br_r"] + branch["br_x"]im)
	tr, ti = PowerModels.calc_branch_t(branch)
	t = tr + ti*im
	lc_fr = branch["g_fr"] + branch["b_fr"]im
	lc_to = branch["g_to"] + branch["b_to"]im
	println("$f_bus --> $t_bus")
	println("	y = $y")
	println("	tr = $tr; ti = $ti")
	println("	t = $t")
	println("	lc_fr = $lc_fr")
	println("	lc_to = $lc_to")
end
