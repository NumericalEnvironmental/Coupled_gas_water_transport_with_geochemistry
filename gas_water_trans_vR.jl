####################################################################################################
#
# gas_water_trans_vR.jl - a julia-language script to model the (1) coupled advection of gas and
# water, and (2) diffusive-dispersive component transport in partially water-saturated porous media,
# and (3) geochemical reactions involving water, gas, and mineral phases
#
# model is not currently set up to handle fully water-saturated conditions
#
# model definition can be multi-dimensional, but cross-term on dispersivity tensor is not
# considered (i.e., only longitudinal dispersion is simulated for aqueous-phase transport)
#
# density-driven gas-phase (and water-phase) advection is not considered
#
####################################################################################################


### constants
const P0 = 1.013e+5 		# reference atmospheric pressure (Pa)
const g = 9.807             # gravitational acceleration (m/sec^2)
const R = 8.316				# universal gas constant (J K-1 mol-1)
const T = 298. 				# ambient temperature (25 C)
const rho_w = 1000.         # density of water (kg/m^3)
const u_w = 8.9e-4          # viscosity of water (Pa*sec)
const u_g = 1.84e-5         # gas viscosity


########################
# data types
########################


type Mineral
	names::Array{AbstractString, 1} 				# names of mineral phases
	moles::Array{Float64, 2}	 					# moles of each mineral present in each node (a num_nodes x num_minerals 2-D array)
end


type Node
	x::Float64 										# element's node point location
	y::Float64
	z::Float64
	vol::Float64 									# volume of element
	k::Float64 										# absolute permeability
	kr_w::Float64 									# relative permeability, water abd gas phases
	kr_g::Float64
	phi::Float64 									# porosity
	n::Float64 										# Van Genuchten relative permeability curve factors
	m::Float64
	alpha::Float64 									# Van Genuchten capillary pressure curve factor
	Sr::Float64 									# residual water saturation
	Q_index::Int64 									# water source index number
	Q::Float64 										# external water volumetric flux
	F_index::Int64 									# gas source index number
	F::Float64 										# external total (molar) gas mass flux
	alpha_L::Float64 								# longitudinal dispersivity
	P::Float64 										# absolute pressure
	Pc::Float64 									# capillary pressure
	S::Float64 										# water saturation (calculated)
	Xg::Array{Float64, 1} 							# mole fractions of volatile components in gas phase
	C_g::Array{Float64, 1} 							# concentration (moles/vol) of gas components in gas phase
	C_aq::Array{Float64, 1} 						# concentration (moles/vol) of components in aqueous phase
	connect_list::Array{Tuple{Int64, Int64}} 		# list of connected elements/nodes (connection number, connecting node index number)
	sigma_w::Float64 								# water and gas total conductance (updated each time step)
	sigma_g::Float64
	sigma_c::Float64 								# aqueous component transport total conductance (updated each time step)
	sigma_x::Array{Float64, 1}	 					# gas component(s) transport total conductance (updated each time step)
	tot_Qw_in::Float64 								# total flow balances for water from volume element by advection, per time step
	tot_Qw_out::Float64	
	tot_Qg_in::Float64 								# total flow balances for gas from volume element by advection, per time step	
	tot_Qg_out::Float64
end


type Connect
	node_1::Int64 									# index numbers for connected nodes
	node_2::Int64
	dx::Float64										# inter-node distance (for computing gradients)
	area::Float64 									# volume element connection interfacial area
	conduct_w::Float64 								# connection conductances for water and gas phases
	conduct_g::Float64
	conduct_c::Float64 								# connection conductances for aqueous component transport
	conduct_x::Array{Float64, 1}	 				# connection conductances for gas component(s) transport
	Qw::Float64 									# advective water flux
	Qg::Float64 									# advective gas flux	
end


type Component
	name::AbstractString 							# name of component
	volatile::Bool 									# if true, component is present in gas phase
	D::Float64 										# gaseous diffusion coefficient
	MW::Float64 									# molecular weight (for quantifying gas-phase source contribution)
	master::AbstractString 							# corresponding name of master species for PHREEQC
end


type Knobs
	gamma::Float64 									# time-stepping weighting factor for implicit solution scheme
	dt_init::Float64 								# initial, minimum, and maximum time step size
	dt_min::Float64
	dt_max::Float64
	dP_max::Float64 								# maximum change in absolute (gas) pressure, per time step
	dPc_max::Float64 								# maximum change in capillary pressure, per time step
	dt_decrease::Float64							# time step reduction and increase factors
	dt_increase::Float64
	upw_w::Float64									# upstream weighting factors for water and gas component transport
	upw_g::Float64	
	chemrxt::Bool
end

include("geochem_module.jl")						# geochem module employs some of the data types listed above; can't include it until here


#############################
# input functions
#############################


function ReadMinerals()::Mineral
	# read mineral assignments to nodes
	infile = "minerals.txt"
	data, header = readdlm(infile, '\t', header=true)
	names = header[2:end] 							# delete first string in mineral 'names', which is just the node # column label
	moles = data[:, 2:end] 								# moles = [node(row), phase(column)]
	mineral = Mineral(names, moles)
	return mineral
end


function ReadSource(source_flag::Int64)::Array{Float64, 2}
	# read source file (for gases or aqueous components) and return distribution array
	if source_flag == 1
		infile = "source_gas.txt"
	else
		infile = "source_aqueous.txt"
	end
	data = readdlm(infile, '\t', header=true)
	nrows, ncols = size(data[1])
	source = zeros(Float64, nrows, ncols-1)
	for i = 1:nrows
		for j = 2:ncols
			source[i, j-1] = Float64(data[1][i, j])
		end
	end
	return source
end


function ReadComps()
	# read in component properties, return separate items for gases and aqueous components
	comp = Component[]
	gas = Component[]
	data = readdlm("components.txt", '\t', header=true)
	for i = 1:size(data[1], 1)
		name = data[1][i, 1]
		volatile = Bool(data[1][i, 2])
		D = Float64(data[1][i, 3])
		MW = Float64(data[1][i, 4])
		master = data[1][i, 5]
		push!(comp, Component(name, volatile, D, MW, master))		
		if volatile == true
			push!(gas, Component(name, volatile, D, MW, master))
		end
	end
	num_comps = length(comp)
	num_gas = length(gas)
	println("Read chemical component properties.")
	return comp, gas, num_comps, num_gas
end


function GetKnobs()::Knobs
	# read numerical model "knobs" from file
	data = readdlm("knobs.txt", '\t', header=false)
	gamma = Float64(data[1, 2])
	dt_init = Float64(data[2, 2])
	dt_min = Float64(data[3, 2])
	dt_max = Float64(data[4, 2])
	dP_max = Float64(data[5, 2])
	dPc_max = Float64(data[6, 2])	
	dt_decrease = Float64(data[7, 2])
	dt_increase = Float64(data[8, 2])
	upw_w = Float64(data[9, 2])
	upw_g = Float64(data[10, 2])	
	chemrxt = Bool(data[11, 2])
	knobs = Knobs(gamma, dt_init, dt_min, dt_max, dP_max, dPc_max, dt_decrease, dt_increase, upw_w, upw_g, chemrxt)
	println("Read in computational knobs.")
	return knobs
end


function ReadNodes(num_comps::Int64, num_gas::Int64, S_gMW::Array{Float64,1})

	# read in nodes file and populate node type array
	node = Node[]
	num_nodes = 0
	sigma_x = zeros(Float64, num_gas) 					# initialize total gas component transport conductance array (with respect to each gas phase)
	data = readdlm("nodes.txt", '\t', header=true)
	for i = 1:size(data[1], 1)
		x = Float64(data[1][i, 1])
		y = Float64(data[1][i, 2])
		z = Float64(data[1][i, 3])
		num = Int64(data[1][i, 4])
		x_step = Float64(data[1][i, 5])	
		y_step = Float64(data[1][i, 6])	
		z_step = Float64(data[1][i, 7])		
		vol = Float64(data[1][i, 8])	
		k = Float64(data[1][i, 9])
		phi = Float64(data[1][i, 10])
		n = Float64(data[1][i, 11])
		m = 1. - 1/n
		alpha = Float64(data[1][i, 12])
		Sr = Float64(data[1][i, 13])
		Q_index = Int64(data[1][i, 14])
		Q = Float64(data[1][i, 15])
		F_index = Int64(data[1][i, 16])
		F_blended_mass = Float64(data[1][i, 17])
		F = F_blended_mass/S_gMW[F_index]
		alpha_L = Float64(data[1][i, 18])
		P = Float64(data[1][i, 19]) 					# initial condition: pressure
		S = Float64(data[1][i, 20]) 					# initial condition: water saturation
		Xg = zeros(Float64, num_gas)
		C_g = zeros(Float64, num_comps)
		C_aq = zeros(Float64, num_comps)
		for j = 1:num_gas  								# initial gas composition (mole fraction and molar concentration)
			Xg[j] = Float64(data[1][i, 20+j])
			C_g[j] = Xg[j]*P/(R*T)
		end
		for j = 1:num_comps 							# initial aqueous composition
			C_aq[j] = Float64(data[1][i, 20+num_gas+j])
		end		
		for j = 1:num
			# for each like node; connecting node list and sigma terms will be created later from connections		
			push!(node, Node(x+(j-1)*x_step, y+(j-1)*y_step, z+(j-1)*z_step, vol, k, 0., 0., phi, n, m, alpha, Sr, Q_index,
				Q, F_index, F, alpha_L, P, 0., S, copy(Xg), copy(C_g), copy(C_aq), Tuple{Int64, Int64}[], 0., 0., 0., copy(sigma_x), 0., 0., 0., 0.))
			node[end].Pc = PCap(node[end])
			num_nodes += 1
		end
	end
	
	# write to node geometry summary file (k and S also listed as proxies for properties and initial conditions)
	fname = "nodes_geometry.csv"
	csvfile = open(fname,"w")
	line_out = "node" * "," * "x" * "," * "y" * "," * "z" * "," * "volume" * "," * "permeability" *
		"," * "saturation"
	println(csvfile, line_out)	
	for i = 1:num_nodes
		line_out = string(i) * "," * string(node[i].x) * "," * string(node[i].y) * "," * string(node[i].z) *
			"," * string(node[i].vol) *	"," * string(node[i].k) * "," * string(node[i].S)
		println(csvfile,line_out)
	end
	close(csvfile)		
	println("Processed nodes.")	
	
	return node, num_nodes

end


function ReadConnects(node::Array{Node,1}, num_nodes::Int64, gas::Array{Component,1}, num_gas::Int64, knobs::Knobs)

	# read in connections file and populate node type array
	connect = Connect[]
	num_connects = 0
	data = readdlm("connects.txt", '\t', header=true)
	for i = 1:size(data[1], 1)
		node_1 = Int64(data[1][i, 1])
		node_2 = Int64(data[1][i, 2])
		num = Int64(data[1][i, 3])
		step_1 = Int64(data[1][i, 4])	
		step_2 = Int64(data[1][i, 5])	
		dx = Float64(data[1][i, 6])		
		area = Float64(data[1][i, 7])	
		for j = 1:num
			conduct_w = Conduct_w(node, node_1, node_2, dx, area)									# water-phase conductance
			conduct_g = Conduct_g(node, node_1, node_2, dx, area) 									# gas-phase conductance			
			conduct_c = Conduct_c(conduct_w, node_1, node_2, node, dx, area, knobs, num_nodes) 		# aqueous chemical transport conductance
			conduct_x = Conduct_x(node, i, j, dx, area, gas, num_gas) 								# gaseous chemical transport conductance (array of size num_gas)
			# for each like connection; note that advective fluxes are assigned later (set to zero initially, here)
			push!(connect, Connect(node_1+(j-1)*step_1, node_2+(j-1)*step_2, dx, area, conduct_w, conduct_g, conduct_c, conduct_x, 0., 0.)) 	
			num_connects += 1
		end
	end
	
	# write to connection summary file ...
	fname = "connects_geometry.csv"
	csvfile = open(fname,"w")
	line_out = "connection" * "," * "node_1" * "," * "node_2" * "," * "dx" * "," * "area"
	println(csvfile, line_out)	
	for i = 1:num_connects
		line_out = string(i) * "," * string(connect[i].node_1) * "," * string(connect[i].node_2) *
			"," * string(connect[i].dx) * "," * string(connect[i].area)
		println(csvfile,line_out)
	end
	close(csvfile)		
	println("Processed connections.")	
	
	return connect, num_connects

end


#############################
# unsaturated flow functions
#############################


function SatEff(node_i::Node)::Float64
	# effective water saturation
	return (node_i.S - node_i.Sr)/(1. - node_i.Sr)
end


function PCap(node_i::Node)::Float64
	# capillary pressure as a function of water saturation (used to calculate initial capillary pressures from saturations)
	Se = SatEff(node_i)
	return -exp(log(exp(-log(Se)/node_i.m) - 1.0)/node_i.n - log(node_i.alpha))
end


function Sat(node_i::Node)::Float64
	# water saturation as a function of capillary pressure (update after solving pressure equations)
	Se = (1 + abs(node_i.alpha * node_i.Pc)^node_i.n)^(-node_i.m)
	return Se * (1.0 - node_i.Sr) + node_i.Sr	
end


function dSdP(node_i::Node)::Float64
	# change in water saturation per change in capillary pressure
	S_prime = -(node_i.alpha*abs(node_i.Pc))^node_i.n * ((node_i.alpha*abs(node_i.Pc))^node_i.n + 1.0)^(-node_i.m-1.0) * node_i.m * node_i.n/node_i.Pc
	return S_prime * (1.0 - node_i.Sr)
end


function PermRel_w(node_i::Node)::Float64
	# water relative permeability
	Se = SatEff(node_i)
	return Se^0.5 * (1. - (1. - Se^(1.0/node_i.m))^node_i.m)^2
end


function PermRel_g(node_i::Node)::Float64
	# gas relative permeability
	Se = SatEff(node_i)
	return (1. - Se)^0.5 * (1. - Se^(1.0/node_i.m))^(2. * node_i.m)
end


######################################
# utility functions
######################################


function RemNegConcs(node::Array{Node, 1}, num_comps::Int64, num_gas::Int64)::Array{Node, 1}
	# remove any artifact negative aqueous or gas concentrations prior to PHREEQC run
	for nd in node
		for icmp = 1:num_comps
			nd.C_aq[icmp] *= (nd.C_aq[icmp]>=0.)										
		end
		for igas = 1:num_gas
			nd.C_g[igas] *= (nd.C_g[igas]>=0.)
		end			
	end
	return node
end


function GasPartition(node::Array{Node, 1}, num_gas::Int64)::Array{Node, 1}
	# calculate gas pressures and gas mole fractions in equilibrium with gas concentrations
	PP = zeros(Float64, num_gas) 											# temporary array for gas partial pressures, per node	
	for nd in node
		nd.P = 0.
		for igas = 1:num_gas 												# calculate gas partial pressures from ideal gas law
			PP[igas] = nd.C_g[igas]*R*T 									# from P = (n/V)*R*T
			nd.P += PP[igas]										
		end
		for igas = 1:num_gas 												# update gas mole fractions (for model output)
			nd.Xg[igas] = PP[igas]/nd.P
		end			
	end
	return node
end


function GetMeanGasMW(gas::Array{Component, 1}, num_gas::Int64, S_gas::Array{Float64, 2})::Array{Float64, 1}
	# compute the weighted mean molecular weight for the gas associated with each source term
	MW_gas = zeros(Float64, num_gas)				# gas molecular weight array
	for (i, gs) in enumerate(gas)
		MW_gas[i] = gs.MW
	end
	num_sources = size(S_gas, 1)
	S_gMW = zeros(Float64, num_sources)	
	for i = 1:num_sources
		S_gMW[i] = dot(MW_gas, S_gas[i, :])
	end
	return S_gMW
end


function MapNodeConnects(node::Array{Node, 1}, connect::Array{Connect, 1})::Array{Node, 1}
	# create list of node connections per each volume element; will be used to construct matrix of flow equations
	for (i, cn) in enumerate(connect)
		push!(node[cn.node_1].connect_list, (i, cn.node_2))
		push!(node[cn.node_2].connect_list, (i, cn.node_1))
	end
	return node
end


function UpdateRelPerm(node::Array{Node,1})::Array{Node, 1}
	# update relative permeabilities of both water and gases phases across all nodes
	for nd in node
		nd.kr_w = PermRel_w(nd)
		nd.kr_g = PermRel_g(nd)		
	end
	return node
end


function Conduct_w(node::Array{Node, 1}, i::Int64, j::Int64, dx::Float64, area::Float64)::Float64
	# compute/update conductance (with respect to unsaturated water flow) between two nodes
	return mean([node[i].kr_w * node[i].k, node[j].kr_w * node[j].k]) / u_w * area / dx
end


function Conduct_g(node::Array{Node, 1}, i::Int64, j::Int64, dx::Float64, area::Float64)::Float64
	# compute/update conductance (with respect to flow of the gas phase) between two nodes
	return mean([node[i].kr_g * node[i].k * node[i].P, node[j].kr_g * node[j].k * node[j].P]) / (u_g * R * T) * area / dx
end


function Conduct_c(conduct_w::Float64, i::Int64, j::Int64, node::Array{Node, 1}, dx::Float64, area::Float64, knobs::Knobs, num_nodes::Int64)::Float64
	# compute/update chemical transport conductance (with respect to the water phase) between two nodes
	dP = zeros(Float64, 2*num_nodes) 																				# empty placeholder for the WaterDarcy function
	v = WaterDarcy(conduct_w, i, j, node, dP, knobs, num_nodes) / (area * mean([node[i].phi, node[j].phi])) 		# pore velocity across interface
	D = abs(v) * mean([node[i].alpha_L, node[j].alpha_L]) 															# (longitudinal) dispersion coefficient
	return mean([node[i].phi * node[i].S, node[j].phi * node[j].S]) * D * area / dx
end


function Conduct_x(node::Array{Node, 1}, i::Int64, j::Int64, dx::Float64, area::Float64, gas::Array{Component, 1}, num_gas::Int64)::Array{Float64, 1}
	# compute/update gas-phase diffusive conductance between two nodes (returns an array of size num_gas)
	cx = zeros(Float64, num_gas)
	for k = 1:num_gas
		cx[k] = mean([node[i].phi * (1. - node[i].S), node[j].phi * (1. - node[j].S)]) * gas[k].D * area / dx
	end
	return cx
end


function UpdateTotConduct(node::Array{Node, 1}, connect::Array{Connect, 1}, num_gas::Int64)::Array{Node, 1}
	# update node total conductances (i.e., across all connections) for water and gas phases
	for nd in node
		nd.sigma_w = 0. 									# clear old values
		nd.sigma_g = 0.
		nd.sigma_c = 0.
		for i = 1:num_gas
			nd.sigma_x[i] = 0.
		end
	end
	for cn in connect
		node[cn.node_1].sigma_w += cn.conduct_w
		node[cn.node_2].sigma_w += cn.conduct_w
		node[cn.node_1].sigma_g += cn.conduct_g
		node[cn.node_2].sigma_g += cn.conduct_g	
		node[cn.node_1].sigma_c += cn.conduct_c
		node[cn.node_2].sigma_c += cn.conduct_c	
		for i = 1:num_gas
			node[cn.node_1].sigma_x[i] += cn.conduct_x[i]
			node[cn.node_2].sigma_x[i] += cn.conduct_x[i]	
		end
	end
	return node
end


function UpdateTotAdvect(node::Array{Node, 1}, connect::Array{Connect, 1})::Array{Node, 1}
	# update node total conductances (i.e., across all connections) for water and gas phases
	for nd in node
		nd.tot_Qw_in = 0. 									# clear old values
		nd.tot_Qw_out = 0.		
		nd.tot_Qg_in = 0.
		nd.tot_Qg_out = 0.		
	end
	for cn in connect 										# flow is always defined relative to node_1 to node_2

		# water
		if cn.Qw > 0.
			node[cn.node_2].tot_Qw_in += cn.Qw 				# an inflow (positive value)
			node[cn.node_1].tot_Qw_out -= cn.Qw 			# an outflow (negative value)
		else
			node[cn.node_2].tot_Qw_out += cn.Qw 			# an outflow ...
			node[cn.node_1].tot_Qw_in -= cn.Qw 				# an inflow
		end
		
		# gas
		if cn.Qg > 0.
			node[cn.node_2].tot_Qg_in += cn.Qg
			node[cn.node_1].tot_Qg_out -= cn.Qg
		else
			node[cn.node_2].tot_Qg_out += cn.Qg
			node[cn.node_1].tot_Qg_in -= cn.Qg
		end			
		
	end
	return node
end


function UpdateAdvect(node::Array{Node, 1}, connect::Array{Connect, 1}, dP::Array{Float64, 1}, knobs::Knobs, num_nodes::Int64)::Array{Connect, 1}
	# calculate advective fluxes (mean values per time step, as constrained by knobs.gamma)
	for cn in connect
		cn.Qw = WaterDarcy(cn.conduct_w, cn.node_1, cn.node_2, node, dP, knobs, num_nodes)		
		cn.Qg = GasDarcy(cn, node, dP, knobs)
	end
	return connect
end


function GasDarcy(cn::Connect, node::Array{Node, 1}, dP::Array{Float64, 1}, knobs::Knobs)::Float64
	# mean gas volumetric flux across volume element interface (from node_1 into node_2) per current time step
	i = cn.node_1
	j = cn.node_2
	tot_P_i = node[i].P + knobs.gamma*dP[i]
	tot_P_j = node[j].P + knobs.gamma*dP[j]
	return R*T*cn.conduct_g/mean([node[i].P, node[j].P]) * (tot_P_i - tot_P_j)
end


function WaterDarcy(conduct_w::Float64, i::Int64, j::Int64, node::Array{Node, 1}, dP::Array{Float64, 1}, knobs::Knobs, num_nodes::Int64)::Float64
	# mean Darcy flux across volume element interface (from node i into node j) per current time step
	# this function is structured differently than gas-phase equivalent because pore velocity is needed elsewhere in the code
	tot_h_i = node[i].P + knobs.gamma*dP[i] + node[i].Pc + knobs.gamma*dP[i+num_nodes] + rho_w*g*node[i].z
	tot_h_j = node[j].P + knobs.gamma*dP[j] + node[j].Pc + knobs.gamma*dP[j+num_nodes] + rho_w*g*node[j].z
	return conduct_w * (tot_h_i - tot_h_j)
end


######################################
# functions to write to output files
######################################


function WriteFinalState(node::Array{Node, 1}, num_nodes::Int64, gas::Array{Component, 1}, num_gas::Int64, comp::Array{Component, 1}, num_comps::Int64, minerals::Mineral)

	# model-wide output file: header
	fname = "final_state.csv"
	csvfile = open(fname,"w")
	line_out = "node" * "," * "x" * "," * "y" * "," * "z" * "," * "P" * "," * "Pc" * "," * "S"
	for i = 1:num_gas
		line_out = line_out * "," * "X_" * gas[i].name
	end
	for i = 1:num_comps
		line_out = line_out * "," * "C_" * comp[i].name
	end	
	for i = 1:length(minerals.names)
		line_out = line_out * "," * minerals.names[i]
	end		
	println(csvfile, line_out)	

	# node-by-node summary
	for i = 1:num_nodes
		line_out = string(i) * "," * string(node[i].x) * "," * string(node[i].y) * "," * string(node[i].z) * "," * string(node[i].P) * "," * string(node[i].Pc) * 
			"," * string(node[i].S)
		for j = 1:num_gas
			line_out = line_out * "," * string(node[i].Xg[j])
		end
		for j = 1:num_comps
			line_out = line_out * "," * string(node[i].C_aq[j])
		end				
		for j = 1:length(minerals.names)
			line_out = line_out * "," * string(minerals.moles[i, j])
		end			
		println(csvfile, line_out)
	end

	close(csvfile)		
	println("Wrote final model state to output file.")
	return
end	


function MonHeader(csvfile::IOStream, gas::Array{Component, 1}, num_gas::Int64, comp::Array{Component, 1}, num_comps::Int64, minerals::Mineral)
	line_out = "time" * "," * "P" * "," * "Pc" * "," * "S"
	for i = 1:num_gas
		line_out = line_out * "," * "X_" * gas[i].name
	end
	for i = 1:num_comps
		line_out = line_out * "," * "C_" * comp[i].name
	end		
	for i = 1:length(minerals.names)
		line_out = line_out * "," * minerals.names[i]
	end	
	println(csvfile, line_out)
end
	
	
function Monitor(t::Float64, node_i::Node, node_index::Int64, minerals::Mineral, num_comps::Int64, num_gas::Int64, csvfile::IOStream)	
	# append current state of monitor node to monitoring well file
	line_out = string(t) * "," * string(node_i.P) * "," * string(node_i.Pc) * "," * string(node_i.S)
	for j = 1:num_gas
		line_out = line_out * "," * string(node_i.Xg[j])
	end
	for j = 1:num_comps
		line_out = line_out * "," * string(node_i.C_aq[j])
	end		
	for j = 1:length(minerals.names)
		line_out = line_out * "," * string(minerals.moles[node_index, j])
	end		
	println(csvfile, line_out)
	return
end


############################################
# linear algebra functions: flow/pressure
############################################


function LHSmatrixFlow(connect::Array{Connect,1}, node::Array{Node,1}, knobs::Knobs, dt::Float64, num_nodes::Int64)

	# fill out the LHS of the flow equation matrix and record row-column index positions

	row_index = Int64[] 					# indexing system for sparse matrix
	col_index = Int64[]
	data = Float64[]
	
	# diagonal elements, or those elements involving the same node but phase interactions
	for (i, nd) in enumerate(node) 

		push!(row_index, i) 								
		push!(col_index, i)				
		push!(data, knobs.gamma*nd.sigma_g + nd.phi*(1.-nd.S)*nd.vol/(R*T*dt))			# gas-phase
		
		push!(row_index, i + num_nodes) 								
		push!(col_index, i + num_nodes)				
		push!(data, knobs.gamma*nd.sigma_w + nd.phi*nd.vol*dSdP(nd)/dt) 				# water-phase
		
		push!(row_index, i + num_nodes) 								
		push!(col_index, i)				
		push!(data, knobs.gamma * nd.sigma_w)	 										# influence of gas-phase on water	
		
	end	
	
	# non-diagonal elements
	for (i, nd) in enumerate(node) 
		for cn in nd.connect_list
		
			push!(row_index, i) 								
			push!(col_index, cn[2])
			push!(data, -knobs.gamma * connect[cn[1]].conduct_g) 			# gas-gas (upper left)
			
			push!(row_index, i + num_nodes) 								
			push!(col_index, cn[2] + num_nodes)
			push!(data, -knobs.gamma * connect[cn[1]].conduct_w) 			# water-water (lower right)	

			push!(row_index, i + num_nodes) 
			push!(col_index, cn[2])				
			push!(data, -knobs.gamma * connect[cn[1]].conduct_w) 			# gas-water	(lower left)
			
			# no equations for upper right, since gas mass balance does not (directly) depend on capillary pressure
			
		end
	end
	
	return data, row_index, col_index
	
end


function UpdateFlowLHS(connect::Array{Connect,1}, node::Array{Node,1}, knobs::Knobs, dt::Float64, num_nodes::Int64)

	# fill out the LHS of the flow equation matrix (row-column index positions already have been recorded)

	data = Float64[]

	# diagonal elements, or those elements involving the same-node phase interactions
	for (i, nd) in enumerate(node) 
		push!(data, knobs.gamma*nd.sigma_g + nd.phi*(1.-nd.S)*nd.vol/(R*T*dt))			# gas-phase
		push!(data, knobs.gamma*nd.sigma_w + nd.phi*nd.vol*dSdP(nd)/dt) 				# water-phase
		push!(data, knobs.gamma * nd.sigma_w)	 										# influence of gas-phase on water	
	end	
	
	# non-diagonal elements
	for (i, nd) in enumerate(node) 
		for cn in nd.connect_list
			push!(data, -knobs.gamma * connect[cn[1]].conduct_g) 			# gas-gas (upper left)
			push!(data, -knobs.gamma * connect[cn[1]].conduct_w) 			# water-water (lower right)	
			push!(data, -knobs.gamma * connect[cn[1]].conduct_w) 			# gas-water	(lower left)
		end
	end
	
	return data
	
end


function RHSvectorFlow(connect::Array{Connect,1}, node::Array{Node,1}, num_nodes::Int64)::Array{Float64, 1}

	# construct explicit flow matrix (run for each time step)
	b = zeros(Float64, 2*num_nodes)

	for i = 1:num_nodes
		b[i] = node[i].F - node[i].sigma_g * node[i].P	
		b[i+num_nodes] = node[i].Q - node[i].sigma_w * (rho_w*g*node[i].z + node[i].P + node[i].Pc)
	end	
	
	for cn in connect

		# gas balance
        b[cn.node_1] += cn.conduct_g * node[cn.node_2].P
        b[cn.node_2] += cn.conduct_g * node[cn.node_1].P
	
		# water balance
        b[cn.node_1 + num_nodes] += cn.conduct_w *
			(rho_w*g*node[cn.node_2].z + node[cn.node_2].P + node[cn.node_2].Pc)
        b[cn.node_2 + num_nodes] += cn.conduct_w *
			(rho_w*g*node[cn.node_1].z + node[cn.node_1].P + node[cn.node_1].Pc)
		
	end	
	
	return b

end


###############################################
# linear algebra functions: aqueous transport
###############################################


function LHSmatrixTrans(connect::Array{Connect,1}, node::Array{Node,1}, knobs::Knobs, dt::Float64, num_nodes::Int64)

	# fill out the LHS of the advection-dispersion equation matrix and record row-column index positions

	row_index = Int64[] 					# indexing system for sparse matrix
	col_index = Int64[]
	data = Float64[]
	
	# diagonal elements
	for (i, nd) in enumerate(node) 
		push!(row_index, i) 								
		push!(col_index, i)		
		push!(data, nd.phi*nd.S*nd.vol/dt + knobs.gamma*(nd.sigma_c - nd.tot_Qw_in + knobs.upw_w*nd.tot_Qw_in - knobs.upw_w*nd.tot_Qw_out 
		+ (node[i].tot_Qw_in + node[i].tot_Qw_out)))		
	end	

	# non-diagonal elements
	for (i, nd) in enumerate(node) 
		for cn in nd.connect_list
			push!(row_index, i) 								
			push!(col_index, cn[2])
			term = -knobs.gamma * connect[cn[1]].conduct_c
			if (((connect[cn[1]].Qw > 0.) & (connect[cn[1]].node_2 == i)) | ((connect[cn[1]].Qw < 0.) & (connect[cn[1]].node_1 == i)))
				# flow is from node j into node i: this is an inflow-type connection for i
				term += -knobs.gamma * knobs.upw_w * abs(connect[cn[1]].Qw) 			
			else
				# flow is from node i into node j: this is an outflow-type connection for i
				term += -knobs.gamma * -abs(connect[cn[1]].Qw)  * (1 - knobs.upw_w) 			
			end
			push!(data, term)
		end
	end
	
	return data, row_index, col_index
	
end


function UpdateTransAqLHS(connect::Array{Connect,1}, node::Array{Node,1}, knobs::Knobs, dt::Float64, num_nodes::Int64)

	# fill out the LHS of the advection-dispersion equation matrix and record row-column index positions

	data = Float64[]
	
	# diagonal elements
	for (i, nd) in enumerate(node) 
		push!(data, nd.phi*nd.S*nd.vol/dt + knobs.gamma*(nd.sigma_c - nd.tot_Qw_in + knobs.upw_w*nd.tot_Qw_in - knobs.upw_w*nd.tot_Qw_out 
		+ (node[i].tot_Qw_in + node[i].tot_Qw_out)))			
	end	

	# non-diagonal elements
	for (i, nd) in enumerate(node) 
		for cn in nd.connect_list
			term = -knobs.gamma * connect[cn[1]].conduct_c
			if (((connect[cn[1]].Qw > 0.) & (connect[cn[1]].node_2 == i)) | ((connect[cn[1]].Qw < 0.) & (connect[cn[1]].node_1 == i)))
				# flow is from node j into node i: this is an inflow-type connection for i
				term += -knobs.gamma * knobs.upw_w * abs(connect[cn[1]].Qw) 			
			else
				# flow is from node i into node j: this is an outflow-type connection for i
				term += -knobs.gamma * -abs(connect[cn[1]].Qw)  * (1 - knobs.upw_w) 			
			end
			push!(data, term)
		end
	end
	
	return data
end


function RHSvectorTransAq(connect::Array{Connect,1}, node::Array{Node,1}, S_aq::Array{Float64, 2}, num_nodes::Int64, icomp::Int64, knobs::Knobs)::Array{Float64, 1}

	# construct explicit advective-dispersive matrix for aqueous species chemical transport (run for each time step)
	b = zeros(Float64, num_nodes)

	for i = 1:num_nodes 			# terms written with respect to concentration in node i
		b[i] = node[i].Q*S_aq[node[i].Q_index, icomp] + node[i].C_aq[icomp]*(-node[i].sigma_c + (1-knobs.upw_w)*node[i].tot_Qw_in + knobs.upw_w*node[i].tot_Qw_out
			- (node[i].tot_Qw_in + node[i].tot_Qw_out)) 								# special provision for change in water saturation
	end	

	for cn in connect 				# terms written with respect to concentration in node j
		if cn.Qw > 0.  
			# flow is from node_1 into node_2 
			b[cn.node_1] += (cn.conduct_c + (1-knobs.upw_w)*(-cn.Qw)) * node[cn.node_2].C_aq[icomp] 			# only outflows (negative) from node i (node_1)
			b[cn.node_2] += (cn.conduct_c + knobs.upw_w*cn.Qw) * node[cn.node_1].C_aq[icomp] 					# only inflows (positive) into node i (node_2)
		else
			# flow is from node_2 into node_1
			b[cn.node_1] += (cn.conduct_c + knobs.upw_w*(-cn.Qw)) * node[cn.node_2].C_aq[icomp] 				# only inflows (positive) into node i (node_2)			
			b[cn.node_2] += (cn.conduct_c + (1-knobs.upw_w)*cn.Qw) * node[cn.node_1].C_aq[icomp] 				# only outflows (negative) from node i (node_1)
		end
	end	
	
	return b

end


###############################################
# linear algebra functions: gas transport
###############################################


function UpdateTransGasLHS(connect::Array{Connect,1}, node::Array{Node,1}, knobs::Knobs, dt::Float64, num_nodes::Int64, igas::Int64)

	# fill out the LHS of the advection-dispersion equation matrix for a gas component

	data = Float64[]
	
	# diagonal elements
	for (i, nd) in enumerate(node) 
		push!(data, nd.phi*(1.-nd.S)*nd.vol/dt + knobs.gamma*(nd.sigma_x[igas] - nd.tot_Qg_in + knobs.upw_g*nd.tot_Qg_in - knobs.upw_g*nd.tot_Qg_out
		- (node[i].tot_Qw_in + node[i].tot_Qw_out)))		
	end	
	
	# non-diagonal elements
	for (i, nd) in enumerate(node) 
		for cn in nd.connect_list
			term = -knobs.gamma * connect[cn[1]].conduct_x[igas]
			if (((connect[cn[1]].Qg > 0.) & (connect[cn[1]].node_2 == i)) | ((connect[cn[1]].Qg < 0.) & (connect[cn[1]].node_1 == i)))
				# flow is from node j into node i: this is an inflow-type connection w.r.t. node i
				term += -knobs.gamma * knobs.upw_g * abs(connect[cn[1]].Qg) 			
			else
				# flow is from node i into node j: this is an outflow-type connection w.r.t. node i
				term += -knobs.gamma * -abs(connect[cn[1]].Qg)  * (1 - knobs.upw_g) 			
			end
			push!(data, term)
		end
	end
	
	return data
	
end


function RHSvectorTransGas(connect::Array{Connect,1}, node::Array{Node,1}, S_gas::Array{Float64, 2}, num_nodes::Int64, gas::Array{Component, 1}, igas::Int64, knobs::Knobs)::Array{Float64, 1}

	# construct explicit advective-dispersive matrix for gas species chemical transport (run for each time step)
	b = zeros(Float64, num_nodes)
	
	for i = 1:num_nodes
		b[i] = node[i].F*S_gas[node[i].F_index, igas] + node[i].C_g[igas]*(-node[i].sigma_x[igas] + (1-knobs.upw_g)*node[i].tot_Qg_in + knobs.upw_g*node[i].tot_Qg_out
		+ (node[i].tot_Qw_in + node[i].tot_Qw_out))
	end	
	
	for cn in connect
		if cn.Qg > 0. 
			# flow is from node_1 into node_2
			b[cn.node_1] += (cn.conduct_x[igas] + (1-knobs.upw_g)*(-cn.Qg)) * node[cn.node_2].C_g[igas] 		# only outflows (negative) from node i (node_1)
			b[cn.node_2] += (cn.conduct_x[igas] + knobs.upw_g*cn.Qg) * node[cn.node_1].C_g[igas] 				# only inflows (positive) into node i (node_2)
		else
			# flow is from node_2 into node_1
			b[cn.node_1] += (cn.conduct_x[igas] + knobs.upw_g*(-cn.Qg)) * node[cn.node_2].C_g[igas] 				# only inflows (positive) into node i (node_2)			
			b[cn.node_2] += (cn.conduct_x[igas] + (1-knobs.upw_g)*cn.Qg) * node[cn.node_1].C_g[igas] 				# only outflows (negative) from node i (node_1)
		end			
	end	
	
	return b

end


######################################
# main script
######################################


function gas_water_trans(t_end::Float64, mon_well::Int64)

	knobs = GetKnobs() 																# read in computational parameters
	comp, gas, num_comps, num_gas = ReadComps() 									# read chemical component properties
	S_gas = ReadSource(1)															# read sources
	S_aq = ReadSource(2)
	S_gMW = GetMeanGasMW(gas, num_gas, S_gas)										# weighted average gas molecular weight (array, one per source)
	node, num_nodes = ReadNodes(num_comps, num_gas, S_gMW)							# read nodes file and organize
	node = UpdateRelPerm(node) 														# populate relative permeabilities for both water and gas phases
	connect, num_connect = ReadConnects(node, num_nodes, gas, num_gas, knobs)		# read connections file and organize
	node = MapNodeConnects(node, connect) 											# create lists of connected nodes, per volume element
	node = UpdateTotConduct(node, connect, num_gas)									# populate total conductances, per node
	minerals = ReadMinerals() 														# read mineral assignments to nodes (regardless of whether reaction modeling is actually employed)
	
	if knobs.chemrxt == true 														# initiate input_string for PHREEQC with header keyword blocks
		include("header_phreeqc.jl")
	end
	
	t = 0.
	dt = knobs.dt_init
	
	# open monitor file and write header
	fname = "monitor.csv"
	csvfile = open(fname,"w")
	MonHeader(csvfile, gas, num_gas, comp, num_comps, minerals)
	
	# set up matrix indexing system for left-hand-side of flow/pressure balance equations, "p_data" is really a placeholder
	p_data, p_row_index, p_col_index = LHSmatrixFlow(connect, node, knobs, dt, num_nodes)
	
	# set up matrix indexing system for left-hand-side of transport equations, "c_data" is really a placeholder
	# same indexing system is used for aqueous component as well as gaseous component transport equations
	c_data, c_row_index, c_col_index = LHSmatrixTrans(connect, node, knobs, dt, num_nodes)	
	
    while (t < t_end)

		### solve the flow/pressure portion ###
	
        advect_complete = false

        while advect_complete == false

            p_data = UpdateFlowLHS(connect, node, knobs, dt, num_nodes) 				# update LHS elements (new dt, new conductances)
 			A = sparse(p_row_index, p_col_index, p_data, 2*num_nodes, 2*num_nodes)		# update sparse flow equation matrix
            b = RHSvectorFlow(connect, node, num_nodes)									# construct explicit vector
			global dP = \(A, b)	 														# solve equation set (1st half of dP is gas pressure, 2nd is capillary pressure)	

            # check maximum pressure criterion (applied to both water and gas phases) at this time step size
			sum_complete = 0
			for i = 1:num_nodes
				sum_complete += (abs(dP[i]) > knobs.dP_max) + (abs(dP[i+num_nodes]) > knobs.dPc_max)
			end
			advect_complete = 1 - sign(sum_complete)
            if advect_complete == false 				# reduce time step size if pressure change criterion not satisfied
				dt *= knobs.dt_decrease
				assert(dt > knobs.dt_min)
			end
			
		end	
		
		# update advective flows; implemented here to use most current P and dP estimates
		connect = UpdateAdvect(node, connect, dP, knobs, num_nodes)
		node = UpdateTotAdvect(node, connect)
		
		### solve the advection-diffusion equation: gaseous components ###			
		for igas = 1:num_gas 																# note: LHS for gas equations is component-specific because of differing diffusion coefficients
			c_data = UpdateTransGasLHS(connect, node, knobs, dt, num_nodes, igas) 			# update LHS elements (new dt, new conductances)
			A = sparse(c_row_index, c_col_index, c_data, num_nodes, num_nodes)				# update sparse transport equation matrix
			b = RHSvectorTransGas(connect, node, S_gas, num_nodes, gas, igas, knobs)		# construct explicit vector
			dC = \(A, b)
			# update concentrations
			for (i, nd) in enumerate(node)
				nd.C_g[igas] += dC[i] 														# gas concentrations are updated here, but not pressures
			end
		end		

		### solute concentrations are computed based on mode, as indicated by knobs.chemrxt
		if knobs.chemrxt == false
			# solve the advection-dispersion equation, ignoring geochemical reactions
			c_data = UpdateTransAqLHS(connect, node, knobs, dt, num_nodes) 						# update LHS elements (new dt, new conductances)
			A = sparse(c_row_index, c_col_index, c_data, num_nodes, num_nodes)					# update sparse transport equation matrix
			for icomp = 1:num_comps
				b = RHSvectorTransAq(connect, node, S_aq, num_nodes, icomp, knobs)				# construct explicit vector
				dC = \(A, b)
				# update concentrations
				for (i, nd) in enumerate(node)
					nd.C_aq[icomp] += dC[i]
				end
			end
		else
			# call reactive transport routine to model node-by-node water chemistry using PHREEQC's mixing capability
			node = RemNegConcs(node, num_comps, num_gas)		# suppress negative concs (e.g., oscillations) that would crash PHREEQC run
			node = GasPartition(node, num_gas) 					# PHREEQC uses gas partial pressures as input, not gas concentration
			node, minerals = ManageGeochem(input_string, node, connect, comp, gas, minerals, dt)
		end

        # apply updates after concluding current time step
        t += dt 																			# increment simulation time
		for (i, nd) in enumerate(node)
			nd.Pc += dP[i+num_nodes] 														# update capillary pressure
			nd.S = Sat(nd)		 															# update water saturation implied by capillary pressure
		end
		node = GasPartition(node, num_gas) 													# update gas pressures/concentrations as a result in changes in component concentrations
		node = UpdateRelPerm(node) 															# relative permeabilities for water and gas phases
		Monitor(t, node[mon_well], mon_well, minerals, num_comps, num_gas, csvfile) 		# append to monitoring well report
		println("Time = " * "\t" * string(t))
		
		# update conductances (prior to next time step)
		for cn in connect
			cn.conduct_w = Conduct_w(node, cn.node_1, cn.node_2, cn.dx, cn.area)
			cn.conduct_g = Conduct_g(node, cn.node_1, cn.node_2, cn.dx, cn.area)
			cn.conduct_x = Conduct_x(node, cn.node_1, cn.node_2, cn.dx, cn.area, gas, num_gas) 						
		end
		if knobs.chemrxt == false 						# don't run this step for reactive transport mode (i.e., dispersive transport not considered)
			for cn in connect
				cn.conduct_c = Conduct_c(cn.conduct_w, cn.node_1, cn.node_2, node, cn.dx, cn.area, knobs, num_nodes)
			end
		end
		
		# update total conductances, per node
		node = UpdateTotConduct(node, connect, num_gas)
		
        # update time step
        dt *= knobs.dt_increase
        dt = min(dt, knobs.dt_max, t_end - t)			
		
	end

	close(csvfile) 																	# close monitor file
	WriteFinalState(node, num_nodes, gas, num_gas, comp, num_comps, minerals) 		# write final state of system to output file

	println("Done.")
	
end


### run script as gas_water_trans(t_end, mon_well)


# t_end = model end-time
# mon_well = node index number corresponding to monitor location (node no.)

gas_water_trans(86400., 1)