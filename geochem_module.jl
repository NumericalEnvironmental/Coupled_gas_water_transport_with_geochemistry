# load phreeqc COM object to allow access under windows; run via python calls
using PyCall
@pyimport win32com.client as opsys 		
dbase = opsys.Dispatch("IPhreeqcCOM.Object") 	
dbase[:LoadDatabase]("phreeqc.dat")

function ManageGeochem(input_string::AbstractString, node::Array{Node, 1}, connect::Array{Connect, 1}, comp::Array{Component, 1}, gas::Array{Component, 1}, minerals::Mineral, dt::Float64)

	# function steps through each node and conducts a PHREEQC equilibration and mixing model w.r.t. surrounding nodes
	C_matrix_aq = zeros(Float64, length(node), length(comp)) 		# create temporary 2-D arrays to hold new results
	C_matrix_gas = zeros(Float64, length(node), length(gas))

	# set up temporary copies of aqueous composition and gases for each node

	for i = 1:length(node)
	
		if node[i].vol <= 1e+20 									# do not model any reactions for fixed pressure- and fixed-concentration nodes
	
			# set up PHREEQC input and run
			input_string = WriteMixingInput(input_string, i, node, connect, comp, gas, minerals, dt)
			dbase[:RunString](input_string)
			phreeqc_result = dbase[:GetSelectedOutputArray]()
			
			# process PHREEQC output and update concentrations
			header = [] 											# process output column headings to allow search for aqueous components, gases, and minerals
			for name in phreeqc_result[1]
				push!(header, replace(name, "(mol/kgw)", ""))
			end
			for (j, cmp) in enumerate(comp)
				indexnum = findfirst(header, cmp.master)
				C_matrix_aq[i, j] = phreeqc_result[end][indexnum]
			end
			for (j, gs) in enumerate(gas)
				indexnum = findfirst(header, "g_" * gs.name * "(g)")
				C_matrix_gas[i, j] = phreeqc_result[end][indexnum] * 1000.0 		# multiplication by 1000 converts mol/L to mol/m^3			
			end
			for (j, ms) in enumerate(minerals.names)
				indexnum = findfirst(header, ms)
				minerals.moles[i, j] = phreeqc_result[end][indexnum]
			end
		
		end
		
	end

	# update aqueous component and gas concentrations
	for (i, nd) in enumerate(node)
		if nd.vol <= 1e+20
			for j= 1:length(comp)
				nd.C_aq[j] = C_matrix_aq[i, j]
			end
			for j = 1:length(gas)
				nd.C_g[j] = C_matrix_gas[i, j]
			end
		end
	end	
	
	return node, minerals
	
end


function WriteSolution(input_string::AbstractString, sol_num::Int64, sol_node::Node, comp::Array{Component, 1}, gas::Array{Component, 1})::AbstractString
	# add to solutions definition in PHREEQC input
	input_string *= "SOLUTION" * "\t" * string(sol_num) * "\n"
	input_string *= "temp" * "\t" * string(T-273) * "\n" 				# convert ambient temperature to C from K
	input_string *= "redox" * "\t" * "pe" * "\n"
	input_string *= "units" * "\t" * "mol/kgw" * "\n"
	input_string *= "density" * "\t" * "1.0" * "\n"	
	for (j, cmp) in enumerate(comp)
		input_string *= cmp.master * "\t" * string(sol_node.C_aq[j]) * "\n"
	end
	input_string *= "-water" * "\t" * "1.0 # kg" * "\n" 		
	return input_string
end


function WriteMixingInput(input_string::AbstractString, node_index::Int64, node::Array{Node, 1}, connect::Array{Connect, 1}, comp::Array{Component, 1}, gas::Array{Component, 1}, minerals::Mineral, dt::Float64)::AbstractString

	# this function prepares the input string that will inform a complete instance of PHREEQC for node[node_index]

	# write out solution compositions
	input_string = WriteSolution(input_string, 1, node[node_index], comp, gas) 			# composition of solution associated with node_index
	for (i, cn) in enumerate(node[node_index].connect_list)								# cn[1] is the connection number; cn[2] is the index number of the connecting node
		input_string = WriteSolution(input_string, i+1, node[cn[2]], comp, gas)
	end

	# specify mixing calculations
	vol_w = [node[node_index].S * node[node_index].phi * node[node_index].vol] 			# volume of water in node_index at start of time step
	for (i, cn) in enumerate(node[node_index].connect_list)
		if connect[cn[1]].node_1 == node_index 				# flow is defined as out of node_index
			Q_in = -connect[cn[1]].Qw * dt
		else
			Q_in = connect[cn[1]].Qw * dt
		end
		if Q_in < 0. 									# only influx water is used in mixing calcs; outflux is already mixed!
			Q_in = 0.
		end
		push!(vol_w, Q_in)
	end
	sum_water = sum(vol_w) 							# to normalize initial water volume + influx volumes to 1 kg/1 L
	input_string *= "MIX 1" * "\n"
	for i = 1:length(vol_w)
		input_string *= "\t" * string(i) * "\t" * string(vol_w[i]/sum_water) * "\n"
	end

	# equilibrate with equilibrium phases
	input_string *= "EQUILIBRIUM_PHASES 1" * "\n"
	for (i, phase) in enumerate(minerals.names)
		input_string *= "\t" * phase * "\t" * string(0) * "\t" * string(minerals.moles[node_index, i]) * "\n"
	end
	
	# equilibrate with gas phases
	P_atm = node[node_index].P/P0 									# total pressure in node, units = atm (at start of time step, not updated yet)
	S_g = 1. - node[node_index].S 									# gas-phase saturation
	if S_g > 0.
		input_string *= "GAS_PHASE 1" * "\n"	
		input_string *= "-fixed_volume" * "\n"
		input_string *= "-pressure" * "\t" * string(P_atm) * "\n"		
		input_string *= "-volume" * "\t" * string(S_g/node[node_index].S) * "\n"		# volume of gas relative to 1 L of water
		input_string *= "-temperature" * "\t" * string(T-273) * "\n"
		# specify initial (start-of-time-step) partial pressures of each gas component
		for (i, gs) in enumerate(gas)
			input_string *= gs.name * "(g)" * "\t" * string(node[node_index].Xg[i]*P_atm) * "\n"
		end		
	end
	
	# write to selected output
	input_string *= "SELECTED_OUTPUT 1" * "\n"
    input_string *= "-reset" * "\t" * "false" * "\n"
	input_string *= "-pH" * "\t" * "true" * "\n"
	input_string *= "-pe" * "\t" * "true" * "\n"	
	input_string *= "-totals"
	for cmp in comp
		if (cmp.master != "pH") & (cmp.master != "pe")
			input_string *= " " * cmp.master * " "
		end
	end
	input_string *= "\n"
	input_string *= "-equilibrium_phases "	
	for mins in minerals.names
		input_string *= " " * mins * " "	
	end
	input_string *= "\n"	
	input_string *= "-gases"	
	for gs in gas
		input_string *= " " * gs.name * "(g)" * " "	
	end	
	input_string *= "\n"
	
	input_string *= "END" * "\n"	
	return input_string
	
end
