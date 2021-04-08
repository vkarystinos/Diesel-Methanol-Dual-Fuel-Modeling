### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 04161a70-72a6-11eb-3255-fb86cd33f09c
begin
	# Constants
	using DelimitedFiles
	using Plots
	using Dierckx
	using PlutoUI
	using NumericalIntegration
	plotly()
end

# ╔═╡ 355886a6-8c94-11eb-37ee-df3c2bf3d161
md"# In cylinder Pressure and Heat Release Rate  Calculation

_Created by V. Karystinos_

_Version f.1.0_

This notebook performs the appropriate calculations in order to export the Net Heat Release Rate $\frac{dQ_n}{dθ}$ and the Heat Loss Rate $\frac{dQ_w}{dθ}$. Also, in-cylinder temperature is calculated and the parameters of _Woschni_ model are calculated. Finally, data are smoothen and exported for further processing
"

# ╔═╡ 6ad766e4-8ccb-11eb-0679-f5ca883233fa
begin
	FilesHRR = readdir("Data/Raw_Data/HRR_Data")
	FilesPress = readdir("Data/Raw_Data/Pressure_Data")
end;

# ╔═╡ 424e90aa-8ccb-11eb-0190-21a63a36b544
md"### _Data Import_

Choose the files containing the combustion data and upload the files. The file types should be in .csv or.txt form and the stracture should be two columns, the fisrt reffering to CAD an the second to the appropriate combustion. Also you need to choose the appropriate units of your data

**_Combustion Data_**\
_File_: $(@bind filename_hr Select(FilesHRR))\
_Units_: $(@bind data_units_hr Select([\"J/CA\", \"kJ/CA\", \"BMF/CAD\"]))


**_In Cylinder Pressure Dat_**\
_File_: $(@bind filename_pr Select(FilesPress))\
_Units_: $(@bind data_units_pr Select([\"MPa\", \"kPa\",\"bar\"]))

"

# ╔═╡ 770f0f50-8ccb-11eb-3d1b-75efc2785534
md"##### Engine and Condition Data

Now, the following table regarding engine characteristics and operation conditions should be filled accordingly

| Property | Value | 
| ------ | ------ |
| _2 or 4 stroke_ | $(@bind K NumberField(2:2:4, default=4)) |
| _Bore [mm]_ | $(@bind B NumberField(10:1:1000, default=126)) |
| _Stroke [mm]_ | $(@bind S NumberField(10:1:1000, default=130)) |
| _Connecting rod length [mm]_ | $(@bind l NumberField(10:1:1000, default=219)) |
| _Number of Cylinders_ | $(@bind z NumberField(0:1:100, default=6)) |
| _Displacment [L]_ | $(@bind V NumberField(0:0.001:20, default=9.726)) |
| _Engine Speed [rpm]_ | $(@bind speed NumberField(0:1:9999, default=1900)) |
| _Fuel 1 mass [kg/h]_ | $(@bind mf1 NumberField(0:0.01:100, default=40)) |
| _Fuel 2 mass [kg/h]_ | $(@bind mf2 NumberField(0:0.01:100, default= 0)) |
| _Overall air fuel ratio_ | $(@bind AF_ratio NumberField(0:0.01:100, default=30)) |
"

# ╔═╡ 7d9cfe12-89c7-11eb-1ba9-3ffa0be6477d
begin
	# Data Load
	comb_data = readdlm(join(["Data/Raw_Data/HRR_Data/", filename_hr],""), ',','\n')
	press_data = readdlm(join(["Data/Raw_Data/Pressure_Data/", filename_pr],""), ',', '\n')
	
	if data_units_pr == "MPa"
		cupr = 1e6
	elseif data_units_pr == "kPa"
		cupr = 1e3
	else
		cupr = 1e5
	end
	
	
	if data_units_pr == "J/CA"
		cuhrr = 1
	elseif data_units_pr == "kJ/CA"
		cuhrr = 1e3
	else
		cuhrr = 1
	end
	
	theta_dat_pr = press_data[:,1]
	meas_dat_pr = cupr*press_data[:,2].+10^5
	
	theta_dat_hrr = comb_data[:,1];
	meas_dat_hrr = comb_data[:,2];
	
	
	#Plotting
	plt1 = plot(theta_dat_pr, meas_dat_pr,xlabel="CAD",ylabel="Pressure [Pa]", linewidth = 2, color = "red")
	plt2 = plot(theta_dat_hrr, meas_dat_hrr,xlabel="CAD",ylabel="Heat Release Rate [J/CAD]", linewidth = 2, color = "green")
	plot(plt1, plt2, layout = (2,1), legend = false)

end

# ╔═╡ d53d40be-8cd3-11eb-2d79-d3d4635dc6f0
spl_theta = collect(ceil(theta_dat_pr[1]):0.001:floor(theta_dat_pr[length(theta_dat_pr)]));

# ╔═╡ 07ef7a5e-7429-4e42-9bd7-628f21c98bd9
md"### _Smoothening Pressure_

In order to calculate the derivatives, pressure data need to be smoothened. In order to this to be done, data are replaced by a 1-D 3rd order spline, with the following filtering parameters:

_**In cylinder Pressure Filt. Param.**_\
$(@bind flt_pr NumberField(0:1e9:1e11, default=1e10))\

"

# ╔═╡ e23dfc7c-89c7-11eb-0aa1-31d6a8017654
begin
	spl_pr = Spline1D(theta_dat_pr, meas_dat_pr; w=ones(length(theta_dat_pr)), k=4, bc="nearest", s=flt_pr);
	
	#Pressure Data
	p_sm = spl_pr(spl_theta);
	dp_sm = derivative(spl_pr, spl_theta);

end;

# ╔═╡ 0a353b5a-89c8-11eb-25f8-a194045f8dea
begin
	scatter(theta_dat_pr, meas_dat_pr, xlabel="CAD", ylabel="Pressure [Pa]", title = "Pressure Smoothening", label = "Raw Data", color = "red")
	plot!(spl_theta, spl_pr(spl_theta), label = "Smooth Spline", linewidth = 2, color = "blue")
end

# ╔═╡ 9ca897b6-8cd4-11eb-2e56-8b698ce7d500
md"### _Net Heat Release Rate Calulation_

In order to calculate the Heat Release Rate, the 1st thermodynamic law for the cylinder is employed and formulated as follows

$$\frac{dQ_n}{dφ} = \frac{γ}{γ-1}\cdot p \cdot \frac{dV}{dφ}+ \frac{1}{γ-1}\cdot V \cdot \frac{dp}{dφ}$$

where $γ$, $p$ and $V$ is the ratio of specific heat, cylinder pressure and cylinder volume, respectively. Below you must choose a constant value fo the specific heat ratio.


_**Ratio of Specific Heat. γ**_\
$(@bind γ NumberField(1:0.01:2, default=1.35))

The instantaneous cylinder volume and its derivative are calculated via the following formulas:

$$V(φ) = V_c +A_{pc}\cdot r\cdot [(1-cosφ)+λ^{-1}-\sqrt{λ^{-2}-sin^2φ}]$$

$$\frac{dV(φ)}{dφ} = A_{pc}\cdot r\cdot \left[sinφ(1+\frac {λcosφ}{\sqrt {1-λ^2sin^2φ}})\right]$$

"

# ╔═╡ b3db40c4-8c83-11eb-0e9a-8febe74e7d51
function V_pist(θ)

	Vc = 0.0006080625/6 #[m^3]
	
	S = 130*0.001	    #[m]
	B = 126*0.001		#[m]
	l = 219*0.001		#[m] 
	
		
	Apc = (pi/4)*B^2
	r = S/2
	λ = r/l
	
	
	φ = deg2rad(θ)

	
	V = Vc + Apc*r*((1-cos(φ))+λ^-1-sqrt(λ^-2-sin(φ)^2))
	
	return V
end

# ╔═╡ 502cd07a-8c87-11eb-0a02-9932ca070bee
function dV_pist(θ)

	Vc = 0.0006080625/6 

	S = 130*0.001	    #[m]
	B = 126*0.001		#[m]
	l = 219*0.001		#[m] 
	
		
	Apc = (pi/4)*B^2
	r = S/2
	λ = r/l
	
	
	φ = deg2rad(θ)

	
	dV = (π/180)*Apc*r*sin(φ)*(1+(λ*cos(φ)/(sqrt(1-(λ^2)*sin(φ)^2))))
	
	return dV
end

# ╔═╡ ab36606a-8c19-11eb-3acc-317c073b1ba6
begin
	dQn = zeros(length(p_sm))

	
	for i = 1:length(p_sm)
		θ = spl_theta[i]
		
		dp = dp_sm[i]
		dV = dV_pist(θ)
		
		dQn[i] = (γ/(γ-1))*p_sm[i]*dV+(1/(γ-1))*V_pist(θ)*dp
	end
		
	
end

# ╔═╡ 04235cf6-8c1b-11eb-1ba6-cb630268097a
begin
	plot(spl_theta, dQn,label="Net Heat Release Rate (calc.)", linewidth = 2)
	plot!(theta_dat_hrr, meas_dat_hrr,label="Gross Heat Release Rate (given)", linewidth = 2, legend=:bottomright)
	plot!(xlabel="CAD",ylabel="Heat Release Rate [J/CAD]")
end

# ╔═╡ 094a73b3-5f92-4a61-b8cc-1fb5201a38d8
begin
	spl_hrr = Spline1D(spl_theta, dQn; w=ones(length(spl_theta)), k=5,bc="nearest", s=10000)
	hrr =  spl_hrr(spl_theta)
end;

# ╔═╡ 9b519213-7fbe-4a43-831b-4836a8e86427
begin
	plot(spl_theta, hrr)
	plot!(spl_theta, dQn)
end

# ╔═╡ 6c93f6e0-8eca-11eb-2fc5-fd42130d5ac9
md"### _In cylinder Temperature_

In order to calulate in cylinder temperature, the gas state equation is employed:

$$T = \frac{pV}{mR}$$

where R is the gas constant.

_**Gas Constant R**_\
$(@bind R NumberField(150:0.01:350, default=287))


The mass m is the total mass inside the cylider and is calculated via air/fuel mass ratio and fuel mass.
"


# ╔═╡ dbc9de38-8ecc-11eb-20a1-632d5440b0d0
begin
	#Air mas calculation
	mair = AF_ratio*(mf1+mf2)
	mair_cyl = mair/60/(speed/(K/2))/z
	
	mf1_cyl =  mf1/60/(speed/(K/2))/z
	mf2_cyl =  mf2/60/(speed/(K/2))/z
	
	m_cyl = (mair_cyl+mf1_cyl+mf2_cyl)
end;

# ╔═╡ 292eaa0a-0a9b-44ed-96dc-41c5ca78dc52
begin
	LHV1 = 42.5
	LHV2 = 19.7
	
	LHV = (mf1*LHV1+mf2*LHV2)/(mf1+mf2)
	v = speed/(30*K);
	Bm = (mf1+mf2)/(z*v);
	Qcyl = Bm*LHV*10^6/3600
	
	for i = 1:1:length(spl_theta)
		if hrr[i] <0
			hrr[i]=0
		end
	end
	
	if data_units_hr == "J/CA"
		dxb = hrr/Qcyl;
		xb = cumul_integrate(spl_theta,dxb);
	elseif data_units_hr == "kJ/CA"
		dxb = hrr/Qcyl*10^3;
		xb = cumul_integrate(spl_theta,dxb);
	end
	
end
	

# ╔═╡ fe946179-86d3-4311-852f-aacf78858c0a
plot(spl_theta, xb)

# ╔═╡ 73d5428b-9ef7-4a8d-ba58-ff33e23e0d00
plot(spl_theta, dxb)

# ╔═╡ 5479a62e-3fcb-4328-ad67-576c5f41ee9c
begin
	
	open("BMF_1900rpm_154_meth.txt", "w") do io
		writedlm(io, [spl_theta xb])
	end
	
	open("BMF_rate_1900rpm_154_meth.txt", "w") do io
		writedlm(io, [spl_theta dxb])
    end
	
	open("Press_1900rpm_154_meth.txt", "w") do io
		writedlm(io, [spl_theta p_sm])
    end
end

# ╔═╡ a444064f-50d4-4729-8cd2-a3f50990597c
T_cyl = (p_sm.*V_pist.(spl_theta))./(m_cyl*R)

# ╔═╡ 030adb20-8ed5-11eb-2d8c-793eea777e72
begin
	plot(spl_theta,T_cyl,xlabel = "CAD", ylabel = "In cylinder temperature [K]", legend = false, linewidth = 2)		
end

# ╔═╡ 425850e6-8f3e-11eb-1c22-6911eddd09ce
begin
	
	MWCH3OH = 32
	MWC12H23 = 167
	MWCO2 = 44
	MWH2O = 18
	MWO2 = 32
	MWN2 = 28
	MWair = 0.21*MWO2+0.79*MWN2
	
	
	cd = 12*12/MWC12H23
	hd = 23*1/MWC12H23
	
	cm = 1*12/MWCH3OH
	hm = 4*1/MWCH3OH
	om = 16*1/MWCH3OH
	
	AFstd = 71/4*(MWO2+3.76*MWN2)/MWC12H23
	AFstm = 3/2*(MWO2+3.76*MWN2)/MWCH3OH
	
	mCO2d = cd/12*MWCO2*mf1_cyl
	mH2Od = hd*0.5*MWH2O*mf1_cyl
	
	mCO2m = cm/12*MWCO2*mf2_cyl
	mH2Om = hm*0.5*MWH2O*mf2_cyl
	
	mCO2 = mCO2d+mCO2m
	mH2O = mH2Od+mH2Om
	
	
	mO2d = (cd/12+hd/4)*MWO2*mf1_cyl
	mO2m = (cm/12+hm/4-om/32)*MWO2*mf2_cyl
	
	mO2e = 0.232*mair_cyl-mO2d-mO2m
	mN2e = 0.768*mair_cyl
	
	
	mO2i = 0.232*mair_cyl
	mN2i = mN2e
	#mOe = 1(cd/12+hd/4)*(λtot-1)*MWO2*mf1_cyl
	
	#mN2 = 0.79*100/21*(cd/12+hd/4)*(λtot)*MWN2*mf1_cy
	

end

# ╔═╡ 9bb17f40-93b3-11eb-1e7f-17e79373a0ae
begin
	wC12H23 = mf1_cyl/m_cyl
	wCH3OH = mf2_cyl/m_cyl
	wO2i = mO2i/m_cyl
	wN2i = mN2i/m_cyl
	wAir = mair_cyl/m_cyl
	
	MWmi =  ((wC12H23/MWC12H23)+(wCH3OH/MWCH3OH)+(wO2i/MWO2)+(wN2i/MWN2))^-1
	
	xC12H23 =(wC12H23/MWC12H23)*MWmi
	xCH3OH =(wCH3OH/MWCH3OH)*MWmi
	xO2i = (wO2i/MWO2)*MWmi
	xN2i = (wN2i/MWN2)*MWmi
	
	xair = (wAir/MWair)*MWmi
	
	w_O_dual = xair/(4.76*(xair+xCH3OH))
	
	xair_1 = ((mair_cyl+mf2_cyl)/m_cyl)/MWair*MWmi
	
	w_O_diesel = xair_1/(4.76*(xair_1))
	
	r = w_O_dual/w_O_diesel
	
end

# ╔═╡ 4829d38c-93b0-11eb-2738-63b1bc2aa4a5
begin
	RO2 = 259.84
	RN2 = 296.80
	RCO2 = 188.92
	RH2O = 461.52
	RC12H23 = 8.31446261815324*1000/MWC12H23
	RCH3OH = 259.5
	

	R_in = (RO2*mO2i+RN2*mN2i+RC12H23*mf1_cyl+RCH3OH*mf2_cyl)/m_cyl
	R_out = (RO2*mO2e+RN2*mN2e+RCO2*mCO2+RH2O*mH2O)/m_cyl
	
end

# ╔═╡ Cell order:
# ╟─355886a6-8c94-11eb-37ee-df3c2bf3d161
# ╟─04161a70-72a6-11eb-3255-fb86cd33f09c
# ╠═6ad766e4-8ccb-11eb-0679-f5ca883233fa
# ╟─424e90aa-8ccb-11eb-0190-21a63a36b544
# ╟─770f0f50-8ccb-11eb-3d1b-75efc2785534
# ╠═7d9cfe12-89c7-11eb-1ba9-3ffa0be6477d
# ╠═d53d40be-8cd3-11eb-2d79-d3d4635dc6f0
# ╟─07ef7a5e-7429-4e42-9bd7-628f21c98bd9
# ╠═e23dfc7c-89c7-11eb-0aa1-31d6a8017654
# ╠═0a353b5a-89c8-11eb-25f8-a194045f8dea
# ╟─9ca897b6-8cd4-11eb-2e56-8b698ce7d500
# ╟─b3db40c4-8c83-11eb-0e9a-8febe74e7d51
# ╟─502cd07a-8c87-11eb-0a02-9932ca070bee
# ╠═ab36606a-8c19-11eb-3acc-317c073b1ba6
# ╠═04235cf6-8c1b-11eb-1ba6-cb630268097a
# ╠═094a73b3-5f92-4a61-b8cc-1fb5201a38d8
# ╠═9b519213-7fbe-4a43-831b-4836a8e86427
# ╟─6c93f6e0-8eca-11eb-2fc5-fd42130d5ac9
# ╠═dbc9de38-8ecc-11eb-20a1-632d5440b0d0
# ╠═292eaa0a-0a9b-44ed-96dc-41c5ca78dc52
# ╠═fe946179-86d3-4311-852f-aacf78858c0a
# ╠═73d5428b-9ef7-4a8d-ba58-ff33e23e0d00
# ╠═5479a62e-3fcb-4328-ad67-576c5f41ee9c
# ╠═a444064f-50d4-4729-8cd2-a3f50990597c
# ╠═030adb20-8ed5-11eb-2d8c-793eea777e72
# ╠═425850e6-8f3e-11eb-1c22-6911eddd09ce
# ╠═9bb17f40-93b3-11eb-1e7f-17e79373a0ae
# ╠═4829d38c-93b0-11eb-2738-63b1bc2aa4a5
