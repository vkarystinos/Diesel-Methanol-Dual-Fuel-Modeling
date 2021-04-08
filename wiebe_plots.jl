### A Pluto.jl notebook ###
# v0.12.21

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

# ╔═╡ 776e40f2-3ec5-11eb-34bd-9f43fe7f1703
begin
	using PlutoUI
	using Plots
	using DelimitedFiles
	using LaTeXStrings
	using NumericalIntegration
	pyplot()
end

# ╔═╡ 7a9502b0-3eb2-11eb-3793-d71407022bb5
md"# Wiebe Function Plots

_Created by Vasileios Karystinos_

_Version p.0.1_

This notebook plots the single and double wiebe function, using the predifined coeeficients from the user. Also additional data from the user can be provided in order to compare the the shape of wiebe function with the actual heat realese fraction.
"

# ╔═╡ 14bd33dc-3eb4-11eb-2ccc-d1b462dfe746
md" ## _Wiebe Function_

Wiebe function purpose is to model the mass fraction burned as a function of engine position [1]. The Wiebe function is expressed as:

$$x_b= \frac{m_b}{m_{tot}} = 1 - exp \left[ 
-\alpha\left( \frac{\theta-\theta_0} {\Delta\theta} \right) ^{m+1} 
\right]$$

where:

-  $x_b$ is the mass fraction burned
-  $\theta$ is the crank angle
-  $\theta_0$ indicates the start of combustion
-  $\Delta\theta$  is the total combustion duration.
-  $m$ s the form factor and it determines the shape of the combustion process curve 
-  $\alpha$ is the efficiency parameter and it controls the duration of the combustion process.

The Wiebe function curve has a characteristic S-shaped curve and is commonly used to
characterize the combustion process. The mass fraction burned profile grows from zero, where zero mass fraction burn indicates the start of combustion and then tends exponentially to one indicating the end of combustion. The difference between those two ends is known as the duration of combustion.

In many cases, the above model is not sufficient to describe the combustion process. So the double-Wiebe function is introduced:

$$x_b  = \lambda \left( 1 - exp \left[ 
-\alpha_1\left( \frac{\theta-\theta_{0,1}} {\Delta\theta_1} \right) ^{m_1+1} 
\right] \right) + (1-\lambda) \left( 1 - exp \left[ 
-\alpha_2\left( \frac{\theta-\theta_{0,2}} {\Delta\theta_2} \right) ^{m_2+1} 
\right] \right)$$

Where $\lambda$ is the Amplitude Correction Factor of the Wiebe function that represents the contribution of each combustion stage to the heat release. Here there are seven optimazation factors $\alpha_1, \alpha_2, m_1, m_2, \lambda, \Delta\theta_1, \Delta\theta_2$. 

The derivative of the above is the burning rate:

$$\frac{dx_b}{d\theta} =\left( \frac{a(m+1)}{\Delta \theta}\right) \left( \frac{\theta-\theta_0} {\Delta\theta} \right) ^{m} exp \left[ 
-\alpha\left( \frac{\theta-\theta_0} {\Delta\theta} \right) ^{m+1} 
\right]$$

and 

$$\begin{align}
\frac{dx_b}{d\theta} =  &&\lambda \left( \left( \frac{a_1(m_1+1)}{\Delta \theta_1}\right) \left( \frac{\theta-\theta_{0,1}} {\Delta\theta_1} \right) ^{m_1} exp \left[ 
-\alpha_1\left( \frac{\theta-\theta_{0,1}} {\Delta\theta_1} \right) ^{m_1+1} 
\right] \right) + \\
 &&  (1-\lambda)\left( \left( \frac{a_2(m_2+1)}{\Delta \theta_2}\right) \left( \frac{\theta-\theta_{0,2}} {\Delta\theta_2} \right) ^{m_2} exp \left[ 
-\alpha_2\left( \frac{\theta-\theta_{0,2}} {\Delta\theta_2} \right) ^{m_2+1} 
\right] \right)
\end{align}$$

"

# ╔═╡ 924ac326-3fdf-11eb-2d11-89173d5e4a47
md" ## Plotting Wiebe Functions

In this section wiebe functions are plotted, all you need is to import data and choose the appropriate wiebe coeeficients"


# ╔═╡ 75d81fb0-3fce-11eb-1147-915f695d8ea7
md" _Loading Libraries_"

# ╔═╡ 6e76cfd2-3fe7-11eb-2f98-f979044ec90b
md"### _Data Import_

Please enter file name, containing your data"

# ╔═╡ 9163ea42-403a-11eb-078f-63c5bf4f1755
@bind filename TextField(default = "HHR_1.3_Mpa_1900_rpm.txt")

# ╔═╡ d2aa353c-4066-11eb-393c-41f6aba4df2e
md"##### _Insert data type_"

# ╔═╡ 607cca16-4060-11eb-228d-85b5f431da32
@bind data_type Select(["Heat Realese Rate", "Burnt Fraction"])

# ╔═╡ 162362f2-4067-11eb-0dfa-ebb81c4bee02
md" ##### _Insert data units_"

# ╔═╡ aaa41fcc-4060-11eb-0a1e-5d8fab35f3b4
begin
	if data_type == "Heat Realese Rate"
		@bind data_units Select(["J/CA", "kJ/CA", "BMF/CA"])
	elseif data_type == "Burnt Fraction"
		@bind data_units Select(["J","kJ","BMF"])
	end
end

# ╔═╡ 47851714-4067-11eb-0962-99d41d49aaa6
md"##### _Insert required engine data units_"

# ╔═╡ 9b2a3cd2-4076-11eb-37a0-e971da6eb6dd
md"_**Number of Cylinders**_"

# ╔═╡ a87c2ea4-4076-11eb-2200-351da994414e
@bind z NumberField(0:1:100, default=6)

# ╔═╡ 04e7124e-4077-11eb-211f-fb9a3accca04
md"_**4-Stroke or 2 Stroke**_"

# ╔═╡ 05d08634-4077-11eb-21ac-491bc6304f9b
@bind K  Slider(2:2:4, default=4)

# ╔═╡ 7084e7b2-4077-11eb-2f44-e79c0595d30b
Print("$K Stroke")

# ╔═╡ eb181e86-406b-11eb-2b26-c1c70b39ca73
md"_**Engine Speed [rpm]**_"

# ╔═╡ b556300c-406c-11eb-3a93-2ff417ac389a
@bind speed NumberField(0:10000, default=1000)

# ╔═╡ ec85d8f8-406b-11eb-2979-3fd9036ef70b
md"_**Fuel Consumption [kg/h]**_"

# ╔═╡ c7d9e034-406c-11eb-2c02-e3fc022c4627
@bind mb NumberField(0:500, default=40)

# ╔═╡ 14c1eb68-406c-11eb-3d08-e1241c7e75d2
md"_**Lower Heating Value [MJ/kg]**_"

# ╔═╡ f2140ea6-4076-11eb-0e7a-7f22c36c1423
@bind LHV NumberField(0:100, default=42.8)

# ╔═╡ f7310372-4075-11eb-23c3-211900d790fb
begin
	# Data Load
	comb_data = readdlm(filename, '\t', Float64, '\n')
	
	theta_dat = comb_data[:,1];
	meas_dat = comb_data[:,2];
	
	# Data Processing
	v = speed/(30*K);
	B = mb/(z*v);
	Qcyl = B*LHV*10^6/3600;
	
	for i = 1:1:length(theta_dat)
		if meas_dat[i] <0
			meas_dat[i]=0
		end
	end
	
	
	if data_type == "Heat Realese Rate"
		if data_units == "J/CA"
			dxb = meas_dat/Qcyl;
			xb = cumul_integrate(theta_dat,dxb);
		elseif data_units == "kJ/CA"
			dxb = meas_dat/Qcyl*10^3;
			xb = cumul_integrate(theta_dat,dxb);
		else
			dxb = meas_dat;
			xb = cumul_integrate(theta_dat,dxb);
		end
	elseif data_type == "Burnt Fraction"
		if data_units == "J"
			xb = meas_dat/Qcyl;
		elseif data_units == "kJ"
			xb = meas_dat/Qcyl*10^3;
		else
			xb = meas_dat;
		end
	end
end

# ╔═╡ 69b9a816-40a5-11eb-0e25-73d7ac8480ed
datapl1 = plot(theta_dat, meas_dat,xlabel="CAD",ylabel="$data_type  ($data_units)",title = "Data Inserted")

# ╔═╡ 6b02c32e-40a5-11eb-1b31-277a155a0503
	datapl2 = plot(theta_dat, xb,xlabel="CAD",ylabel="$data_type (BMF)",title = "Burnt Mass Fraction")

# ╔═╡ 6bb3753e-40a5-11eb-3d07-8b5fab0298da
if data_type == "Heat Realese Rate"
	datapl3 = plot(theta_dat, dxb,xlabel="CAD",ylabel="$data_type $data_units",title = "Burnt Mass Fraction Rate")
end

# ╔═╡ 5cb52d20-3fe2-11eb-2057-c3ef7cd9c88a
md"_Define Wiebe Functions in julia_"

# ╔═╡ 4d763146-40a1-11eb-3043-e959256b6372
md"##### _Insert Wiebe Function Parameters_"

# ╔═╡ 8688fe14-40a1-11eb-18f6-e367bdc6589f
md"_**Calculate Single Wiebe Function**_"

# ╔═╡ 6efcb54c-40a1-11eb-0c42-7b1b2b527aed
@bind check_s CheckBox(default=true)

# ╔═╡ 89559620-40a1-11eb-2f8a-8d32af0d9f3b
md"_**Calculate Double Wiebe Function**_"

# ╔═╡ 6e87dab0-40a1-11eb-3a35-f7c5317cc4f1
@bind check_d CheckBox(default=true)

# ╔═╡ 10568dba-40a6-11eb-1244-f584a08b24de
md"###### _Single Wiebe Parameters_"

# ╔═╡ c05ca01e-40a1-11eb-3681-b792e7a535b2
md"_**Form Factor $m$**_"

# ╔═╡ 402c08a8-40a6-11eb-020c-3d6beb807d5c
@bind m0 NumberField(0:100, default=2)

# ╔═╡ 44e5912a-40a1-11eb-0213-257019cf1286
md"_**Efficiency Parameter $\alpha$**_"

# ╔═╡ 763d643c-40a6-11eb-31ef-1b2da2561813
@bind a0 NumberField(0:100, default=6)

# ╔═╡ 85827694-40a6-11eb-32da-d78c4c4b107c
md"_**Start of Combustion**_ $\theta_0$"

# ╔═╡ 8550cfe0-40a6-11eb-0f65-37ec392219e5
@bind theta00 NumberField(0:100, default=0)

# ╔═╡ db82fdca-40a6-11eb-2cc1-6bda16a70ec8
md"_**Combustion Duration**_ $\Delta\theta$"

# ╔═╡ dd04973a-40a6-11eb-1428-4b53d3a056f2
@bind Dtheta0 NumberField(0:100, default=50)

# ╔═╡ 22ad599a-3fe2-11eb-0471-153426bb0ea4
#Single Wiebe Function
wiebe_1(theta) = (1 - exp(-a0*(((theta-theta00)/Dtheta0)*(theta-theta00>=0))^(m0+1)))

# ╔═╡ c967589a-41d0-11eb-1b89-87db1e9e2539
#Single Derivative Wiebe Function
dwiebe_1(theta) = ((a0*(m0+1))/(Dtheta0))*(((theta-theta00)/(Dtheta0)*(theta-theta00>=0))^m0)*exp(-a0*(((theta-theta00)/Dtheta0)*(theta-theta00>=0))^(m0+1))

# ╔═╡ c397a120-6090-11eb-0f3e-0bd491cc7935
theta = collect(-20:0.01:50)

# ╔═╡ 8e2739dc-6091-11eb-36fc-616033ef5a3b
if check_s
	xbw1 = wiebe_1.(theta)
	plot(theta,xbw1,xlabel = "CAD",ylabel = "BMF", title = "Single Wiebe Function vs Data", labels = "Double Wiebe Function",legend = :bottomright)
	plot!(theta_dat, xb, labels= "Data")
end

# ╔═╡ 00048dd8-41eb-11eb-1e63-258821c45848
if check_s
	dxbw1 = dwiebe_1.(theta)
	plot(theta,dxbw1,xlabel = "CAD",ylabel = "Rate of BMF", title = "Single Wiebe Function vs Data", labels = "Single der. Wiebe Function",legend = :bottomright)
	plot!(theta_dat,dxb, labels = "Data" )
end

# ╔═╡ 7dd97b6c-40a7-11eb-14f5-2f96654f5bc7
md"###### _Double Wiebe Parameters_"

# ╔═╡ 7a690ea0-415b-11eb-2b47-1b328d905b8f
md"
_**Form Factor**_ $m_1$"

# ╔═╡ 8cfb60ce-415b-11eb-2552-59db5843dc52
@bind m1 NumberField(0:100, default=2)

# ╔═╡ 8df487be-415b-11eb-3813-59dae073f979
md"
_**Form Factor**_ $m_2$"

# ╔═╡ 8dc61528-415b-11eb-10fc-797a9f39c488
@bind m2 NumberField(0:100, default=2)

# ╔═╡ 9c97d76a-415b-11eb-27b0-b3504dfaa4c5
md"
_**Efficiency Parameter**_ $\alpha_1$"

# ╔═╡ 58c6686a-417b-11eb-0e77-4f42f276c054
@bind a1 NumberField(0:100, default=6)

# ╔═╡ 5334812a-417b-11eb-1421-3b83576ebbfb
md"
_**Efficiency Parameter**_ $\alpha_2$"

# ╔═╡ 59a14926-417b-11eb-2a6a-99891f24e4a2
@bind a2 NumberField(0:100, default=6)

# ╔═╡ ca826dd2-417b-11eb-2e03-eb701100a7f7
md"_**Start of 1st mode of Combustion**_ $\theta_1$"

# ╔═╡ e4048506-417b-11eb-3122-cbbe7c539ae2
@bind theta01 NumberField(0:100, default=0)

# ╔═╡ dec45c72-417b-11eb-14b1-0d8c920dd9a2
md"_**Start of 2nd mode of Combustion**_ $\theta_2$"

# ╔═╡ e496941c-417b-11eb-3dff-05a10b25267f
@bind theta02 NumberField(0:100, default=0)

# ╔═╡ 045d5dbe-417c-11eb-2bbe-89357a75dd20
md"_**1st mode Combustion Duration**_ $\Delta\theta_1$"

# ╔═╡ 23d02d0e-417c-11eb-3cfd-0147c076e69e
@bind Dtheta1 NumberField(0:100, default=50)

# ╔═╡ 12ae08b4-417c-11eb-09c4-017f25b6b632
md"_**2nd mode Combustion Duration**_ $\Delta\theta_2$"

# ╔═╡ 2455f946-417c-11eb-2b39-fb33d4f23513
@bind Dtheta2 NumberField(0:100, default=50)

# ╔═╡ 42b11a92-417c-11eb-0599-472931d4d466
md"_**Amplitude Correction Factor**_ $\lambda$"

# ╔═╡ b3d49920-41d0-11eb-2bef-8b06db8364b9
@bind lambda NumberField(0:1, default=0.3)

# ╔═╡ 2a3f6e50-3fe2-11eb-26b2-f1273fc7c302
#Double Wiebe Function
wiebe_2(theta) = lambda*(1 - exp(-a1*(((theta-theta01)/Dtheta1)*(theta-theta01>=0))^(m1+1)))+(1-lambda)*(1 - exp(-a2*(((theta-theta02)/Dtheta2)*(theta-theta02>=0))^(m2+1)))

# ╔═╡ ca0531a0-41d0-11eb-2325-dfe886f7c363
#Single Derivative Wiebe Function
dwiebe_2(theta) = lambda*(((a1*(m1+1))/(Dtheta1))*(((theta-theta01)/(Dtheta1)*(theta-theta01>=0))^m1)*exp(-a1*(((theta-theta01)/Dtheta1)*(theta-theta01>=0))^(m1+1)))+(1-lambda)*(((a2*(m2+1))/(Dtheta2))*(((theta-theta02)/(Dtheta2)*(theta-theta02>=0))^m2)*exp(-a2*(((theta-theta02)/Dtheta2)*(theta-theta02>=0))^(m2+1)))

# ╔═╡ e6662f96-41eb-11eb-3a01-0f5a4c7183db
if check_d
	xbw2 = wiebe_2.(theta)
	plot(theta,xbw2,xlabel = "CAD",ylabel = "BMF", title = "Single Wiebe Function vs Data", labels = "Double Wiebe Function",legend = :bottomright)
	plot!(theta_dat, xb, labels= "Data")
end

# ╔═╡ dfa376d6-41ec-11eb-3d86-fdb48f8b090d
if check_d
	dxbw2 = dwiebe_2.(theta)
	plot(theta,dxbw2,xlabel = "CAD",ylabel = "Rate of BMF", title = "Single Wiebe Function vs Data", labels = "Single der. Wiebe Function",legend = :bottomright)
	plot!(theta_dat,dxb, labels = "Data" )
end

# ╔═╡ 605896b4-60cc-11eb-2f5f-0fd8a39a09fa
if check_d
	p1 = plot(theta,xbw2,xlabel = "CAD",ylabel = "BMF", title = "Burnt Mass Fraction Simulation", labels = "Wiebe Function",color = "red", legend = :topleft)
	plot!(theta_dat, xb, labels= "Data",color = "blue")
	

	p2 = plot(theta,dxbw2,xlabel = "CAD",ylabel = "HRR", title = "Heat Realese Rate Simulation", labels = "Wiebe Function",color = "red", legend = :topleft)
	plot!(theta_dat,dxb, labels = "Data",color = "blue")
	
	fn = plot(p1,p2, layout = (2,1) )
 	
end

# ╔═╡ 139297ac-614f-11eb-1faa-7dd6b0361014
begin	
# Data Load
	comb_data1 = readdlm("Pressure.txt", '\t', Float64, '\n')
	
	theta1_dat = comb_data1[:,1].-480;
	meas1_dat = comb_data1[:,2];
	
	meas1_dat = meas1_dat.-0.5*theta1_dat.+10
	
	comb_data2 = readdlm("press2.txt", '\t', Float64, '\n')
	theta2_dat = comb_data2[:,1];
	meas2_dat = comb_data2[:,2];

end

# ╔═╡ 137e8dca-614f-11eb-1944-8980ce03daf6
begin	
# Data Load
	 plot(theta1_dat[1880:2040], meas1_dat[1880:2040], xlabel = "CAD",ylabel = "BAR", title = "In cylinder pressure Simulation", labels = "Simulation Result",color = "red", legend = :bottomleft)
	dn = plot!(theta2_dat,meas2_dat, labels= "Data",color = "blue")

end

# ╔═╡ dfc87f2a-6150-11eb-3c0a-39fbece38317
theta1_dat[1880:2040]

# ╔═╡ 7756306e-60e1-11eb-18ad-7fdb8c071680
gui(dn)

# ╔═╡ 9ad6ea6a-3fd0-11eb-34f3-851bd36d9c8f
md" ### Bibliography
[1] Heywood, J. B., _Internal Combustion Engine Fundamentals._ New York: McGraw-Hill, 1988.

[2] Wiebe, Ivan Ivanovitch. _Progress in engine cycle analysis: Combustion rate and cycle processes_. 1962
"

# ╔═╡ 66391718-8a50-11eb-3417-d58afee9bf00


# ╔═╡ Cell order:
# ╟─7a9502b0-3eb2-11eb-3793-d71407022bb5
# ╟─14bd33dc-3eb4-11eb-2ccc-d1b462dfe746
# ╟─924ac326-3fdf-11eb-2d11-89173d5e4a47
# ╟─75d81fb0-3fce-11eb-1147-915f695d8ea7
# ╠═776e40f2-3ec5-11eb-34bd-9f43fe7f1703
# ╟─6e76cfd2-3fe7-11eb-2f98-f979044ec90b
# ╟─9163ea42-403a-11eb-078f-63c5bf4f1755
# ╟─d2aa353c-4066-11eb-393c-41f6aba4df2e
# ╟─607cca16-4060-11eb-228d-85b5f431da32
# ╟─162362f2-4067-11eb-0dfa-ebb81c4bee02
# ╟─aaa41fcc-4060-11eb-0a1e-5d8fab35f3b4
# ╟─47851714-4067-11eb-0962-99d41d49aaa6
# ╟─9b2a3cd2-4076-11eb-37a0-e971da6eb6dd
# ╟─a87c2ea4-4076-11eb-2200-351da994414e
# ╟─04e7124e-4077-11eb-211f-fb9a3accca04
# ╟─05d08634-4077-11eb-21ac-491bc6304f9b
# ╟─7084e7b2-4077-11eb-2f44-e79c0595d30b
# ╠═eb181e86-406b-11eb-2b26-c1c70b39ca73
# ╠═b556300c-406c-11eb-3a93-2ff417ac389a
# ╟─ec85d8f8-406b-11eb-2979-3fd9036ef70b
# ╠═c7d9e034-406c-11eb-2c02-e3fc022c4627
# ╟─14c1eb68-406c-11eb-3d08-e1241c7e75d2
# ╠═f2140ea6-4076-11eb-0e7a-7f22c36c1423
# ╠═f7310372-4075-11eb-23c3-211900d790fb
# ╟─69b9a816-40a5-11eb-0e25-73d7ac8480ed
# ╟─6b02c32e-40a5-11eb-1b31-277a155a0503
# ╟─6bb3753e-40a5-11eb-3d07-8b5fab0298da
# ╟─5cb52d20-3fe2-11eb-2057-c3ef7cd9c88a
# ╟─4d763146-40a1-11eb-3043-e959256b6372
# ╠═22ad599a-3fe2-11eb-0471-153426bb0ea4
# ╟─2a3f6e50-3fe2-11eb-26b2-f1273fc7c302
# ╠═c967589a-41d0-11eb-1b89-87db1e9e2539
# ╟─ca0531a0-41d0-11eb-2325-dfe886f7c363
# ╟─8688fe14-40a1-11eb-18f6-e367bdc6589f
# ╟─6efcb54c-40a1-11eb-0c42-7b1b2b527aed
# ╟─89559620-40a1-11eb-2f8a-8d32af0d9f3b
# ╟─6e87dab0-40a1-11eb-3a35-f7c5317cc4f1
# ╟─10568dba-40a6-11eb-1244-f584a08b24de
# ╟─c05ca01e-40a1-11eb-3681-b792e7a535b2
# ╟─402c08a8-40a6-11eb-020c-3d6beb807d5c
# ╟─44e5912a-40a1-11eb-0213-257019cf1286
# ╟─763d643c-40a6-11eb-31ef-1b2da2561813
# ╟─85827694-40a6-11eb-32da-d78c4c4b107c
# ╟─8550cfe0-40a6-11eb-0f65-37ec392219e5
# ╟─db82fdca-40a6-11eb-2cc1-6bda16a70ec8
# ╟─dd04973a-40a6-11eb-1428-4b53d3a056f2
# ╠═c397a120-6090-11eb-0f3e-0bd491cc7935
# ╠═8e2739dc-6091-11eb-36fc-616033ef5a3b
# ╠═00048dd8-41eb-11eb-1e63-258821c45848
# ╟─7dd97b6c-40a7-11eb-14f5-2f96654f5bc7
# ╟─7a690ea0-415b-11eb-2b47-1b328d905b8f
# ╟─8cfb60ce-415b-11eb-2552-59db5843dc52
# ╟─8df487be-415b-11eb-3813-59dae073f979
# ╟─8dc61528-415b-11eb-10fc-797a9f39c488
# ╟─9c97d76a-415b-11eb-27b0-b3504dfaa4c5
# ╟─58c6686a-417b-11eb-0e77-4f42f276c054
# ╟─5334812a-417b-11eb-1421-3b83576ebbfb
# ╟─59a14926-417b-11eb-2a6a-99891f24e4a2
# ╟─ca826dd2-417b-11eb-2e03-eb701100a7f7
# ╟─e4048506-417b-11eb-3122-cbbe7c539ae2
# ╟─dec45c72-417b-11eb-14b1-0d8c920dd9a2
# ╟─e496941c-417b-11eb-3dff-05a10b25267f
# ╟─045d5dbe-417c-11eb-2bbe-89357a75dd20
# ╟─23d02d0e-417c-11eb-3cfd-0147c076e69e
# ╟─12ae08b4-417c-11eb-09c4-017f25b6b632
# ╟─2455f946-417c-11eb-2b39-fb33d4f23513
# ╟─42b11a92-417c-11eb-0599-472931d4d466
# ╟─b3d49920-41d0-11eb-2bef-8b06db8364b9
# ╠═e6662f96-41eb-11eb-3a01-0f5a4c7183db
# ╠═dfa376d6-41ec-11eb-3d86-fdb48f8b090d
# ╠═605896b4-60cc-11eb-2f5f-0fd8a39a09fa
# ╠═139297ac-614f-11eb-1faa-7dd6b0361014
# ╠═137e8dca-614f-11eb-1944-8980ce03daf6
# ╠═dfc87f2a-6150-11eb-3c0a-39fbece38317
# ╠═7756306e-60e1-11eb-18ad-7fdb8c071680
# ╟─9ad6ea6a-3fd0-11eb-34f3-851bd36d9c8f
# ╠═66391718-8a50-11eb-3417-d58afee9bf00
