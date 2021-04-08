### A Pluto.jl notebook ###
# v0.14.1

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

# ╔═╡ 38168156-4200-11eb-1428-df9a090087c2
begin
	using PlutoUI
	using Plots
	using DelimitedFiles
	using LaTeXStrings
	using NumericalIntegration
	using LsqFit
	using ControlSystems
	using Dierckx
	using Roots
	using Peaks
	using LaTeXStrings
	using ControlSystems
	#using CurveFit
	plotly()
end;

# ╔═╡ 7a243fdc-405e-11eb-3e37-cfda2218520d
md"# Wiebe Function Fitting

_Created by V. Karystinos_

_Version f.1.0_

This notebook performs fitting of the coefficeints of the single and double wiebe function, using data provided by the user. The fitting algorythms are described below.
"

# ╔═╡ 79d3def6-405f-11eb-0e3c-43105efaef03
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

# ╔═╡ a7461f68-41ed-11eb-0a48-c1f7c41834a5
md" ## _Fitting Wiebe Functions_

In this section wiebe functions are fitted, all you need is to import data and choose an estimation for wiebe coeeficients"

# ╔═╡ 2d2cf1bc-4200-11eb-3771-61e3cbcedf90
md" _Loading Libraries_"

# ╔═╡ fed92118-729a-11eb-11a5-e5267ce84380
begin
	FilesHRR = readdir("Data/Smoothen_Data/HRR_Data")
	FilesPress = readdir("Data/Smoothen_Data/Pressure_Data")
end;

# ╔═╡ 5763f37c-4200-11eb-2c4c-01586566124c
md"### _Data Import_

##### Measurment Data


Choose the files containing the combustion data and upload the files. The file types should be in .csv or.txt form and the stracture should be two columns, the fisrt reffering to CAD an the second to the appropriate combustion. Also you need to choose the appropriate units of your data


**_Burnt Mass Fraction Data_**\
_File_: $(@bind filename_bmf Select(FilesHRR))\

**_Burnt Mass Fraction Rate Data_**\
_File_: $(@bind filename_hrr Select(FilesHRR))\

**_In Cylinder Pressure Dat_**\
_File_: $(@bind filename_pr Select(FilesPress))\
_Units_: $(@bind data_units_pr Select([\"Mpa\", \"kPa\",\"bar\"]))

"

# ╔═╡ 5985dcd0-4201-11eb-3637-4ba209c60c01
begin
	# Data Load
	bmf_data = readdlm(join(["Data/Smoothen_Data/HRR_Data/", filename_bmf],""), '\t','\n')
	hrr_data = readdlm(join(["Data/Smoothen_Data/HRR_Data/", filename_hrr],""), '\t','\n')
	press_data = readdlm(join(["Data/Smoothen_Data/Pressure_Data/", filename_pr],""), '\t','\n')
	
	θ = bmf_data[:,1];
	
	xb = bmf_data[:,2];
	dxb = hrr_data[:,2];
	press = press_data[:,2];	
end;

# ╔═╡ 42af1e4c-634e-11eb-0a08-bd55c6c2f851
md" #### _Data Plots_"

# ╔═╡ 52a7eacc-4201-11eb-0e52-0d5343c651af
function plot_data()
	plt1 = plot(θ, xb,xlabel="CAD",ylabel="Burnt Mass Fraction")
	plt2 = plot(θ, dxb,xlabel="CAD",ylabel="Burnt Mass Rate")
	plt3 = plot(θ, press,xlabel="CAD",ylabel="In cylinder pressure $(data_units_pr)")
	plt = plot(plt1,plt2,plt3, layout=(3,1), legend = false, size=(600,600), width =2)
	return plt
end

# ╔═╡ db421954-6e41-11eb-1933-e13d8c281d97
plot_data()

# ╔═╡ 3c2f3e1e-45d3-11eb-2ced-0913e04fef87
md"### _Wiebe Fit Algorithm_

#### 1. Estimation of Start of Ignition (SOI)

There are a few methodologies for locating the start of ignition. 

###### A. Locating the zeroing of pressure second derivative

Here, SOI is defined as the CAD where the second derivative of pressure changes from negative to zero, immediately after the fuel injection timing which is $(@bind Inj_Tim NumberField(-20:0.01:20, default=-5.5)) (CAD)
"

# ╔═╡ 1e25f538-6de3-11eb-2387-c99790a03f43
function funarrzero(array)
	root = -1
	eps = 10^-7
	for i =1:length(array)-1
		if abs(array[i])<=eps
			root = i
			break
		else
			if array[i]*array[i+1]<0
				if abs(array[i])<=abs(array[i+1])
					root = i
				else
					root = i+1
				end
			end
		end
	end
	return root
end;

# ╔═╡ 1f3541a0-45d3-11eb-3e7d-c7d551cdf958
begin
	spl_pr =  Spline1D(θ, press; w=ones(length(θ)), k=4, bc="nearest", s=0);
	#Calculation of 2nd Derivative of pressure
	dp1 = derivative(spl_pr, θ)
	dp2 = derivative(spl_pr, θ; nu = 2)
	
	# Find where the 2nd der. is zero
	SOI_pr_idx_min = findmin(abs.(θ.-Inj_Tim))[2]; # SOI is after injection
	SOI_pr_idx_max = findmax(dp2)[2] # SOI is before max dp2
	
	th0_pr_idx = argminima(abs.(dp2[SOI_pr_idx_min:SOI_pr_idx_max]))[2]+SOI_pr_idx_min
	th0_pr = θ[th0_pr_idx]
end;

# ╔═╡ d0c1bd1c-45cf-11eb-147d-85df770e4cc3
function plot_pr_der()
	plt1 = plot(θ, dp1,xlabel="CAD",ylabel="$data_units_pr /CAD",title = "Pressure 1st derivative")
	plt2 = plot(θ, dp2,xlabel="CAD", ylabel="$data_units_pr /CAD^2",title = "Pressure 2nd derivative")
	scatter!([th0_pr], [0.0])
	plt = plot(plt1, plt2, layout=(2,1), legend = false, width =2 )
	return plt 
end;

# ╔═╡ ab8c336a-708b-11eb-0fc4-074924b77e01
plot_pr_der()

# ╔═╡ cbcd7732-6de9-11eb-39ce-670abe47ef15
md" ###### B. Locating the first positve Heat Realese Rate Value

It is obvious that the fisrt positve point of the heat realese rate curve is the Start of Ignition (SOI). 
"

# ╔═╡ 90c0c57a-6df0-11eb-00dc-b75aa76b5498
begin
	th0_hrr_idx = -1
	th0_hrr = -1
	for i=1:length(dxb)
		if dxb[i]>0
			th0_hrr_idx = i-1
			th0_hrr = θ[i-1]
			break
		end
	end
	th0_hrr
end

# ╔═╡ b629557e-708a-11eb-1d0c-6fa2fd31003d
md" Using the above methods the followeing vlaues occured which can be also seen at the graph bellow. Select which are to be taken into account.

**SOI pressure method** $\theta_0$ = $(th0_pr) $(@bind check_SOI_pr CheckBox(default=true)) 


**SOI HRR method** $\theta_0$ = $(th0_hrr) $(@bind check_SOI_hrr CheckBox(default=true))"



# ╔═╡ 6bd2fb94-708d-11eb-0583-a5ac59a10447
θ0 = (th0_pr*check_SOI_pr + th0_hrr*check_SOI_hrr)/(check_SOI_pr+check_SOI_hrr)

# ╔═╡ 44580476-6df2-11eb-02eb-55090c83d2b6
function plot_SOI()
	plot(θ,dxb,xlabel="CAD",ylabel="Heat Release Rate",title = "Burnt Mass Fraction Rate", linewidth =2)
	scatter!([th0_hrr], [0.0])
	scatter!([θ0], [0.0])
	plt = scatter!([th0_pr], [0.0])
	return plt
end;

# ╔═╡ 8ca5d232-6e42-11eb-3222-45b159f8b693
plot_SOI()

# ╔═╡ ddde8d88-cd57-4a5d-94d1-118635532b91
md"#### 2. Single Wiebe linearization and fitting 

Now the single wiebe function is linearized as follows:

$$x_b=  1 - exp \left[ 
-\alpha\left( \frac{\theta-\theta_0} {\Delta\theta} \right) ^{m+1} 
\right] \Rightarrow$$

$$G = (m+1)^{-1}[H-ln(a)]+ln(Δθ)$$
$$G = ln(φ-φ_0)$$
$$Η = ln[-ln(1-x_b(φ))]$$

Assuming that 

$$A = m+1$$
$$B = lnα-(m+1)ln(Δφ)$$
$$x = ln(φ-φ_0)$$
$$y = ln[-ln(1-xb)]$$

We fit the data to the linear equation y=Ax+B

"

# ╔═╡ 80961c54-2686-46fc-a0bb-e927e38f073d
plot(

# ╔═╡ aaea78cc-4ba3-4f51-9106-76550b49507c
θ[θ.-θ0 .>0]

# ╔═╡ b124a4c9-8a0a-4139-902e-bcfa9a871113


# ╔═╡ db173ccc-7144-11eb-2838-ab394e2f237e
begin
	#Calculation of 2nd and 3rd derivative of pressure
	spl_dxb = Spline1D(θ, dxb; w=ones(length(θ)), k=4, bc="nearest", s=0);
	dxb2 = derivative(spl_dxb, θ)
	dxb3 = derivative(spl_dxb, θ; nu = 2)
	dxb3_max_idx = argmaxima(dxb3)
	dxb3_sel = string.(collect(1:length(dxb3_max_idx)))
	
	spl_hrr_max = findmax(spl_dxb(θ))[1]
	spl_hrr_max_idx = findmax(spl_dxb(θ))[2]
end

# ╔═╡ e7cec92c-7136-11eb-1e59-0d56fc015cca
md"#### 3. Estimation of transition angles. 

The combustion is considered to be happen in multiple phases. These phases can be seperated by locating the local maximum of the third derivative of burnt mass fraction. As it can be shown, there are a few maxima. Select the appropriate point:

$(@bind dxb3_0 Select(dxb3_sel)) 
"

# ╔═╡ df8b1074-71d9-11eb-28ec-6fcc3fa9a294
md" $\theta_{0,2}$ = $(th2_in = θ[dxb3_max_idx[findall(x->x==dxb3_0, dxb3_sel)]])
"

# ╔═╡ 7a02d25c-713c-11eb-01c9-b364b27a165c
begin

	plot(θ,dxb, label = L"\x_b'")

end

# ╔═╡ f9df1a06-7094-11eb-0281-2bad62653c94
md"#### 3. Wiebe shape parameters estimation

Now having split the data into 2 Sections 



"

# ╔═╡ 109c1934-712c-11eb-1c86-43bfced1f906
begin
	#Filtered burnt mass fraction
	spl_xb = cumul_integrate(spl_theta2,spl_hrr(spl_theta2));
	indx1 = dxb3_max_idx[findall(x->x==dxb3_0, dxb3_sel)][1]
	
	#New variables
	psi = log.(spl_theta2.-th0_hrr)[50:end]
	a_in = 6.908
	xi = log.(- 1 ./a_in.*log.(1 .-abs.(spl_xb)))[50:end]
	
	psiD = ones(length(xi)-1)
	for i =1:length(xi)-1
		psiD[i] = (psi[i+1]-psi[i])/(xi[i+1]-xi[i])
	end
end

# ╔═╡ dd7a0a66-71dd-11eb-3326-03225651dc72
begin
	@. wlin(x,p) = p[1]*x+p[2]
	
	lbl = 1.0*[0.01, 0]
	ubl = 1.0*[1, 10]
	pl = 1.0*[0.5, 4]

	#1st Wiebe
	#coefs, converged, iter = nonlinear_fit(model_wiebe_1, theta_dat_hr, xb, p0,)
	fitlin1 =  curve_fit(wlin, xi[1:indx1-1000], psi[1:indx1-1000], pl)
	
	lam, bet = coef(fitlin1)
	
	m1_in = 1/lam-1
	dtheta1_in = exp(bet)
	
	xig1 = xi[1:indx1-100]
	yi1 = lam*xig1 .+bet
	
	#2nd Wiebe
	fitlin2 =  curve_fit(wlin, xi[indx1+100:end], psi[indx1+100:end], pl)
	
	lam, bet = coef(fitlin2)
	
	m2_in = 1/lam-1
	dtheta2_in = exp(bet)
	
	
	xig2 = xi[indx1+100:end]
	yi2 = lam*xig2 .+bet
	
	
	yi_sep =  log(spl_theta2[indx1]-th0_hrr)
	xi_sep =  log(- 1/a_in*log(1-spl_xb[indx1]))
	
	lamin = spl_xb[indx1]
	
end

# ╔═╡ 9af0a6aa-7b47-11eb-2909-052440774dcb
coef(fitlin2)

# ╔═╡ 9adea2cc-712c-11eb-2a8d-d7ab17b66c2a
begin
	plot(xi,psi, label = "Transformed Wiebe")
	plot!(xig1, yi1)
	plot!(xig2, yi2)
	scatter!([xi_sep],[yi_sep],legend=:outerbottomright)
end

# ╔═╡ ee6a7ca8-6345-11eb-0b83-b97b40c731e3
md"### Double Wiebe Function Fitting"

# ╔═╡ 4e606bb4-4202-11eb-228a-5521a8fe08b1
#Single Wiebe Function
wiebe_1(theta) = (1 - exp(-a0*(((theta-theta00)/Dtheta0)*(theta-theta00>=0))^(m0+1)))

# ╔═╡ 5744af30-4202-11eb-277b-f9ef05fd00f2
#Single Derivative Wiebe Function
dwiebe_1(theta) = ((a0*(m0+1))/(Dtheta0))*(((theta-theta00)/(Dtheta0)*(theta-theta00>=0))^m0)*exp(-a0*(((theta-theta00)/Dtheta0)*(theta-theta00>=0))^(m0+1))

# ╔═╡ bd76bef0-4201-11eb-0618-9bdc4c1022e9
md"#####  Wiebe Function Parameters Estimation_
Using the above methodologies, the parameters of double wiebe function were estimated the the ones below.

| Parameter  | Value | 
| ------ | ------ |
| _Start of 1^{st} Combustion Phase \theta_{0,1}_ | $(@bind th01in NumberField(-30:0.001:30, default=theta0)) |
| _Start of 2^{nd} Combustion Phase \theta_{0,2}_ | $(@bind th02in NumberField(-30:0.001:30, default=th2_in[1])) |
| _Efficiency Factor a_1_ | $(@bind a1in NumberField(0:0.0001:30, default=a_in)) |
| _Efficiency Factor a_2_ | $(@bind a2in NumberField(0:0.0001:30, default=a_in)) |
| _Shape Factor m_1_ | $(@bind m1in NumberField(0:0.0001:20, default=m1_in)) |
| _Shape Factor m_2_ | $(@bind m2in NumberField(0:0.0001:20, default=m2_in)) |
| _Duration of Phenomenon 1_ | $(@bind Dth1in NumberField(0:500, default=dtheta1_in)) |
| _Duration of Phenomenon 2_ | $(@bind Dth2in NumberField(0:500, default=dtheta2_in)) |
| _\lambda_ | $(@bind lambdain NumberField(0:0.001:1, default=lamin)) |

"

# ╔═╡ 76d6d118-6099-11eb-15cd-692e5cbe0502
theta = collect(-20:0.01:70)

# ╔═╡ 885ce27c-6343-11eb-1d95-db43ff2e4634
begin

	
	@. mbw2(x,p) = p[5]*(1 - exp(-p[1]*(((x-p[3])/p[4])*(x-p[3]>=0))^(p[2]+1)))+(1-p[5])*(1 - exp(-p[6]*(((x-p[8])/p[9])*(x-p[8]>=0))^(p[7]+1)))
	
	#@. mbw2(x,p) = p[3]*(1 - exp(-6.908*(((x-th01in)/p[2])*(x-th01in>=0))^(p[1]+1)))+(1-p[3])*(1 - exp(-6.908*(((x-p[5])/p[6])*(x-p[5]>=0))^(p[4]+1)))
	
	
	           #a1     m1     θ1  	  Δθ1    λ      α2        m2       θ2    Δθ2
	lb2 = 1.0*[0.001,  0.001,  -10,  10,      0.001,     0.001,  0,   -10,   10]
	ub2 = 1.0*[10,    50,    10, 	  100,     1,        10,     50,   30,   100]
	p02 = 1.0*[a1in,  m1in,  th01in,  Dth1in,  lambdain, a_in,   m2in, th02in,   Dth2in]
	
	
	#lb2 = 1.0*[0.001,  2,      0.001,       0,   -10,   2]
	#ub2 = 1.0*[50,     200,     1,          50,   20,   200]
	#p02 = 1.0*[m1in,   Dth1in,  lambdain,   m2in, th02in,   Dth2in]

	w1=ones(length(spl_xb))
	
	w1[indx1-50:indx1+50].= 1000
	w1[1]= 1000
	
	#coefs, converged, iter = nonlinear_fit(model_wiebe_1, theta_dat_hr, xb, p0,)
	fit2 =  curve_fit(mbw2, spl_theta2, spl_xb, w1,  p02, lower = lb2, upper = ub2)
	
	a1, m1, theta01, Dtheta1, lambda, a2, m2, theta02, Dtheta2  = coef(fit2)
	
	#m1, Dtheta1, lambda, m2, theta02, Dtheta2  = coef(fit2)
	#a1 = 6.908
	#a2 = 6.908
	#theta01 = th01in
	
end

# ╔═╡ 56ce7216-4202-11eb-3a1e-93b09ab25c32
#Double Wiebe Function
wiebe_2(theta) = lambda*(1 - exp(-a1*(((theta-theta01)/Dtheta1)*(theta-theta01>=0))^(m1+1)))+(1-lambda)*(1 - exp(-a2*(((theta-theta02)/Dtheta2)*(theta-theta02>=0))^(m2+1)))

# ╔═╡ 57dbb0f4-4202-11eb-06c6-6d6d395a5588
#Single Derivative Wiebe Function
dwiebe_2(theta) = lambda*(((a1*(m1+1))/(Dtheta1))*(((theta-theta01)/(Dtheta1)*(theta-theta01>=0))^m1)*exp(-a1*(((theta-theta01)/Dtheta1)*(theta-theta01>=0))^(m1+1)))+(1-lambda)*(((a2*(m2+1))/(Dtheta2))*(((theta-theta02)/(Dtheta2)*(theta-theta02>=0))^m2)*exp(-a2*(((theta-theta02)/Dtheta2)*(theta-theta02>=0))^(m2+1)))

# ╔═╡ 4449b2a0-8a4d-11eb-36ac-3d6c73b7197f


# ╔═╡ d6406e62-88c3-11eb-14be-2bd17c85abbf
spl_theta2[indx1]

# ╔═╡ 475df1bc-76df-11eb-26ff-7be9cf36cf1f
md"###### Final Values

| Parameter  | Value | 
| ------ | ------ |
| _Start of 1^{st} Combustion Phase \theta_{0,1}_ | $(theta01) |
| _Start of 2^{nd} Combustion Phase \theta_{0,2}_ | $(theta02) |
| _Efficiency Factor a_1_ | $(a1) |
| _Efficiency Factor a_2_ | $(a2) |
| _Shape Factor m_1_ | $(m1) |
| _Shape Factor m_2_ | $(m2) |
| _Duration of Phenomenon 1_ | $(Dtheta1) |
| _Duration of Phenomenon 2_ | $(Dtheta2) |
| _\lambda_ | $(lambda) |

"

# ╔═╡ 763ca7cc-8a4d-11eb-12ce-11d0b64a26f8
1-0.716

# ╔═╡ e9fef324-6345-11eb-19c0-a137a92b02d1
begin
	xb_m2 = wiebe_2.(theta)
	dxb_m2 = dwiebe_2.(theta)
end

# ╔═╡ 0fa88214-6346-11eb-3a74-835f8a7d28a4
begin
	plot(theta_dat_hr, xb, xlabel="CAD", ylabel="Burned Mass Fraxction", title = "Modelling Results")
	plot!(theta, xb_m2)
end

# ╔═╡ 35613566-6346-11eb-22c2-470e24eebc10
begin
	plot(theta_dat_hr, dxb, xlabel="CAD", ylabel="Burned Mass Fraxction", title = "Modelling Results")
	plot!(theta, dxb_m2)
end

# ╔═╡ 97bb8dc8-4371-11eb-06cc-59ffba3bc103
#begin
#	xx = hcat(theta_dat_hr,xb)
#	
#	@. model_wiebe_1(x,p) = (1-exp(-p[1]*(((x[1]-p[2])/p[3])*((x[1]-p[2])>0))^(p[4]+1)))-x[2]
#	
#	#a = p[1],theta00 = p[2], Dtheta0 = p[3], m0 = [4]
#
#	#lb = 1.0*[-Inf, -Inf, 0.1, 0]
#	#ub = 1.0*[Inf, Inf, Inf, Inf]
#	p0 = 1.0*[3, 0, 50, 2]
#	
#	eps = 0.000001 
#	maxiter=200 
#	
#	#coefs, converged, iter = nonlinear_fit(model_wiebe_1, theta_dat_hr, xb, p0,)
#	coefs, converged, iter =  nonlinear_fit(xx, model_wiebe_1, p0, eps, maxiter)
#	
#	a0, theta00, Dtheta0, m0 = coefs
#end

# ╔═╡ aa6e46e0-6dd3-11eb-1cc9-894b66955ebe
md" ### Bibliography:

[2] _Rakopoulos CD, Antonopoulos KA, Rakopoulos DC. Experimental heat release
analysis and emissions of a HSDI diesel engine fueled with ethanol-diesel fuel
blends. Energy 2007;32(10):1791–808._"

# ╔═╡ 65b8cc60-8d26-489c-ab75-3381f47398fa
xx = xb[θ.-θ0 .>0]

# ╔═╡ 6f419661-7330-4a08-a71c-839581ca17f6
begin
	
	φ = θ[θ.-θ0 .>0]
	xx = xb[θ.-θ0 .>0]
	
	χ = log.(φ .- θ0)
	ψ = log.(-log.(1 .-xb))
	
end
	
	

# ╔═╡ Cell order:
# ╟─7a243fdc-405e-11eb-3e37-cfda2218520d
# ╟─79d3def6-405f-11eb-0e3c-43105efaef03
# ╟─a7461f68-41ed-11eb-0a48-c1f7c41834a5
# ╟─2d2cf1bc-4200-11eb-3771-61e3cbcedf90
# ╟─38168156-4200-11eb-1428-df9a090087c2
# ╟─fed92118-729a-11eb-11a5-e5267ce84380
# ╟─5763f37c-4200-11eb-2c4c-01586566124c
# ╟─5985dcd0-4201-11eb-3637-4ba209c60c01
# ╟─42af1e4c-634e-11eb-0a08-bd55c6c2f851
# ╟─52a7eacc-4201-11eb-0e52-0d5343c651af
# ╟─db421954-6e41-11eb-1933-e13d8c281d97
# ╟─3c2f3e1e-45d3-11eb-2ced-0913e04fef87
# ╟─1e25f538-6de3-11eb-2387-c99790a03f43
# ╟─1f3541a0-45d3-11eb-3e7d-c7d551cdf958
# ╟─d0c1bd1c-45cf-11eb-147d-85df770e4cc3
# ╟─ab8c336a-708b-11eb-0fc4-074924b77e01
# ╟─cbcd7732-6de9-11eb-39ce-670abe47ef15
# ╟─90c0c57a-6df0-11eb-00dc-b75aa76b5498
# ╟─b629557e-708a-11eb-1d0c-6fa2fd31003d
# ╠═6bd2fb94-708d-11eb-0583-a5ac59a10447
# ╠═44580476-6df2-11eb-02eb-55090c83d2b6
# ╠═8ca5d232-6e42-11eb-3222-45b159f8b693
# ╟─ddde8d88-cd57-4a5d-94d1-118635532b91
# ╠═6f419661-7330-4a08-a71c-839581ca17f6
# ╠═80961c54-2686-46fc-a0bb-e927e38f073d
# ╠═65b8cc60-8d26-489c-ab75-3381f47398fa
# ╠═aaea78cc-4ba3-4f51-9106-76550b49507c
# ╠═e7cec92c-7136-11eb-1e59-0d56fc015cca
# ╠═df8b1074-71d9-11eb-28ec-6fcc3fa9a294
# ╠═b124a4c9-8a0a-4139-902e-bcfa9a871113
# ╠═db173ccc-7144-11eb-2838-ab394e2f237e
# ╠═7a02d25c-713c-11eb-01c9-b364b27a165c
# ╠═f9df1a06-7094-11eb-0281-2bad62653c94
# ╠═109c1934-712c-11eb-1c86-43bfced1f906
# ╠═dd7a0a66-71dd-11eb-3326-03225651dc72
# ╠═9af0a6aa-7b47-11eb-2909-052440774dcb
# ╠═9adea2cc-712c-11eb-2a8d-d7ab17b66c2a
# ╟─ee6a7ca8-6345-11eb-0b83-b97b40c731e3
# ╠═4e606bb4-4202-11eb-228a-5521a8fe08b1
# ╠═5744af30-4202-11eb-277b-f9ef05fd00f2
# ╟─56ce7216-4202-11eb-3a1e-93b09ab25c32
# ╠═57dbb0f4-4202-11eb-06c6-6d6d395a5588
# ╠═bd76bef0-4201-11eb-0618-9bdc4c1022e9
# ╟─76d6d118-6099-11eb-15cd-692e5cbe0502
# ╠═885ce27c-6343-11eb-1d95-db43ff2e4634
# ╠═4449b2a0-8a4d-11eb-36ac-3d6c73b7197f
# ╠═d6406e62-88c3-11eb-14be-2bd17c85abbf
# ╟─475df1bc-76df-11eb-26ff-7be9cf36cf1f
# ╠═763ca7cc-8a4d-11eb-12ce-11d0b64a26f8
# ╟─e9fef324-6345-11eb-19c0-a137a92b02d1
# ╟─0fa88214-6346-11eb-3a74-835f8a7d28a4
# ╟─35613566-6346-11eb-22c2-470e24eebc10
# ╟─97bb8dc8-4371-11eb-06cc-59ffba3bc103
# ╟─aa6e46e0-6dd3-11eb-1cc9-894b66955ebe
