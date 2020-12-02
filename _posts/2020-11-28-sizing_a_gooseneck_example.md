# Compressible Flow Example - Sizing a goose neck vent
 
When determining the required venting for aboveground storage tanks it is typical to calculate normal and emergency vent rates using standards such as API 2000, which gives the venting as an equivalent flowrate of air at standard state. Most off-the-shelf vents are sized in terms of pressure drop and SCFH through it, so after calculating the required venting one can simply buy a vent with the needed characteristics. For example for an emergency vent from Morrison brothers:

![image.png](/images/sizing_a_gooseneck_example_files/att1.png)


![image.png](/images/sizing_a_gooseneck_example_files/att2.png)

It's not uncommon, though, for tanks to have gooseneck vents constructed from piping. This is fairly normal for tanks holding water, or other non-volatile substances, where the tank is open to atmosphere and there are no dangerous vapours that need to be managed. The gooseneck itself is merely to keep rain, and wildlife, out of the tank.

Sizing the gooseneck to match the required venting involves performing some simple compressible flow calculations, which is fairly straight-forward to set up in a generalized way such that, beyond this motivating example, it can be extended to lots of other problems. Though it isn't at all uncommon, in this particular case, for the flow calculations to be done assuming an incompressible fluid as, over the length of the gooseneck, the pressure drop is typically slight and the compressibility of the gas (air usually) is not important.

With that in mind, I am going to work through the problem in stages of escalating complexity, very likely the most complicated (fanno flow) is overkill for this specific example but it's worth putting it all down as the same tools can be used for compressible flow calculations through many piping situations.

The primary references I am using are:
+ *TP410M Flow of Fluids*, Crane Co., Stamford, Connecticut, 2013
+ *Perry's Chemical Engineers' Handbook, 8th ed.*


## The Scenario

Suppose an atmospheric storage tank with a normal venting requirement, as calulated from API 2000 or the like, of 200&times;10³ SCFH and a max pressure of 1 psig. We wish to design a gooseneck vent that can handle that level of venting. Suppose that the gooseneck design we have in mind is a vertical length of pipe extending up from the tank roof, two 90° bends, and an exit covered with a mesh screen (to keep birds from nesting in it, yes this is a thing). The gooseneck is a constant diameter of pipe throughout.

For notation, the flow begins at the pipe entrance (1) and ends at the exit (2).

![image.png](/images/sizing_a_gooseneck_example_files/att3.jpg)

For the pipe bends the bend radius to diameter ratio needs to be specified, I'm going to suppose $r/D = 1.5$. Another important parameter is the pipe roughness, $\epsilon$, which for commercial steel is $\epsilon = 0.0457 \mathrm{mm}$[^1]. At this point I could specify a length for the straight section of pipe, a fixed height above the tank roof that is independent of the final chosen diameter of the piping, or I could fix a design and scale the whole vent up and down as required. For simplicity I am going to assume 3ft of piping.

At this point I am going to set up the equations with no knowledge of what the final pipe diameter $D$ will be, then numerically solve for the minimum diameter that meets the requirements. The actual diameter will be the next largest NPS size pipe.

[^1]: *Perry's* pg 6-10 table 6-1


```julia
# Setting up some convenient unit conversions
SCFH = uconvert(u"m^3/s", 1u"ft^3/hr") # unit conversion scfh -> sm^3/s
psi = uconvert(u"Pa", 1u"psi")         # unit conversion psi->Pa
ft = uconvert(u"m", 1u"ft")            # unit conversion ft->m
inch = uconvert(u"m", 1u"inch")        # unit conversion inch->m
mm = uconvert(u"m", 1u"mm")            # unit conversion mm->m

# Given in the scenario, now converted to SI
Q = 200e3SCFH
pₐ = 14.696psi
pₘₐₓ = 1psi + pₐ
L = 3ft
ϵ = 0.0457mm
```


I am assuming, for simplicity, that ambient conditions are standard conditions.


```julia
# Some useful physical properties of air
R = 8.31446261815324u"Pa*m^3/mol/K" # Universal gas constant to more digits than are at all necessary
Tₐ = 298.15u"K"                     # Standard conditions, 15°C 
Mw = 0.02896u"kg/mol"               # Molar weight of air, from Perry's
k = 1.4                             # Ratio of heat capacities, Cp/Cv, ideal gas

ρ(p, T) = (p * Mw)/(R * T)          # density of air, ideal gas law, kg/m^3
μ(T) = 1u"Pa*s"*(1.425e-6*ustrip(u"K",T)^0.5039)/(1 + 108.3/ustrip(u"K",T)) # viscosity of air, from Perry's
```


### Frictional head loss

Regardless of the method of performing the compressible flow calculations, the frictional head loss in the piping needs to be accounted for. I am using the K factor method as it is convenient and K factors are tabulated for most everything in references such as Crane's TP-410. One thing to be very careful with is the difference between the Darcy and Fanning friction factors. I am using Crane's where everything is in terms of the Darcy friction factor, which is 4&times; the Fanning friction factor, but Perry's defaults to the Fanning friction factor.

From Crane's I have the following K factors for each piece of the gooseneck[^2]
+ Entrance - $K_1 = 0.5$
+ Vertical Pipe - $K_2 = f \frac{L}{D}$
+ First 90° bend - $K_3 = 14 f_T$
+ Second 90° bend - $K_4 = 14 f_T$
+ Mesh screen - $K_5 = f_T$
+ Exit to atmosphere - $K_6 = 1.0$

Where $f$ is the Darcy friction factor, $f_T$ is the turbulent friction factor. I am assuming the entrance to the vent is sharp edged, and the K factors for the bends are for bends with $r/D = 1.5$.

For some notational convenience I am going to define the relative roughness $\kappa = { \epsilon \over D }$ and the reduced length $l = {L \over D}$ so that, along with the Reynolds number $Re$, the K factors are in terms of dimensionless numbers only.

The Darcy friction factor generally depends on which regime the flow is in, laminar, transitional, or turbulent. Since I don't want to be referring to a Moody diagram I want to use a single equation that operates over a wide range of Reynolds numbers and potentially in laminar and transitional flow regimes since I don't know *a priori* what the flow in the vent will be. There are equations like the Serghide correlation or Churchill correlation that attempt to fit a Moody diagram but in a more convenient to use manner.

[^2]: See Crane's page A-30 for the K factors for pipe entrances, exits, and pipe bends



![svg](/images/sizing_a_gooseneck_example_files/output_7_0.svg)



The black line conservatively takes the max of either the laminar or turbulent friction factor in the transitional region $2100 \le Re \le 4000$. The Churchill correlation fits the general curve well for both the turbulent and laminar region, and provides reasonable values in the transitional region.

The Churchill correlation is[^3]

    $$ A = \left( 2.457 \ln\left( {1 \over {\left( \frac{7}{Re} \right)^{0.9} + 0.27\kappa} } \right) \right)^{16} $$
    
    $$ B = \left( \frac{37530}{Re} \right)^{16} $$
    
    $$ f = 8 \left( \left( \frac{8}{Re} \right)^{12} + { 1 \over {\left( A+B \right)^{3/2} } } \right)^{1/12} $$

The turbulent friction factor is the friction factor at fully turbulent flow, when *f* is no longer dependent upon the Reynolds number.

$$ f_T = { 0.25 \over { \left( \log \left( \kappa \over 3.7 \right) \right)^2 } } $$

With these defined I can write a function that gives $\sum_j K_j$ for any $\kappa$, *l*, and *Re*

[^3]: *Perry's* pg 6-11 equation 6-40, the equation here is given in terms of the darcy friction factor


```julia
function f(κ, Re)
    A = (2.457 * log(1/((7/Re)^0.9 + 0.27*κ)))^16
    B = (37530/Re)^16
    return 8*((8/Re)^12 + 1/(A+B)^(3/2))^(1/12)
end

fT(κ) = 0.25/log10(κ/3.7)^2

ΣK(κ,l,Re) = 0.5 + f(κ,Re)*l + 14*fT(κ) + 14*fT(κ) + fT(κ) + 1.0
```


### The Reynolds number

Since the vent has a constant cross sectional area, the mass velocity, $G$, is constant throughout

    $$ G = \frac{4 \dot{m} }{\pi D^2} $$

Where $\dot{m}$ is the mass flow rate in kg/s flowing through the vent, which is
$$\dot{m} = Q_{std} \cdot \rho \left(p_{std}, T_{std} \right) $$
with $Q$ the flow at standard conditions and $\rho$ the density at standard conditions.

*Note* The required vent flow is given in SCFH, this is not temperature or pressure dependent and $\dot{m}$ is a constant. If the flow given was a true volumetric flowrate, then $\dot{m}$ would be a function of temperature and pressure in general. This is one of those things that routinely snags inexperienced engineers as a flow in terms of a standard volume *looks* like a volumetric flowrate, it's in units of volume, but really isn't one since the temperature and pressure are set to standard state by definition.

The Reynolds number in terms of the mass velocity is

    $$ Re = \frac{G D}{\mu} $$

The only parameter in the Reynolds number which is not a constant is the viscosity, $\mu$, which is mostly dependent upon temperature and not pressure. So, to a very good first approximation, the Reynolds number is only a function of temperature. Which is very convenient.


```julia
m = Q*ρ(pₐ, Tₐ)      # mass flow, recall Q is at standard state so is not a function of temperature or pressure

G(D) = (4*m)/(π*D^2)

Re(D, T) = upreferred(G(D)*D/μ(T)) # upreferred promotes derived units to combos of base units, e.g. Pa -> kg/m/s^2
```


## Incompressible Flow

For problems, like this vent, where the pressure drop is expected to be small it is not unreasonable to assume the flow is approximately incompressible. This is very often what is done since, typically, it does not require any iterative methods and one can solve for incompressible flow directly. It is also useful to do even if one plans on performing a more complete compressible flow calculation, since this provides something of a sanity check and can be a good place to start compressible flow iterative calculations, i.e. the *initial guess* is from the incompressible case

To check whether or not the incompressible assumption is reasonable, consider the ratio of density inside the tank (at the max allowable pressure) to the density outside the tank, assuming ambient temperature. For an ideal gas this is

    $$ { \rho_1 \over \rho_2 } = { {p_1  Mw} \over {R  T_a} } { {R T_a} \over {p_2 Mw} } = \frac{p_1}{p_2} = \frac{p_{max} }{p_a}$$


```julia
pₘₐₓ/pₐ
```

    1.0680457267283614


Typically flows are considered incompressible if the density varies by less than ~5-10%, so this example (where the density varies by ~7%) is right in that range. You could justify it either way and it's more a function of how accurate the calculations need to be. Since, ultimately, we are solving for the pipe diameter and choosing the next largest pipe size it's probably fine to use an incompressible flow assumption. If anything the incompressible flow assumption will overestimate the pressure drop and thus lead to an over-sized pipe (erring on the side of caution)

The mechanical energy balance for an incompressible fluid is[^4]

    $$ p_1 - p_2 = \alpha_2 \frac{\rho v_2^2}{2} - \alpha_1 \frac{\rho v_1^2}{2} + \rho g \left( h_2 - h_1 \right) + \sum_j K_j \frac{\rho v^2}{2} $$

With the following simplifications
+ given the assumption of incompressible flow and a vent with a constant cross-sectional area $v_1 = v_2 = v$, 
+ the flow is uniform througout $\alpha_1 = \alpha_2 = 1.0$ 
+ the contribution due to hydrostatic pressure is negligible as the gas density is very small, $\rho g \left( h_2 - h_1 \right) \approx 0$

This becomes 

    $$ p_1 - p_2 = \sum_j K_j \frac{\rho v^2}{2} $$

Where the velocity, $v$ is

    $$ v = \frac{Q}{A} = { Q \over { \frac{\pi}{4} D^2 } } $$

Which can be solved algebraically for $D$, where $\rho$ is taken at the average pressure. That said it is easier to solve it numerically.

[^4]: *Perry's* pg 6-16 equation 6-90


```julia
v(D) = Q / ((π/4)*D^2)

Dᵢₙ(Dₗ, Dᵤ) = find_zero( D -> pₘₐₓ - pₐ - 0.5*ΣK(ϵ/D,L/D,Re(D,Tₐ))*ρ((pₘₐₓ + pₐ)/2, Tₐ)*v(D)^2, (Dₗ, Dᵤ), Roots.A42())

D0 = Dᵢₙ(4inch, 12inch)

uconvert(u"inch", D0)
```


    6.443404339132202 inch


At this point one would typically stop for this example, compressible flow calculations are probably unnecessary.

## Compressible Flow

In general compressible flow situations can be very difficult to solve, since the density of the working fluid is a function of the pressure and temperature and the pressure and temperature are varying throughout, which means heat transfer must also be accounted for in some way. There are several key simplifying assumptions that take this from the sort of problem solved with CFD to something actually quite simple. The first is to assume an ideal gas, the second is to examine two *extreme* cases of heat transfer: isothermal flow and adiabatic flow.

At these two extremal cases, the first where heat transfer is instantaneous and the second where it doesn't occur at all, provide bounds on the problem.

### Isothermal Compressible Flow

Isothermal compressible flow of an ideal gas is fairly straight forward. As already mentioned the Reynold's number depends only on temperature, which is constant by definition, so the Reynold's number is a constant. This means the frictional head loss is also constant throughout, and it is a simple matter to calculate the pressure drop.

The assumption that the flow is isothermal is very reasonable in this case. We are assuming normal venting from a tank at thermal equilibrium with it's surroundings, that is that the air flowing through the vent starts and ends in reservoirs of equal temperature. As gases expand the temperature decreases but the pressure drop across the vent is small so this effect should be negligible.

A quick check is to estimate the ratio of temperatures at the start and end of the vent assuming a frictionless adiabatic expansion
    $$ p_1^{1-k} T_1^k = p_2^{1-k} T_2^k = \mathrm{const}$$

    $$ \frac{T_2}{T_1} = \left( p_1 \over p_2 \right)^{ {1-k} \over k}$$


```julia
(pₘₐₓ/pₐ)^((1-k)/k)
```

    0.9813670503935878


So we expect even in the most extreme case the temperature change is ~2%, justifying the assumption that the venting is isothermal.

The isothermal flow of an ideal gas going through a length of piping is[^5]

    $$ p_{1}^{2} = G^{2} \frac{RT}{Mw} \left[ \sum \limits_{j} K_{j} + 2\ln \frac{p_{1} }{p_{2} } \right] + p_{2}^{2} $$

If we assume the system is at thermal equilibrium with the outside air, then $T = T_a$ and $p_2 = p_a$

The only unknown is *p<sub>1</sub>*, which can be solved for numerically.

[^5:] *Perry's* pg 6-23 equation 6-114, this equation also neglects changes in elevation


```julia
pᵢₜ(G,κ,l,Re) = find_zero(p -> p^2 - G^2 * (R*Tₐ/Mw) * (ΣK(κ,l,Re) + 2*log(p/pₐ)) - pₐ^2, pₘₐₓ) # solves the isothermal flow equation with an initial guess p₁=pₘₐₓ
```

At this point we can write a simple function to solve for the minimum diameter that meets our requirement that $p_1 \le p_{max}$.


```julia
Dᵢₜ(Dₗ, Dᵤ) = find_zero( D -> pₘₐₓ - pᵢₜ( G(D), ϵ/D, L/D, Re(D, Tₐ)), (Dₗ, Dᵤ), Roots.A42())

D1 = Dᵢₜ(4inch, 12inch)

uconvert(u"inch", D1)
```


    6.4377172085806516 inch


### Adiabatic (Fanno) Flow


Adiabatic flow of an ideal gas through a pipe, also called Fanno flow, is somewhat more difficult than isothermal flow -- there are more steps in the iterative solution as the temperature along the length of the vent changes and thus the Reynolds number changes. The general process starts by assuming a constant friction factor, calculating the pressure and temperature changes due to the adiabatic expansion of an ideal gas, adjusting the friction factor for the temperature change, and iterating until everything converges.

There are a few ways of setting up the calculations. We could assume the gas exits at ambient conditions -- both ambient temperature and pressure -- or assume the tank starts at thermal equilibrium with the environment but at a higher pressure and the gas exits at ambient pressure and some other temperature -- less than ambient due to adiabatic expansion. The first set of assumptions is in some ways easier to calculate, but the second set of assumptions is more physically realistic, and consistent with the assumptions made when solving the isothermal case.

One thing we should check before proceeding is whether or not the flow will be choked, essentially will the flow velocity reach $Ma = 1$, the following discussion assumes flow remains sub-sonic, and this is easy to check. The critical pressure, at which flow becomes sonic, is given by[^6]

    $$ { p^o \over p_1 } = \left(2 \over k+1 \right)^{k \over {k-1} } $$

with the criteria that flow is sub-sonic if

    $$ { p_2 \over p_1 } > { p^o \over p_1 } $$

[^6]: *Perry's* pg 6-23 equation 6-119


```julia
(pₐ / pₘₐₓ) > (2 / (k+1))^(k/(k-1))
```

    true


The basic relation of Fanno flow that drives the equations is the relationship between the Fanno parameter and the Mach number[^7]

    $$ Fa = \left( \frac{fL^{*} }{D} \right) = \frac{1 - Ma^{2} }{kMa^{2} } + \frac{k+1}{2k} \ln \left( \frac{ \left( k+1 \right) Ma^{2} }{ 2+\left( k+1 \right) Ma^{2} } \right) $$

Where I am defining $Fa$ to be the Fanno parameter. The Fanno parameter is calculated from some point in the flow path through to the critical point, where flow goes sonic. The critical point can be a hypothetical point, assuming the pipe is infinite, or it can be real. In this case I am assuming the flow within the vent will remain subsonic.

It is worth noting that elbows near the exit of a pipe may choke the flow even though the $Ma < 1$ due to the nonuniform velocity profile in the elbow. By the design of this gooseneck we know this is the case and should keep that in mind when evaluating the results.

For two points along a pipe, 1 and 2, the difference between their Fanno parameters is the frictional loss between those two points[^8]

    $$ Fa_1 - Fa_2 = \sum_{j} K_{j} $$

Where the $K_j$ are usually evaluated at the *average* temperature ${ {T_1 + T_2} \over 2}$.

The Mach number at some point *i* along the pipe, for an ideal gas, is given by[^9]

    $$ Ma_{i} = \frac{v}{c} = \frac{G}{p_{i} } \sqrt{ \frac{RT_{i} }{kMw} } $$

and the pressure can be calculated given a Mach number by re-arranging

    $$ p_{i} = \frac{G}{Ma_{i} } \sqrt{ \frac{RT_{i} }{kMw} } $$

and for any two points along the pipe the temperatures are related by[^10]

    $$ T_{1} = T_{2} \frac{2 + \left( k-1 \right) Ma_{2}^{2} }{2 + \left( k-1 \right) Ma_{1}^{2} } $$


Putting all of this together, the procedure for adiabatic ideal gas flow through the piping with a given diameter $D$ is:
1. Given $G$ calculate $Ma_2$ at ambient conditions, this is the initial guess for the exit conditions
1. Calculate $\sum_j K_j$ at ambient conditions, this is the initial guess for the frictional loss
1. Calculate $Fa_2$ with $Ma_2$
1. Calculate $Fa_1$ from $Fa_2$ and $\sum_j K_j$
1. Solve for $Ma_1$ given $Fa_1$, this is done numerically
1. Solve for $T_2$ given $Ma_1$ and letting $Ma_2$ vary with temperature, this is done numerically
1. Recalculate $Ma_2$ given $T_2$ and $\sum_j K_j$ at the average temperature ${ {T_1 + T_2} \over 2}$ and repeat from step 3
1. Continue to iterate until $p_1$ stops changing

While that looks complicated, each step is fairly easy. In my experience, with subsonic flow, this converges very quickly.

[^7]: *Perry's* pg 6-24 equation 6-123

[^8]: *Perry's* pg 6-24 equation 6-124

[^9]: derived for an ideal gas
    
    $$ G = \rho v = { {p  Mw} \over {R T} } v$$
    
    $$ c = \sqrt{ {k R T} \over Mw } $$
    
    $$ Ma = { v \over c } = G { {R T} \over {p  Mw} } \sqrt{ Mw \over {k R T} } = \frac{G}{p} \sqrt{ \frac{RT}{kMw} }$$


[^10]: *Perry's* pg 6-23 equation 6-116, taking two points and canceling out the stagnation temperature


```julia
Fa(Ma) = ((1-Ma^2)/(k*Ma^2)) + ((k+1)/(2k))*log( ((k+1)*Ma^2) / (2 + (k+1)*Ma^2))

Ma(G, p, T) = upreferred((G/p)*√((R*T)/(k*Mw)))  # upreferred(...) ensures units cancel appropriately and the Ma is unitless

T2(T₁, Ma₁, G, p) = find_zero(T₂ -> (2 + (k-1)*Ma₁^2)*T₁ - (2 + (k-1)*Ma(G, p, T₂)^2)*T₂, T₁)


function pfa(D)
    # pre-calculating diameter dependent variables
    Gᵢ = G(D)
    κᵢ = ϵ/D
    lᵢ = L/D

    # initial values
    T₂ = Tₐ
    Ma₁ = Ma(Gᵢ, pₐ, Tₐ)
    p₁ⁿᵉʷ = uconvert(u"Pa",(Gᵢ/Ma₁)*√((R*Tₐ)/(k*Mw)))
    
    # loop until the error is below the given tolerance, but don't loop forever!
    err, i = 1.0, 0
    rtol, max_count = 1e-9, 1e5
    while (err > rtol) && (i < max_count)
        # Starting up
        Tₐᵥ = 0.5*(Tₐ + T₂)
        Reᵢ = Re(D,Tₐᵥ)
        Ma₂ = Ma(Gᵢ, pₐ, T₂)
             
        # Steps 3 - 6
        Fa₂ = Fa(Ma₂)
        Fa₁ = Fa₂ + ΣK(κᵢ,lᵢ,Reᵢ)
        Ma₁ = find_zero(x -> Fa₁ - Fa(x), Ma₁)
        T₂ = T2(Tₐ, Ma₁, Gᵢ, pₐ)
        
        # Check if pressure has converged
        p₁ᵒˡᵈ = p₁ⁿᵉʷ
        p₁ⁿᵉʷ = uconvert(u"Pa",(Gᵢ/Ma₁)*√((R*Tₐ)/(k*Mw)))
        err = abs(p₁ⁿᵉʷ - p₁ᵒˡᵈ)/p₁ᵒˡᵈ
        i += 1
    end
       
    # if the loop failed to converge, let me know
    if i >= max_count
        println("iterations exceeded max count")
        println("remaining error is $err")
    end
    
    return p₁ⁿᵉʷ, T₂
end
```



```julia
Dfa(Dₗ, Dᵤ) = find_zero( D -> pₘₐₓ - pfa(D)[1], (Dₗ, Dᵤ), Roots.A42())

D2 = Dfa(4inch, 12inch)

uconvert(u"inch", D2)
```

    6.431779377501573 inch


## Minimum Diameter

At this point we have solved for the minimum vent diameter in three different ways and, more or less, got the same answer three times. The minimum diameter is ~6.4in ID for all cases and the next largest standard pipe size is 8in so regardless of the method, in this particular example, one arrives at the same final answer.

In general the incompressible model will always *overestimate* the pressure drop across the vent, leading to a larger vent size, and the adiabatic flow will provide an *underestimate*, the true minimum would be somewhere between the two. This is seen much more clearly at vent diameters less than ~5in where the pressure drop is more significant, more analogous to relief piping for a pressure vessel than venting for an atmospheric storage tank. Of course all of this is assuming flow remains subsonic, if the pressure drop leads to sonic flow then things are quite different.



![svg](/images/sizing_a_gooseneck_example_files/output_29_0.svg)



## Concluding Remarks

Often things like compressible flow can be intimidating since these problems, even in the simplified ideal gas case, require iterative solutions and often iterative solutions within iterative solutions. However, once the basic pieces are set up, compressible flow can be fairly simple to deal with. There are some pitfalls here that, if one wanted to create a nice generalized set of code, would have to be dealt with.

The big one being all the `find_zero()` calls rely on the initial guess being a good one, or the bracketed values actually bracketing the answer. It's more than possible to supply a terrible initial guess, especially for pipe diameter, and have the root solver fail outright. Adding some code to both check that guesses are within the domains of functions would be a start, e.g. catching attempts to take `log(0)` and returning `-Inf` or something to ensure that the root-finding algorithms respect function domains. This also presents an opportunity to generate better default values, programatically, prior to solving. As opposed to me just picking reasonable numbers off the top of my head and having everything work out because I'm lucky.

Relatedly there is a lot of room to fiddle around with which root finding algorithm is employed, this intersects with a design decision I made early on to use the library `Unitful` to track units throughout and ensure unit consistency. Unfortunately it does not play nicely with all root finding algorithms, for example the Brent method, as implemented in the `Roots` library, threw a bunch of arcane errors related to unit consistency when I tried using it. To get around this one could add a step to strip out units before solving and then add them back in after.
