---
title: "Air Dispersion Model - Worst Case Meterological Conditions"
toc: true
toc_label: "contents"
toc_sticky: true
comments: true
categories:
  - examples
tags:
  - air dispersion modeling
  - atmospheric stability
excerpt: This example considers the impact of different Pasquill stability classes and windspeeds on a gaussian dispersion model of a plume from an elevated stack.
---

# Air Dispersion Model - Worst Case Meterological Conditions

In a [previous post](https://selenized.github.io/2020/12/05/gaussian_dispersion_example/) I modeled an example of a plume from an elevated stack. In that example I assumed very stable conditions and a low windspeed -- pasquill stability class F and a windspeed of 1.5m/s -- as the worst case. This was an error!

For neutrally buoyant releases at or near ground-level that is a common "worst case", for example when considering the potential impact due to a vapour cloud explosion. But for elevated stacks, releasing a buoyant plume, a class D stability with a moderate windspeed is often recommended. I thought it would be interesting to explore how the maximum concentration at the point of interest -- an elevated work area downwind of the stack -- varies with stability class and windspeed.

I am not going to repeat all of the assumptions and working out of the previous notebook, the important results are in the source code for this post. I have also re-defined some of the functions to be a little more re-useable and to represent other stability cases not covered in the original notebook. See the jupyter notebook [here](https://github.com/selenized/selenized.github.io/tree/main/_notebooks).




## Pasquill Stability

As a refresher, Pasquill stability classes are a qualitative way of describing the atmospheric stability -- the tendency of the atmosphere to resist or enhance vertical motion. Stability is itself related to the temperature gradient with height, wind speed, and various other things. For a simple model such as this the key model parameters are tabulated with respect to the Pasquill stability, which is why it is relevant to this discussion.

#### Pasquill Stability Classes

| Stability Class | Description                | 
|-----------------|----------------------------| 
| A               | Extremely unstable         | 
| B               | Unstable                   | 
| C               | Slightly unstable          | 
| D               | Neutral                    | 
| E               | Slightly stable            | 
| F               | Stable to extremely stable | 

In general the more stable the class the less dispersion, and thus the higher the concentration within the plume. Which is why class F is typically used for a ground level, neutrally buoyant, cloud. However plume rise is also a function of stability and, in general, more stable plumes rise without as much dispersion and thus the *ground level* concentration is lower than if the plume dispersed more. Furthermore the plume rise is a function of windspeed, the greater the windspeed the less the plume rises before leveling off and, again, the greater the ground level concentration.

We can visualize the relationship between windspeed, stability class, and the concentration at the point of interest with the following plot which includes plume rise relationships for unstable, neutral, and stable atmospheres.



![svg](/images/worst_case_weather_files/output_7_0.svg)



We note that, as we expect, lower stability (e.g. A or B) corresponds to a higher groundlevel concentration at low windspeed, but at high windspeed higher stability leads to a greater groundlevel concentration.

At first blush it would appear that class F is still the worst case, however this plot naively assumes atmospheric stability is unrelated to windspeed. This is not true and roughly speaking the stability transitions towards classes C and D as the windspeed increases.


#### Pasquill Stability vs Incoming Solar Radiation

| Windspeed (m/s) | Strong | Moderate | Slight | 
|:---------------:|:------:|:--------:|:------:| 
| < 2             | A      |  A - B   | B      | 
| 2 - 3           | A - B  | B        | C      | 
| 3 - 5           | B      | B - C    | C      | 
| 5 - 6           | C      | C - D    | D      | 
| > 6             | C      | D        | D      | 


#### Pasquill Stability vs Nighttime Cloud Cover

| Windspeed (m/s) | > 4/8 cloud | < 3/8 cloud | 
|:---------------:|:-----------:|:-----------:| 
| < 2             |             |             | 
| 2 - 3           | E           | F           | 
| 3 - 5           | D           | E           | 
| 5 - 6           | D           | D           | 
| > 6             | D           | D           | 

To represent this crudely, the plots can be chopped off at the windspeed limits from the tables above. So, for example, the class F plot would end at 3m/s.

The plots also present an obvious way of finding the worst case for a particular scenario: find the extremal point for the worst stability class. This is found rather simply by setting the derivative to zero. Something that would be complicated to do analytically but is quite straight forward when using the `ForwardDiff` library for automatic differentiation.


```julia
# ForwardDiff doesn't play nicely with unitful
# so values have to be stripped of units first

Fb′ = ustrip(Fb)
x₁′ = ustrip(x₁)
h₁′ = ustrip(h₁)
Q′  = ustrip(Q)
hₛ′ = ustrip(hₛ)

∂Cᵤ(u) = ForwardDiff.derivative(u -> C(u, x=x₁′, 
                                          y=0.0, 
                                          z=h₁′, 
                                          Q=Q′, 
                                          h=hₛ′, 
                                          Δh=(x,u) -> Δhᵣ(x, u, Fb=Fb′, stable=false), 
                                          σy=σy("D"), 
                                          σz=σz("D")), float(u))
```



```julia
# Find the point where ∂C/∂u = 0
# Initial guess of 25 just by eye-ball

u_worst = find_zero(∂Cᵤ, 25)

C_worst = C(x₁, 0.0u"m", h₁, 
            u=u_worst*1u"m/s", 
            Q=Q, 
            h=hₛ, 
            Δh=(x,u) -> Δhᵣ(x, u, Fb=Fb, stable=false), 
            σy=σy("D"), 
            σz=σz("D"))

uconvert(u"mg/m^3", C_worst)
```


    0.2753490886768969 mg m^-3


![svg](/images/worst_case_weather_files/output_11_0.svg)



We find that the worst case is indeed class D but with quite a high windspeed, ~26.4m/s or 95kph, which would be considered a 10 on the [Beaufort scale](https://www.canada.ca/en/environment-climate-change/services/general-marine-weather-information/understanding-forecasts/beaufort-wind-scale-table.html) with trees being uprooted and considerable structural damage. It's unlikely that workers would still be on the platform and it may not even be the case that the scaffolding would still be standing!

Regardless we can look at the contour plots at the work platform elevation and vertically, along the centerline. 

*Note* the colours are scaled to 4mg/m^3, one tenth the occupational limit of 40mg/m^3.

![svg](/images/worst_case_weather_files/output_13_0.svg)



Another interesting impact of elevated releases like this is that the worst concentration for an observer on the ground is often a significant distance *downwind* of the stack. Because the plume must disperse downwards.

Below is a contour plot showing the downwind concentration at a 2m elevation -- the height of a reasonably tall person. Note the scale is set to even lower concentrations. The maximum of the colour bar is 1000x lower than the occupational limit.


![svg](/images/worst_case_weather_files/output_15_0.svg)



## Multiple Concentrations

Previously, in [the discussion of the occupational exposure limit](https://selenized.github.io/2020/12/05/gaussian_dispersion_example/#the-concentration-of-interest) I noted that, in general, one would have to account for the impact of multiple substances in the flue gas, though in that particular example I was only modeling carbon monoxide and I just moved on. I think I left the impression that one would have to model each substance separately and, at least with this simple gaussian dispersion model, that is very much not the case.

Consider for some substance *i* being released with in-stack concentration $C_{s,i}$, we can define a dimensionless "dilution" $\chi$ as

$$ \chi \left( x, y, z \right) = { C_i \left( x, y, z \right) \over C_{s,i} } $$

Assuming the in-stack concentration to be simply the mass emission rate of *i* divided by the volumetric flow-rate[^state]

$$ C_{s,i} = { Q_i \over V_s^o } $$

and recalling the concentration function from the gaussian dispersion model 

$$ C_i \left( x, y, z \right) = {Q_i \over 2 \pi u \sigma_{ye} \sigma_{ze} } f \left( x, y, z \right) $$

where $ f \left( x, y, z \right) $ is the products of the exponentials, and is a function of x, y, z only. Putting all that together we get an expression for the dilution that does not depend upon the substance being released

$$ \chi \left( x, y, z \right) = {V_s^o \over 2 \pi u \sigma_{ye} \sigma_{ze} } f \left( x, y, z \right) $$

If you have already done the modeling for a particular substance then calculate $\chi$ for the points of interest by dividing the concentrations by the in stack concentration, otherwise substitute the formula for $\chi$ given above for the concentration and model that instead.

Then, when evaluating multiple substances, the test[^oel]

$$ \sum_i {C_i \left( x, y, z \right) \over T_i } \lt 1 $$

becomes

$$ \chi \left( x, y, z \right) \cdot \sum_i {C_{s,i} \over T_i } \lt 1 $$

where $T_i$ is the relevant occupational exposure limit.

To recap, instead of calculating the concentrations at the points of interest using a gaussian dispersion model multiple times, calculate a dimensionless *dilution* at the points of interest and apply that to the in stack concentrations of all of the substances of interest. Then combine those as per the relevant rules for occupational hygiene.

Below are a series of contour plots showing the dilution $\chi$, where colours are are from 0-5% -- i.e. the concentration within the yellow region is &ge; 5% the in-stack concentration. 

*Note:* This is backwards to the usual way of defining dilution, where a $\chi$ of 5% would be a 95% dilution.

[^state]: at standard state, because the concentrations given for the occupational exposure limits are given in terms of a volume at standard state. This is also a potential error in the original model as it does not correct the concentrations back to standard state, nor does it really track temperature to make that even possible, especially near the stack.

[^oel]: from [CCOHS](https://www.ccohs.ca/oshanswers/hsprograms/occ_hygiene/occ_exposure_limits.html)


![svg](/images/worst_case_weather_files/output_17_0.svg)



## Thoughts on Code and Reusability

A simple way of taking the work from a previous notebook and adding to it is just to import the notebook. This loads the results into the current notebook including any function definitions and such.

For example

```julia

using NBInclude
@nbinclude("2020-12-05-gaussian_dispersion_example.ipynb")

```

I didn't do that here for two reasons:
1. I want these notebooks to be independent and stand on their own
1. I didn't write the previous notebook in a very extendable or reusable way

The second point is worth going into if one wants to build a library of worked out, generic, models as notebooks. This way an engineer can import previously defined models as needed for a particular analysis while also keeping the documentation for the models *in the model*. In the previous notebook I left most things in the global namespace and defined functions that used those global variables. Which is fine for that particular notebook but it means that the current workspace gets very cluttered when importing things and also those functions are not very re-usable as they were defined for a very particular example.

It's better, I think, to write functions that use keyword arguments for any important parameters that can then be passed as needed, instead of defining those parameters in the global namespace. Unless they are truly constants, like *g* the acceleration due to gravity or *R* the universal gas constant.

Another point is on the use of the library `Unitful`. It is convenient, and a good check, to have units propagate through calculations, however `Unitful` does not play nicely with all libraries. This is especially the case with `Plots` but it can also be real hassle to use with correlations that have lots of parameters. I think this is a good opportunity to take advantage of julia's multiple dispatch.

For example, suppose a correlation of the form $ f \left( x \right) = a \cdot x^b $, this can be written in julia very simply (supposing *a* and *b* are parameters)

```julia

f(x; a, b) = a*x^b

```

But if we pass *x* with some units and don't pass the matching units with the parameter *a* this will throw an error. We could tediously work out the units for each set of parameters *a* and *b* to make the units cancel out properly, or we could use multiple dispatch to manage this for us

```julia

function f(x::Quantity; a, b)::Quantity
    x′ = ustrip(u"expected input unit", x)
    return f(x′, a=a, b=b)*1u"correlation output unit"
end

```

Where we convert the input to the expected units, whatever they may be, evaluate the function in a unitless way, then tack on the expected output units at the end. Now when we use *f(x)* in contexts without units, for example when plotting *f(x)*, it works as expected and if we pass a value of *x* with units attached we get the unit conversion/checking that we want from `Unitful`.

---

