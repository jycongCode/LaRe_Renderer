1. Change the absorption parameters to zero in scenes/volpath_test/volpath_test1.xml. What do you
see? Why?
A: The object becomes brighter, and there is no difference in brightness between the edge and the center of the object. 
This is because when the absorption parameter is 0, the external medium does not absorb energy, and there is no scattering. 
Since the target object is self - luminous, it becomes brighter. 
Also, the transmittance is independent of t because the absorption parameter is 0, causing the difference in brightness to disappear. 

2. In the homework, we assume the volume being not emissive. If you were tasked to modify the pseudo
code above to add volume emission, how would you do it? Briefly describe your approach.
A:We assume that the volume has emittance and is homogeneous. 
Then, according to the radiative transfer equation (RTE), we need to add a term, which is volume_emit * t, after transmittance * Le.

3. In the derivation above, how did we get from p(t) ∝ exp(−σtt) to p(t) = σt exp(σtt)?
A: Because we can set p(t) = C*exp(−σtt), when t -> infinity, we have 1/σt = 1/C, therefore C = σt, then p(t) = σt*exp(−σtt)

4. How was Equation (13) incorporated into the pseudo code above? Why is it done this way?
A: Equation (13) is incorporated into the pseudo code by calculating the trans_pdf and transmittance directly, and if the
hit target is a light source, multiply the emission. It's done this way because for the part inside object, we don't consider
scattering(It's going to be occluded any way), then we only consider the later part of (9), which can be calculated in closed form.

Play with the parameters σs and σa, how do they affect the final image? Why? (we assume monochro-
matic volumes, so don’t set the parameters to be different for each channel.)
A: sigma_a affect's the overall brightness, while sigma_s mainly affects the area inside sphere more significant than the surrounding. It's because sigma_s determines the possiblity of light being scattered in the medium, less scattering makes the light source clearer than the surrounding, while sigma_a determines how much energy is aborped which affects the overall brightness.

5. Change the phase function from isotropic (the default one) to Henyey-Greenstein by uncommenting
the phase function specification in the scene scenes/volpath_test/volpath_test2.xml. Play with the
parameter g (the valid range is (−1, 1)). What does g mean? How does the g parameter change the
appearance? Why?

A: The letter "g" represents the average cosine value of the scattering direction. 
When g is around 0.5, it is the brightest, and there is also light around the object. 
However, when it approaches 1 and -1, there is no light around the object. This is because when it approaches -1 and 1, the scattering direction is concentrated completely backward and completely forward, resulting in a low probability of the light rays around the object sampling the light source.
Therefore, there is no "floodlight - like" effect. When it approaches 0.5, both light sources can be sampled, and since it is slightly forward - biased, the blue light avoids the pink light more obviously and has a larger proportion. 

6. Play with the parameters σs, σa of different volumes, and change max_depth, how do they affect the final image? How does increasing/decreasing σs and σa of medium1 and medium2 affect the appearance,
respectively? Why? Do different σs and σa values affect how high you should set max_depth?
A:  
m1: sigma_a: higher darker over all, lower brighter
    sigma_s: higher darker more noise, lower birghter light source,darker other parameters

m2: sigma_a: higher index-match sphere more solid, lower more transparent
    sigma_s: higher more solid and more noise the index-match sphere

when increasing max_depth, there is more noise, but it's a bit brighter overall.
m1 is the medium outside the shapes, higher sigma_a means more absorption, higher sigma_s means more scattering,
m2 is the medium inside the index-match sphere, therefore it mostly affect the index-match sphere itself.

When we increase max_depth gives higher variance, to make the noise less obvious, we should try to decrease sigma_s

7. Switch to the Henyey-Greenstein phase function again. How does changing the g parameter affect the
appearance? Why?
A: Reduce the g - value of m1 to a negative value to darken the overall picture, 
because the scattering direction is opposite to the incident direction. Whether it is increased or decreased, 
the result is the same. The light source becomes brighter, but the floodlight effect weakens. For the g - value of m2, 
it mainly affects the index - match objects. The overall change is not significant. However, when g is negative, 
the outline of the light source at the rear becomes less obvious. 

8. Propose a phase function yourself (don’t have to describe the exact mathematical form). How would
you design the shape of the phase function? What parameter would you set to control it?
A: I designed a phase function that will have two control parameters. 
One is g, which is similar to the Henyey - Greenstein parameter and represents the average cosine value. 
The other parameter is h, which represents the cosine value of the angle between the target direction and the incident direction. 
That is, h determines the overall orientation of our lobe. 
Then g specifies whether the overall distribution is more towards the direction specified by h or its opposite direction (similar to adding an additional degree of freedom to the Henyey - Greenstein function).

9. When will next event estimation be more efficient than phase function sampling? In our test scenes,
which one is more efficient? Why?
A: When the light source is larger or closer to the object, phase function sampling is more efficient. 
Conversely, next event estimation is more efficient when the light source is smaller or farther away. 
In our test scenario *volpath_test4*, where the light source is small, next event estimation is more efficient. 
However, for the previous test cases with larger light sources, phase function sampling is more efficient.

10. In scenes/volpath_test/volpath_test4_2.xml, we render a scene with an object composed of dense
volume. How does it compare to rendering the object directly with a Lambertian material? Why are
they alike or different?
A:Compared to the Lambertian model, our scene surface appears rougher and contains more noise due to the higher overall sampling variance. 
However, the shadows are softer than those in the Lambertian model. 
This is because Lambertian shading only considers the surface, whereas volume rendering accounts for scattering within the volume, resulting in smoother shadows without distinct boundaries.

11.Jim Kajiya famously has predicted in 1991 that in 10 years, all rendering will be volume rendering.
What do you think that makes him think so? Why hasn’t it happened yet?
Kajiya's reasoning for this prediction was likely that volume rendering could encompass many phenomena that surface rendering cannot easily represent, such as smoke, clouds, and other participating media. 
By integrating these effects, rendering could achieve more realistic simulations of indirect lighting and interactions between different objects.  
However, the reason this transition has not fully happened yet is likely due to computational limitations—real-time volume rendering remains extremely demanding. 
Additionally, surface rendering has performed exceptionally well in practical applications, offering both high efficiency and visual fidelity, making it the preferred choice in most real-world scenarios.

12. Play with the index of refraction parameter of the dielectric interface in
scenes/volpath_test/volpath_test5_2.xml. How does that affect appearance? Why?
First, reducing the IOR causes the object to become a dense sphere at an IOR of 1, making the center invisible. At an IOR of 0.8, the center becomes brighter, the surrounding medium becomes visible, but the outermost area darkens. At an IOR of 0.5, the center becomes even brighter, while the dark outer region expands. This occurs because decreasing the IOR makes the refracted light bend more toward the center, resulting in a brighter center and a darker periphery.  

Increasing the IOR to 1.5 makes the center brighter while the surrounding medium darkens as a whole, without forming a dark outer ring. When raised to 2, the surrounding area becomes darker, and the center becomes even brighter, but there is still no black hollow region. This is because increasing the IOR causes light to refract more outward rather than creating a hollow effect. At the same time, internal light tends to reflect more within the object due to Fresnel’s effect, making the center appear even brighter.

13. In the scene scenes/volpath_test/vol_cbox_teapot.xml, we model the glass teapot as a transparent
glass with blue homogeneous medium inside. What is the difference in terms of appearance between
this approach and just making the color of the glass blue without any medium inside?

If the surface is directly set to blue, the interior of the teapot will appear brighter because there is no absorption and scattering from the medium. Additionally, the teapot's shadow will be less noticeable.

14. For heterogeneous volumes, what kind of distribution of the volume density makes the null scattering
efficient/inefficient? Can you think of a way to improve our current sampling scheme in the inefficient
case?
A: For heterogeneous volumes, if the distribution is more "homogeneous" with less significant variation in volume
density makes the null scattering more efficient, according to sigma_t/sigma_m, we should have higher possiblity of 
hitting the 'real' particles. Thus a way of improving the sampling scheme is to divide the volumes into smaller regions,
each with its own majorant, making the variation less significant in each local area, then apply certain searching algorithms
to search for local majorant efficiently.

15. How do we make the null-scattering work for emissive volumes? Briefly describe a solution.
A: Similar to null scattering, we can define a emissive function sigma_e(t), and find its majorant,
then sample free-sample the light ray using sigma_e/sigmat just like with sigma_t, then add the local radiance to the total
radiance,

16. Why is it important to have an unbiased solution for volume rendering? Would it be sensible to have
something that is biased but faster? How would you do it?
A: Unbiased solution gives more physically correct result, which can be important to various field such as medical simulations and 
CG industry, with biased solution, increasing the sample count may lead to incorrect images. For tradeoff between efficiency and accuracy,
it is something sensible to do in realtime application which puts responsive first. We can do it by decreaing sampling count and max_null_scattering_count,
which may lead to noisy results, but can be alleviated by using various denoising methods.