import marimo

__generated_with = "0.11.22"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Lightwave Explorer
        Nick Karpowicz, Max Planck Institute of Quantum Optics

        ![Lightwave Explorer icon](https://raw.githubusercontent.com/NickKarpowicz/LightwaveExplorer/refs/heads/master/Source/BuildResources/Icon.svg)

        Thanks for using LWE and reading the manual! This Marimo (Python) notebook contains background on the physical system that is actually being solved in the simulation, and instructions on how to interact with it in the software. There are some interactive plots that let you see how some of the internals work: if you want to see the code used to generate the plots, hit the "..." button in the corner and select "Show code".

        Please let me know if something you want to know about is missing - this documentation is still incomplete (since the thing its describing is still a moving target) but I'm happy to help out if you have questions.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Physical model""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        #### Basic equations: The unidirectional propagation equation
        I'm sure we've all derived the wave equation from Maxwell's equations before, so I'm not going to repeat that. Let's just talk about the equations I'm using and where they come from.

        In LWE, there are two modes of propagation. The most commonly-used is the nonlinear, unidirectional wave equation, which has a typical form like this:

        \begin{equation}
        \frac{\partial}{\partial z} E(\mathbf{x} ,\omega) = i(k + \frac{1}{2k}\nabla^2_\perp)E(\mathbf{x} ,\omega) + \frac{i\omega}{2\epsilon_0 \tilde{n}\left(\omega\right)c}\mathbf{P}^{\mathrm{NL}}\left(\mathbf{x},\omega\right)\tag{1}
        \end{equation}

        where propagation of the field is assumed to be primarily along the $z$ axis. $\nabla^2_\perp$ is the transverse Laplacian (_e.g._ $\frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}$ in Cartesian coordinates). $\omega$ is the angular frequency (this equation is in the frequency domain), $\mathbf{P}^{NL}$ is the nonlinear polarization, $\tilde{n}$ is the complex index of refraction, and the other constants have their usual definition. Bold $\mathbf{x}$ is shorthand for either $x$, $y$, and $z$ in Cartesian coordinates, or $\rho$, $\theta$, and $z$ in cylindrical coordinates.

        In deriving this equation, we make the assumption that $\left|k_z^2\right| \gg \left|\frac{\partial^2}{\partial z^2} \tilde{A}(z)\right|$, where $\tilde{A}$ is the complex amplitude of a given frequency component of the electric field. This means that we are assuming that the field evolves gradually on the spatial scale of the wavelength. This is a good assumption when the light is propagating in the presence of "normal" nonlinearities, where the nonlinear signal builds up slowly over time. It's not a good assumption near interfaces, where the reflection from the surface suddenly changes the evolution of the wave, and there is a counterpropagting wave. So, for thin crystals, or surface nonlinearities, we'll want to fully solve Maxwell's equations. You'll notice that the spatial form of this equation looks like the Helmholtz equation. The underlying assumptions are quite similar.

        If $\mathbf{P}^{NL} = 0$, we have linear propagation of light, and for example a Gaussian beam is a solution to the above equation. When it's nonzero is when things get interesting. We'll go over how that happens in a bit.

        #### Basic equations: Finite-difference time-domain (FDTD)

        The other propagation mode we can use is FDTD mode, where Maxwell's equations are solved on a grid. This mode doesn't use the paraxial or slowly-evolving wave equations used in the unidirectional propagation equation above. This gives us a few things: effects at interfaces are handled more correctly, and we can simulate etalon effects inside of crystals (where the light makes multiple reflections inside, resulting in a modulated transmission vs. frequency).

        For this mode, we don't need to derive much, we just start with two equations from Maxwell. Specifically,

        Faraday's law:

        \begin{equation}
        \nabla \times \mathbf{E}=-\frac{\partial \mathbf{B}}{\partial t}\tag{2}
        \end{equation}

        and Ampere's law:

        \begin{equation}
        \nabla \times \mathbf{B} = \mu_0 \left(\mathbf{J} + \epsilon_0\frac{\partial\mathbf{E}}{\partial t}\right)\tag{3}
        \end{equation}

        These are just two coupled first-order partial differential equations. We have both quantities $\mathbf{E}$ and $\mathbf{B}$ definied on a grid in space, and these quantities advance over time. If something interesting is happening, there will also be a material, whose response to the field is $\mathbf{J}$. The physical model in the basic version of LWE for calculating the current is in the form of a series of Lorentzian oscillators. Their basic equations of motion look like this when we're dealing with the linear response:
        \begin{equation}
        \frac{\partial \mathbf{J}}{\partial t} = k_L \mathbf{E} - \gamma \mathbf{J} - \omega_0^2 \mathbf{P}\tag{4}
        \end{equation}

        \begin{equation}
        \frac{\partial\mathbf{P}}{\partial t} = \mathbf{J}\tag{5}
        \end{equation}

        where $k_L = \frac{Ne^2}{m_e}$ is the strength of the response depending of the density of dipoles N and their charge and mass, $\gamma$ describes the relaxation rate of the oscillator and $\omega_0$ is its resonance frequency. So now we have two more coupled, first-order partial differential equations. Since the first one contains $\mathbf{E}$ as a driving term, inside of the material, all of our equations are coupled together! When we start doing nonlinear optics, we just have to modify the driving term to make use of the nonlinear tensor (plus some additional bookkeeping if we're modeling the material's response with an instantaneous component).

        You can see that there are more quantities to take care of when we solve Maxwell's equations for our field. Even if we only want $\mathbf{E}$, we also have to keep track of $\mathbf{B}$, as well as $\mathbf{J}$ and $\mathbf{P}$ for *each* resonance that makes up the material response.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Linear crystal properties

        For every material we're going to deal with, we need to have a way of describing its refractive index. This is typically done with a parameterized equation known generally as a Sellmeier equation. There are many different forms of this that you'll find in literature. So far, I've only implemented two of them here:

        #### General fitting sellmeier equation
        \begin{equation}
        n^2 = a[0] + \frac{a[1] + a[2]\lambda^2}{\lambda^2 + a[3]} + \frac{a[4] + a[5]\lambda^2}{\lambda^2 + a[6]} + \frac{a[7] + a[8]\lambda^2}{\lambda^2 + a[9]} + \frac{a[10] + a[11]\lambda^2}{\lambda^2 + a[12]} + a[13] \lambda^2 + a[14] \lambda^4 + a[15] \lambda^6 + \frac{k a[16]}{a[17] - \omega^2 + i a[18]\omega} + \frac{k a[19]}{a[20] - \omega^2 + i a[21]\omega} \tag{6}
        \end{equation}

        There are quite a few terms. You're not expected to use all of them, they're just there to make it more likely that it can accommodate any equation you find in literature, containing common forms of various elements. And terms not present, you can just set the top value to zero. The constant $k = \frac{e^2}{\epsilon_0 m_e}$ has a value of about 3182.6 in SI units. The last two terms are complex valued; this allows you to add absorption, with a Lorentzian line shape.

        Let's look at a concrete example. The sliders below will let you change the values of a[3] and a[6], adjusting the positions of the UV and IR resonances.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    #Slider to control the UV resonance
    UV = mo.ui.slider(start=0.5, stop=2.0, step=0.1, value=1, label="UV factor")
    UV
    return (UV,)


@app.cell(hide_code=True)
def _(mo):
    #Slider to control the IR resonance
    IR = mo.ui.slider(start=0.5, stop=2.0, step=0.1, value=1, label="IR factor")
    IR
    return (IR,)


@app.cell(hide_code=True)
def _(IR, UV, mo):
    #importing some modules I'll want
    import LightwaveExplorer as lwe
    import io
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'Verdana', 'DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif']


    def showmo():
        """
        Helper function to plot as an svg and have it display in marimo in vector form
        """
        svg_buffer = io.StringIO()
        plt.savefig(svg_buffer, format='svg')
        svg_buffer.seek(0)
        svg_data = svg_buffer.getvalue()
        return mo.Html(svg_data)

    #first we'll make a wavelength grid to work with
    l = np.linspace(0.3,3,1024)

    #next we'll need Sellmeier coefficients, these are for barium fluoride, H. H. Li., J. Phys. Chem. Ref. Data 9, 161-289 (1980)
    a = np.array([1.33973,0,0.81070,-0.010130,0,0.19652,-892.22,0,4.52469,-2896.6,0,0,1,0,0,0,0,0,0,0,0,0])

    #we can get the refractive index for the wavelengths we put in the grid by calling the sellmeier() function
    #from the lightwaveExplorer module, with the equationType set to 0.
    n = lwe.sellmeier(l, a, 0)

    #let's make it so we can adjust the resonances of the oscillators and see how it affects the index

    fig,ax = plt.subplots(1,1, figsize=(5,4))
    a2 = np.array([1.33973,0,0.81070,-0.010130,0,0.19652,-892.22,0,4.52469,-2896.6,0,0,1,0,0,0,0,0,0,0,0,0])
    a2[3] *= UV.value
    a2[6] *= IR.value
    n2 = lwe.sellmeier(l, a2, 0)
    ax.plot(l,np.real(n),label="original",color="blue")
    ax.plot(l,np.real(n2), label = "modified", color = "magenta")
    ax.set_xlabel("Wavelength (microns)")
    ax.set_ylabel("n")
    ax.legend()
    showmo()
    return a, a2, ax, fig, io, l, lwe, n, n2, np, plt, rcParams, showmo


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        You can play with the position of the resonances and see how this affects the value of the refractive index, and its dispersion. Feel free to add sliders or attach them to other parameters to see what they do. Using the function above, you can use a fitting routine to get coefficients in the appropriate format if there isn't a 1-to-1 mapping from literature to the equation.

        Another equation that's implemented in the code is all Lorentzians:

        \begin{equation}
        n^2 = a[0] + \frac{k a[1]}{a[2] - \omega^2 + i a[3]\omega} + \frac{k a[4]}{a[5] - \omega^2 + i a[6]\omega} + ... \tag{7}
        \end{equation}

        all the way up to $a[21]$.

        This one is required for FDTD: it is simply the response of a collection of Lorentzian oscillators, and obeys causality. The previous equation (although more flexible) does not obey causality (there's no absorption corresponding to some of its resonances so by definition it doesn't obey the Kramers-Kronig relations) and doesn't have a purely time-domain equation of motion associated with it.

        The final equation that's implemented is a series of Gaussian resonances, whose associated contribution to the real part of the susceptibility is a Dawson function F (the Dawson function is the Hilbert transform of a Gaussian). This equation does in fact obey causality and could in principle be associated with a pure time-domain response suitable for FDTD, but is not yet implemented there because I haven't found an efficient way to do it (although I'm convinced that one exists). Gaussian functions (i.e. inhomogenously broadened lines) often exhibit better agreement with resonances inside of solids. This model looks like:

        \begin{equation}
        n^2 = a[0] + \sum_{i}{-\frac{1}{\sqrt{\pi}}a[3i]\left(F\left({\omega_s(i)}\right) - e^{-\omega_s^2(i)}\right)}\tag{8}
        \end{equation}

        \begin{equation}
        \omega_s = \frac{\omega - a[1 + 3i]}{\sqrt{2} a[2+3i]}\tag{9}
        \end{equation}

        which is to say, each resonance has a strength given by a parameter in the array, and a scaled frequency axis given by $\omega_s$, which is scaled by two additional parameters.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Nonlinear crystal properties
        -----

        There are currently three different types of nonlinear effects that can happen inside of the simulation: second order, such as second harmonic generation and difference frequency generation, third order, such as the Kerr effect, four-wave mixing, and third harmonic generation, and multiphoton absorption/nonlinear plasma interaction in a classical picture of free carriers.

        The nonlinear tensors are supplied to the program in the form of the contracted notation for $\chi^{(2)}$, as this is the most readily available in the literature, and in the full 81-component form of $\chi^{(3)}$. For $\chi^{(3)}$, a second mode is allowed, where only one value is supplied, $\chi^{(3)}_{1111}$, and the propagation will work under the assumtion of a centrosymmetric medium. This is often all the information available in literature for a given material, even if we know that the tensor has other independent elements. If your measurement is sensitive to these unknown elements, maybe with the right modelling and fitting, you can provide the first experimental values.

        ### Plasma formation
        In many crystals, the damage mechanism will be multiphoton absorption followed by free-carrier absorption. In other cases, these effects will limit the intensity inside of the medium, through intensity claming, or the plasma formation may participate in your nonlinear process. This phenomenon is complex (and interesting!) but a full quantum mechanical calculation is not included presently. Instead, there is a classical approximation that has empiracally been observed to be quite effective: the Drude model.

        This calculation is done in the time domain:

        \begin{equation}
        J^{\mathrm{Drude}}_i(t) = \frac{e^2}{\mu_e}\exp{\left(-\gamma t\right)}\int_{-\infty}^t dt' \exp{\left(\gamma t\right)} N(t) E_i(t)\tag{10}
        \end{equation}

        where $e$ is the electron charge, $\mu_e$ is the reduced effective mass, $\gamma$ is the momentum relaxation rate, and N(t) is the density of free carriers. $E_i$ indicates the field in the $i$-direction, where $E$ with no index indicates the total field magnitude, $E = \sqrt{\sum E_i^2}$


        To perform this calculation, we need the time-dependent carrier density $N$ and must conserve energy in the absorption process. Exciting the system and introducing a pair of carriers costs energy, in the form of $n$ photons. The field must then also see an absorptive nonlinearity, whose associated current is:

        \begin{equation}
        J^\mathrm{abs}_i(t) = \beta E^{2n-2}(t)E_i(t)\tag{11}
        \end{equation}

        where $\beta$ is a nonlinear absorption parameter, which will depend on the system and frequency of the excitation. From this rate of nonlinear absorption, the power transfer from the field to the system will be $\sum J^\mathrm{abs}_i(t)E_i(t)$.

        Accordingly, the number of carriers generated, assuming that they just go to the $\Gamma$ point of the crystal is

        \begin{equation}
        \frac{d N(t)}{d t} = \frac{2}{\Delta_g}\sum J^\mathrm{abs}_i(t)E_i(t)\tag{12}
        \end{equation}

        where $\Delta_g$ is the band gap. The factor of 2 is due to the creation of an electron-hole pair.


        #### Dispersion of the nonlinear coefficients
        In all cases, we take into account the dispersion of the nonlinear coefficients through Miller's rule, or an extrapolation of it to higher orders. Miller's rule is the idea that the nonlinear susceptibilities will scale with the linear ones, so in the case of $\chi^{(2)}$:

        \begin{equation}
        \chi^{(2)}_{ijk}\left(\omega_1+\delta_1,\omega_2+\delta_2;\omega_3+\delta_3\right) \approx \frac{\chi^{(1)}_i(\omega_1 + \delta_1)\chi^{(1)}_j(\omega_2+\delta_2)\chi^{(1)}_k(\omega_3+\delta_3)}{\chi^{(1)}_i(\omega_1)\chi^{(1)}_j(\omega_2)\chi^{(1)}_k(\omega_3)}\chi^{(2)}_{ijk}\left(\omega_1,\omega_2;\omega_3\right)\tag{13}
        \end{equation}

        Doing so allows us to approximate the dispersive behavior of the nonlinear coefficients to a good degree of accuracy, especially if the resonances are Lorentzian or far away from all frequencies. To extend this to higher order nonlinearities, the pattern is identical, but the numbered indicies extend up to higher numbers. 

        Miller's rule strictly speaking is the case up to 3. The higher the nonlinearity, the more likely one is to have components close to resonances, where things can go wrong: for example, infrared-active phonon modes will affect $\chi^{(3)}$, instead of the Raman-active ones. For calculating resonant behavior, a different calculation should be added to the nonlinear polarization tailored to the specific problem.

        As a result, the applied nonlinearities are not instantaneous: the response to a $\delta$-function excitation will not be another $\delta$-function. Practically, this helps avoid overestimates e.g. of infrared nonlinearities based on coefficients taken from visible measurements. It also means that when configuring a crystal, you need to supply the frequencies at which the measurement that provided the elements of the nonlinear tensor was performed. Without these, the response _will_ be treated as instantaneous.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        #### Crystal properties in FDTD mode
        In FDTD mode, the crystal response is contained in two places: a frequency-independent modification of the dielectric constant, e.g. $\epsilon_\infty$ obtained from the first element of the Sellmeier coefficient array, and the current $\mathbf{J}$. The former is straight-forward: in Ampere's law, we just replace $\epsilon_0$ with $\epsilon_0\epsilon_\infty$. The rest of the response, including all of the nonlinear optical responses, comes through $\mathbf{J}$, and requires more care.

        Ideally, the response we end up with should be consistent with the dispersion of the nonlinear effects obtained through Miller's rule above, such that the different modes of propagation have similar approximations to the local polarization response. We just have to translate the frequency-domain picture into a purely time-domain one such that we can solve it in FDTD.

        It turns out to be fairly straightforward. Applying Miller's rule is equivalent to using a scaled version of the nonlinear coefficient, and instead of mixing the fields $E_1(t)$ and $E_2(t)$, we mix the polarizations $P_1(t) = \chi^{(1)}(t) \ast E_1(t)$ and $P_2(t) = \chi^{(1)}(t) \ast E_2(t)$. Here I'm using the symbol $\ast$ to denote convolution. Satisfying the above formulation of Miller's rule requires us to calculate the nonlinear polarization as (for the example of a second order nonlinearity):

        \begin{equation}
        P_k^{(2)}(t) = \chi^{1}(t) \ast \chi^{(2)*}_{ijk}P_i(t)P_k(t)\tag{14}
        \end{equation}

        where $P_k$ is a given polarization component, and $\chi^{(2)*}_{ijk}$ is a given scaled component of the non-contracted nonlinear tensor. Essentially, we just have to do an additional convolution of the nonlinear response with the linear response function.

        First, we have to scale the nonlinear tensor. What I mean by that is that typically in literature, we'll have a single value of the tensor component, e.g. a number in pm/V, which was measured in a given way. For example, if it was measured using second harmonic generation of a 1064-nm laser to 532-nm, it corresponds specifically to the components of $\chi^{(2)}$(1064 nm, 1064 nm, 532 nm). Using the above equations formulating the response in terms of polarizations, we can extend this to an arbitrary combination of wavelengths, but we must first apply the denominator of the Miller's rule expression to get the scaled tensor components:

        \begin{equation}
        \chi^{(2)*} = \frac{\chi^{(2)}(\omega_1,\omega_2;\omega_3)}{\chi^{(1)}(\omega_1)\chi^{(1)}(\omega_2)\chi^{(1)}(\omega_3)}\tag{15}
        \end{equation}

        within the Miller's rule model, these scaled tensor components are frequency-independent. The frequency-dependence of the resulting polarization response comes through the application of the first-order response function to the polarizations serving as inputs to this equation, and the final convolution of the resulting nonlinear polarization it contains.

        We have in fact already set up the required numerical system for performing this convolution in the time domain -- it's just the oscillators whose equations of motion we are solving to obtain the current required for linear propagation. We just have to replace the electric field in the driving term with a more generalized force $\mathbf{F}$:

        \begin{equation}
        \frac{1}{e}\mathbf{F} = \mathbf{E} + \mathbf{P}^{N*}\tag{16}
        \end{equation}

        where $\mathbf{P}^{N*}$ is the sum of all of the nonlinear polarizations calculated using linear polarization as the input and the scaled nonlinear tensors.

        Simply adding this additional driving term to the equations of motion for the oscillators allows the fully-time domain calculation of the nonlinear polarization within the Miller's rule approximation, but is missing one component: the *instantaneous* response in some of the Sellmeier equations, i.e. when a[0] doesn't equal 1. This part of the linear response is not associated with an oscillator, so it has to be handled independently.

        Since in calculating the change in field through Ampere's law we need the associated *current* and not the polarization, it's a bit tricky!

        In addition to the nonlinear driver term $\mathbf{P}^{N*}$, we also need to find its time-derivative $\mathbf{J}^{N*}$ because the contribution of the instantaneous response to the nonlinear current will be $(\epsilon_\infty - 1)\mathbf{J}^{N*}$. Calculating this requires having access to the derivative of the polarization. The parts of it associated with the oscillators are easy, but due to the presence of an instantaneous term, that bit that we eliminated from explicit calculation by the substitution $\epsilon_0$ to $\epsilon_0\epsilon_\infty$, we need to obtain the time-derivative of the field so that we can calculate the full time derivative:

        \begin{equation}
        \frac{\partial\mathbf{P}}{\partial t} = (\epsilon_\infty - 1) \frac{\partial \mathbf{E}}{\partial t} + \mathbf{J^{(1)}}\tag{17}
        \end{equation}

        We can obtain $\frac{\partial \mathbf{E}}{\partial t}$ and $\mathbf{J^{(1)}}$ (the linear response dipole current) through the solution of Ampere's law that we've already done through the calculations for however!

        We then just need to calculate the contributions to $\mathbf{J}$ from the instantaneous responses using these values and the product rule for derivatives, e.g. for the example above:

        \begin{equation}
        J_k^{\infty} = \frac{\partial P_k^{(2)}(t)}{\partial t} = (\epsilon_\infty - 1) \chi^{(2)*}_{ijk}\left(\frac{\partial P_i(t)}{\partial t}P_k(t) + P_i(t)\frac{\partial P_k(t)}{\partial t}\right)\tag{18}
        \end{equation}

        Adding these current terms (remember we have to do this product rule for each combination of fields corresponding to a tensor component!) we have the complete nonlinear response. And we've done it purely in the time domain, which is required for FDTD.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Matrix math for rotating into and out of the crystal frame

        A confusing aspect of nonlinear optics is getting into and out of the crystal coordinate system, starting from the coordinate system of the beam.
        Here, I'm assuming that the beam is pointing in the z-direction, and we have a set of crystal angles, $\theta$ and $\phi$ which give the orientation of the crystal surface relative to the internal coordinates of the crystal.

        We need to move from the beam coordinates $(x, y, z)$ to the crystal coordinates $(x', y', z')$. Essentially the crystal angles give the orientation of the beam's propagation vector ($\hat{z}$) relative to the $\hat{z}'$ axis of the crystal. This takes the form of several rotations, and we have to do them in the correct order and direction.

        We'll start out in a simple case, a uniaxial crystal, and add complications as necessary for dealing with the more complex crystal systems.

        Say we have our fields in the beam coordinate system. We want to move to the crystal coordinates in order to perform a calculation. Because we are already rotated relative to the principle axes of the crystal, we effectively have to undo the rotation we've accumulated. The correct order for this is to first rotate by $-\theta$ around the y'-axis, and then by $-\phi$ around the z'-axis. The rotation matrices take the usual form:

        $$
        \begin{equation}
        R_y\left(\theta\right) = \left[\begin{matrix}\cos{\left(\theta \right)} & 0 & \sin{\left(\theta \right)}\\0 & 1 & 0\\- \sin{\left(\theta \right)} & 0 & \cos{\left(\theta \right)}\end{matrix}\right]\tag{19}
        \end{equation}
        $$

        $$
        \begin{equation}
        R_z\left(\phi\right) = \left[\begin{matrix}\cos{\left(\phi \right)} & - \sin{\left(\phi \right)} & 0\\\sin{\left(\phi \right)} & \cos{\left(\phi \right)} & 0\\0 & 0 & 1\end{matrix}\right]\tag{20}
        \end{equation}
        $$

        The combined forward matrix $R_f\left(\theta,\phi\right) = R_z\left(-\phi\right) R_y\left(-\theta\right)$ is then:

        $$
        \begin{equation}
        R_f\left(\theta,\phi\right) = \left[\begin{matrix}\cos{\left(\phi \right)} \cos{\left(\theta \right)} & \sin{\left(\phi \right)} & - \sin{\left(\theta \right)} \cos{\left(\phi \right)}\\- \sin{\left(\phi \right)} \cos{\left(\theta \right)} & \cos{\left(\phi \right)} & \sin{\left(\phi \right)} \sin{\left(\theta \right)}\\\sin{\left(\theta \right)} & 0 & \cos{\left(\theta \right)}\end{matrix}\right]\tag{21}
        \end{equation}
        $$

        To get back out of the crystal system, we have the opposite sign of the angles and apply the rotations in reverse order, so the backwards matrix $R_b\left(\theta,\phi\right) =  R_y\left(\theta\right) R_z\left(\phi\right)$

        $$
        \begin{equation}
        R_b\left(\theta,\phi\right) = \left[\begin{matrix}\cos{\left(\phi \right)} \cos{\left(\theta \right)} & - \sin{\left(\phi \right)} \cos{\left(\theta \right)} & \sin{\left(\theta \right)}\\\sin{\left(\phi \right)} & \cos{\left(\phi \right)} & 0\\- \sin{\left(\theta \right)} \cos{\left(\phi \right)} & \sin{\left(\phi \right)} \sin{\left(\theta \right)} & \cos{\left(\theta \right)}\end{matrix}\right]\tag{22}
        \end{equation}
        $$

        In the current version of the code, this pair of rotations is done at each propagation step, rotating from the beam frame to the crystal frame, calculating the nonlinear polarization along the principle axes and then rotating the polarization back into the beam frame. You may notice that this could be optimized by pre-calculating an effective tensor in the beam frame. This optimization has been left for a later date, so that the first correctness checks are done in an easier-to-diagnose form, where the a priori known tensors are in the inner loop. Since this rotation is much faster than the Fourier transforms needed to obtain the polarization in any case, it is not expected to be a big deal.

        These same matrices naturally give us expressions for the refractive indices for the two principle polarization components $E_x, E_y$. In the typical manner of the refractive index ellipsoid, we have to deal with the diagonal matrix N in the cystal frame:

        $$
        \begin{equation}
        N = \left[\begin{matrix}\frac{1}{n_{x}^{2}} & 0 & 0\\0 & \frac{1}{n_{y}^{2}} & 0\\0 & 0 & \frac{1}{n_{z}^{2}}\end{matrix}\right]\tag{23}
        \end{equation}
        $$

        The effective refractive indices for the beam's principle axes are then found along the diagonal of the matrix $N_e = R_b\left(\theta,\phi\right) N R_f\left(\theta,\phi\right)$

        $$
        \begin{equation}
        N_e = \left[\begin{matrix}\frac{\sin^{2}{\left(\theta \right)}}{n_{z}^{2}} + \frac{\sin^{2}{\left(\phi \right)} \cos^{2}{\left(\theta \right)}}{n_{y}^{2}} + \frac{\cos^{2}{\left(\phi \right)} \cos^{2}{\left(\theta \right)}}{n_{x}^{2}} & - \frac{\sin{\left(\phi \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)}}{n_{y}^{2}} + \frac{\sin{\left(\phi \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)}}{n_{x}^{2}} & \frac{\sin{\left(\theta \right)} \cos{\left(\theta \right)}}{n_{z}^{2}} - \frac{\sin^{2}{\left(\phi \right)} \sin{\left(\theta \right)} \cos{\left(\theta \right)}}{n_{y}^{2}} - \frac{\sin{\left(\theta \right)} \cos^{2}{\left(\phi \right)} \cos{\left(\theta \right)}}{n_{x}^{2}}\\- \frac{\sin{\left(\phi \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)}}{n_{y}^{2}} + \frac{\sin{\left(\phi \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)}}{n_{x}^{2}} & \frac{\cos^{2}{\left(\phi \right)}}{n_{y}^{2}} + \frac{\sin^{2}{\left(\phi \right)}}{n_{x}^{2}} & \frac{\sin{\left(\phi \right)} \sin{\left(\theta \right)} \cos{\left(\phi \right)}}{n_{y}^{2}} - \frac{\sin{\left(\phi \right)} \sin{\left(\theta \right)} \cos{\left(\phi \right)}}{n_{x}^{2}}\\\frac{\sin{\left(\theta \right)} \cos{\left(\theta \right)}}{n_{z}^{2}} - \frac{\sin^{2}{\left(\phi \right)} \sin{\left(\theta \right)} \cos{\left(\theta \right)}}{n_{y}^{2}} - \frac{\sin{\left(\theta \right)} \cos^{2}{\left(\phi \right)} \cos{\left(\theta \right)}}{n_{x}^{2}} & \frac{\sin{\left(\phi \right)} \sin{\left(\theta \right)} \cos{\left(\phi \right)}}{n_{y}^{2}} - \frac{\sin{\left(\phi \right)} \sin{\left(\theta \right)} \cos{\left(\phi \right)}}{n_{x}^{2}} & \frac{\cos^{2}{\left(\theta \right)}}{n_{z}^{2}} + \frac{\sin^{2}{\left(\phi \right)} \sin^{2}{\left(\theta \right)}}{n_{y}^{2}} + \frac{\sin^{2}{\left(\theta \right)} \cos^{2}{\left(\phi \right)}}{n_{x}^{2}}\end{matrix}\right]\tag{24}
        \end{equation}
        $$

        If $n_x = n_y$ (i.e. we have a uniaxial crystal), the off-diagonal values at the 0,1 and 1,0 indices of the matrix are zero, and we can simply read off the first two diagonal elements to get the extraordinary and ordinary refractive index values for the beam. 


        Now, let's look at more complicated crystals, biaxial ones which still have orthogonal axes. This includes orthorhombic systems. $n_x$, $n_y$, and $n_z$ are independent, so we have to perform an additional rotation around the beam's z-axis to orient the beam x and y to the refractive index ellipsoid. 

        First, we redefine our forward and backward rotation matrices to incorporate the rotation around the beam's z-axis, $\delta$:

        $$
        \begin{equation}
        R_f\left(\theta,\phi,\delta\right) = \left[\begin{matrix}- \sin{\left(\delta \right)} \sin{\left(\phi \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)} & \sin{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)} + \sin{\left(\phi \right)} \cos{\left(\delta \right)} & - \sin{\left(\theta \right)} \cos{\left(\phi \right)}\\- \sin{\left(\delta \right)} \cos{\left(\phi \right)} - \sin{\left(\phi \right)} \cos{\left(\delta \right)} \cos{\left(\theta \right)} & - \sin{\left(\delta \right)} \sin{\left(\phi \right)} \cos{\left(\theta \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)} & \sin{\left(\phi \right)} \sin{\left(\theta \right)}\\\sin{\left(\theta \right)} \cos{\left(\delta \right)} & \sin{\left(\delta \right)} \sin{\left(\theta \right)} & \cos{\left(\theta \right)}\end{matrix}\right]\tag{25}
        \end{equation}
        $$

        $$
        \begin{equation}
        R_b\left(\theta,\phi,\delta\right) = \left[\begin{matrix}- \sin{\left(\delta \right)} \sin{\left(\phi \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)} & - \sin{\left(\delta \right)} \cos{\left(\phi \right)} - \sin{\left(\phi \right)} \cos{\left(\delta \right)} \cos{\left(\theta \right)} & \sin{\left(\theta \right)} \cos{\left(\delta \right)}\\\sin{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)} + \sin{\left(\phi \right)} \cos{\left(\delta \right)} & - \sin{\left(\delta \right)} \sin{\left(\phi \right)} \cos{\left(\theta \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)} & \sin{\left(\delta \right)} \sin{\left(\theta \right)}\\- \sin{\left(\theta \right)} \cos{\left(\phi \right)} & \sin{\left(\phi \right)} \sin{\left(\theta \right)} & \cos{\left(\theta \right)}\end{matrix}\right]\tag{26}
        \end{equation}
        $$

        And the (big) matrix giving the refractive index values is $N_e = R_b\left(\theta,\phi,\delta\right) N R_f\left(\theta,\phi,\delta\right)$:

        $$
        \begin{equation}
        \left[\begin{matrix}\frac{\sin^{2}{\left(\theta \right)} \cos^{2}{\left(\delta \right)}}{n_{z}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \cos{\left(\phi \right)} - \sin{\left(\phi \right)} \cos{\left(\delta \right)} \cos{\left(\theta \right)}\right)^{2}}{n_{y}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \sin{\left(\phi \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)}\right)^{2}}{n_{x}^{2}} & \frac{\sin{\left(\delta \right)} \sin^{2}{\left(\theta \right)} \cos{\left(\delta \right)}}{n_{z}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \cos{\left(\phi \right)} - \sin{\left(\phi \right)} \cos{\left(\delta \right)} \cos{\left(\theta \right)}\right) \left(- \sin{\left(\delta \right)} \sin{\left(\phi \right)} \cos{\left(\theta \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)}\right)}{n_{y}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \sin{\left(\phi \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)}\right) \left(\sin{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)} + \sin{\left(\phi \right)} \cos{\left(\delta \right)}\right)}{n_{x}^{2}} & \frac{\sin{\left(\theta \right)} \cos{\left(\delta \right)} \cos{\left(\theta \right)}}{n_{z}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \cos{\left(\phi \right)} - \sin{\left(\phi \right)} \cos{\left(\delta \right)} \cos{\left(\theta \right)}\right) \sin{\left(\phi \right)} \sin{\left(\theta \right)}}{n_{y}^{2}} - \frac{\left(- \sin{\left(\delta \right)} \sin{\left(\phi \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)}\right) \sin{\left(\theta \right)} \cos{\left(\phi \right)}}{n_{x}^{2}}\\\frac{\sin{\left(\delta \right)} \sin^{2}{\left(\theta \right)} \cos{\left(\delta \right)}}{n_{z}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \cos{\left(\phi \right)} - \sin{\left(\phi \right)} \cos{\left(\delta \right)} \cos{\left(\theta \right)}\right) \left(- \sin{\left(\delta \right)} \sin{\left(\phi \right)} \cos{\left(\theta \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)}\right)}{n_{y}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \sin{\left(\phi \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)}\right) \left(\sin{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)} + \sin{\left(\phi \right)} \cos{\left(\delta \right)}\right)}{n_{x}^{2}} & \frac{\sin^{2}{\left(\delta \right)} \sin^{2}{\left(\theta \right)}}{n_{z}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \sin{\left(\phi \right)} \cos{\left(\theta \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)}\right)^{2}}{n_{y}^{2}} + \frac{\left(\sin{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)} + \sin{\left(\phi \right)} \cos{\left(\delta \right)}\right)^{2}}{n_{x}^{2}} & \frac{\sin{\left(\delta \right)} \sin{\left(\theta \right)} \cos{\left(\theta \right)}}{n_{z}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \sin{\left(\phi \right)} \cos{\left(\theta \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)}\right) \sin{\left(\phi \right)} \sin{\left(\theta \right)}}{n_{y}^{2}} - \frac{\left(\sin{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)} + \sin{\left(\phi \right)} \cos{\left(\delta \right)}\right) \sin{\left(\theta \right)} \cos{\left(\phi \right)}}{n_{x}^{2}}\\\frac{\sin{\left(\theta \right)} \cos{\left(\delta \right)} \cos{\left(\theta \right)}}{n_{z}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \cos{\left(\phi \right)} - \sin{\left(\phi \right)} \cos{\left(\delta \right)} \cos{\left(\theta \right)}\right) \sin{\left(\phi \right)} \sin{\left(\theta \right)}}{n_{y}^{2}} - \frac{\left(- \sin{\left(\delta \right)} \sin{\left(\phi \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)}\right) \sin{\left(\theta \right)} \cos{\left(\phi \right)}}{n_{x}^{2}} & \frac{\sin{\left(\delta \right)} \sin{\left(\theta \right)} \cos{\left(\theta \right)}}{n_{z}^{2}} + \frac{\left(- \sin{\left(\delta \right)} \sin{\left(\phi \right)} \cos{\left(\theta \right)} + \cos{\left(\delta \right)} \cos{\left(\phi \right)}\right) \sin{\left(\phi \right)} \sin{\left(\theta \right)}}{n_{y}^{2}} - \frac{\left(\sin{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)} + \sin{\left(\phi \right)} \cos{\left(\delta \right)}\right) \sin{\left(\theta \right)} \cos{\left(\phi \right)}}{n_{x}^{2}} & \frac{\cos^{2}{\left(\theta \right)}}{n_{z}^{2}} + \frac{\sin^{2}{\left(\phi \right)} \sin^{2}{\left(\theta \right)}}{n_{y}^{2}} + \frac{\sin^{2}{\left(\theta \right)} \cos^{2}{\left(\phi \right)}}{n_{x}^{2}}\end{matrix}\right]\tag{27}
        \end{equation}
        $$

        The value of $\delta$ that sets the off-diagonals at 0,1 and 1,0 to zero is:

        $$
        \begin{equation}
        \tan{2\delta} = \frac{\sin{2\phi}\cos{\theta}}{\Omega\sin^2\theta - \cos^2\phi\cos^2\theta + \sin^2\phi}\tag{28}
        \end{equation}
        $$

        where

        $$
        \begin{equation}
        \Omega = \frac{\frac{1}{n^2_y} - \frac{1}{n^2_z}}{\frac{1}{n^2_x} - \frac{1}{n^2_y}}.\tag{29}
        \end{equation}
        $$

        We can now write down expressions for the effective refractive indices for the principle field axes in the beam coordinates:

        $$
        \begin{equation}
        n_x^e = \sqrt{\frac{n_{x}^{2} n_{y}^{2} n_{z}^{2}}{n_{x}^{2} n_{y}^{2} \sin^{2}{\left(\theta \right)} \cos^{2}{\left(\delta \right)} + n_{x}^{2} n_{z}^{2} \left(\sin{\left(\delta \right)} \cos{\left(\phi \right)} + \sin{\left(\phi \right)} \cos{\left(\delta \right)} \cos{\left(\theta \right)}\right)^{2} + n_{y}^{2} n_{z}^{2} \left(\sin{\left(\delta \right)} \sin{\left(\phi \right)} - \cos{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)}\right)^{2}}}\tag{30}
        \end{equation}
        $$

        $$
        \begin{equation}
        n_y^e = \sqrt{\frac{n_{x}^{2} n_{y}^{2} n_{z}^{2}}{n_{x}^{2} n_{y}^{2} \sin^{2}{\left(\delta \right)} \sin^{2}{\left(\theta \right)} + n_{x}^{2} n_{z}^{2} \left(\sin{\left(\delta \right)} \sin{\left(\phi \right)} \cos{\left(\theta \right)} - \cos{\left(\delta \right)} \cos{\left(\phi \right)}\right)^{2} + n_{y}^{2} n_{z}^{2} \left(\sin{\left(\delta \right)} \cos{\left(\phi \right)} \cos{\left(\theta \right)} + \sin{\left(\phi \right)} \cos{\left(\delta \right)}\right)^{2}}}\tag{31}
        \end{equation}
        $$

        Note that for the case of a uniaxial crystal, where $n_x = n_y$ (so $\delta = 0$), this reduces significantly. $n_x^e$ simplifies to the expected behavior of the extraordinary axis, and $n_y^e$ follows that of the ordinary axis (and no longer depends on angle, just becoming the un-rotated value of $n_y$).

        The angle $\delta$ may depend on frequency. Accordingly, the optical axes (the orientation of the refractive index ellipse) may be different for different wavelength components of the light. This has a significant implication for the linear and nonlinear propagation of the light field: there has to be an additional frequency-dependent rotation of the field vector (by $\delta$ around the z-axis) after Fourier transformation to enter the basis of plane waves where the propagator is diagonalized. This system no longer corresponds to the lab frame and has to be exited before results are interpreted.

        In LWE, the propagation in a biaxial crystal will be done using the crystal system so that the effect of the wavelength-dependence of the angle $\delta$ will be present, but the result will always be in the lab frame. This means that the portions of the propagation which take place in the frequency basis will be rotated by $\delta$, but the result will not be. Accordinly, the polarization angle of the light might not automatically align with the optical axes (since there may not be a well-defined optical axis for a broadband field). A pair of helper functions, described in the section on sequences, called rotateIntoBiaxial() and rotateFromBiaxial() make this easier.

        In the FDTD mode, the oscillators are always in the crystal frame. The electric field vector is rotated into this frame, using $R_f\left(\theta,\phi\right)$, and then the results are rotated back to the field frame using $R_b\left(\theta,\phi\right)$.

        Now, let's add the final level of complication. When working with the refractive index tensor above, we assumed it was diagonal, but this is not always the case. Let's look at the triclinic crystal.

        In a triclinic system, there are 9 different nonzero values in the linear susceptibility tensor:

        $$
        \begin{equation}
        \chi^{(1)} = \left[\begin{matrix}xx & xy & xz \\ yx & yy & yz \\ zx & zy & zz\end{matrix}\right]\tag{32}
        \end{equation}
        $$

        In a monoclinic system, it's a bit simpler, with only one pair of off-diagonal elements:

        $$
        \begin{equation}
        \chi^{(1)} = \left[\begin{matrix}xx & 0 & xz \\ 0 & yy & 0 \\ zx & 0 & zz\end{matrix}\right]\tag{33}
        \end{equation}
        $$

        In both of these systems, we might not be able to use the assumption we made above that it's possible in principle to get the refractive index matrix into a diagonal form. In the worst case, we would have to solve this system not with two electric field components coupled only by the nonlinearity like we were assuming for the UPPE propagation model, but instead treat even the linear propagation as a coupled differential equation.

        However, if we can assume that the tensors are invariant with respect to permutation of the axes (e.g. $zx=xz$), which should hold far away from resonances, they'll be symmetric, and we can define yet another rotation to diagonalize these matrices. Since the off-diagonal components will still have a frequency-dependent dispersion, the angle will still be wavelength dependent, however. Within the approximation of a fixed angle, we can rotate the tensor into an effective orthorhombic system. It is only to this level of approximation that monoclinic and triclinic systems are currently treated within LWE, although a further addition is planned to fully accommodate them.

        Note that when you do this, you will also have to rotate any nonlinear tensors associated with the system, since they are usually given with their indices defined by the crystal axes, not the optical ones (and the symmetry that you can use to check their validity is also in the crystal frame of reference). There are functions in the LWE python module to apply this to the chi(2)/d tensors, called chi2axisSwap() and chi2rotate() to make this easier.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Pulse parameters

        We also need to describe the initial condition of the electric field, as that is what's going to evolve over the course of the simulation. Currently, you can easily add 2 fields, and go up to an arbitrary number with a bit of effort. These fields maybe be input as a number of pulse parameters on the interface, or you may load the results of a FROG or EOS measurement. To use a FROG, you need to provide a .Speck.dat file from the FROG software. I haven't decided on the EOS format yet, but it should likely be a two-column ascii file, first column time, second, field.

        If you use a synthesized field, it will be defined in the frequency domain, in the following form:

        \begin{equation}
        \tilde{E}(\omega) = A_o e^{-(\omega-\omega_o)^{N_s}/\Delta_\omega^{N_s} - i\phi(\omega)}\tag{34}
        \end{equation}

        where $\omega_o$ is the central frequency of the pulse, $\Delta_\omega$ is the bandwidth, $N_s$ is the order of the supergaussian, and $\phi(\omega)$ is the spectral phase, which is:

        \begin{equation}
        \phi(\omega) = \phi_{\mathrm{ce}} + \omega\tau + \frac{1}{2}\phi_2(\omega-\omega_o)^2 + \frac{1}{6}\phi_3(\omega-\omega_o)^3 + \phi_\mathrm{m}(\omega)\tag{34}
        \end{equation}

        These are:
        - $\phi_\mathrm{ce}$: carrier-envelope phase
        - $\tau$: delay (the actual value contains an offset of half the time grid length, so that $\tau = 0$ means the pulse is centered)
        - $\phi_2$: Group delay dispersion (GDD), a.k.a. chirp
        - $\phi_3$: Third order dispersion (TOD)
        - $\phi_\mathrm{m}(\omega)$ = $\frac{n_\mathrm{m}(\omega) - n_\mathrm{m}({\omega_o})}{c}\omega L$, the phase acquired via linear propagation through a user-selected medium of length L.

        The interactive plot below will let you see what these parameters do to the field as defined in the code, and tell you the full-width-at-half-maximum intensity (FWHM) of the resulting pulse.
        """
    )
    return


@app.cell
def _(mo):
    f0 = mo.ui.slider(start=50, stop=1000, step=10, value=400, label="Frequency")
    f0
    return (f0,)


@app.cell
def _(mo):
    bandwidth = mo.ui.slider(start=2, stop=300, value=80, label="Bandwidth")
    bandwidth
    return (bandwidth,)


@app.cell
def _(mo):
    Ns = mo.ui.slider(start=2, stop=32, step=2, value=2, label="SG order")
    Ns
    return (Ns,)


@app.cell
def _(mo):
    cep = mo.ui.slider(start=0.0, stop=2.0, step=0.1, value=0, label="CEP/pi")
    cep
    return (cep,)


@app.cell
def _(mo):
    tau = mo.ui.slider(start=-20, stop=20, step=1, value=0, label="Delay")
    tau
    return (tau,)


@app.cell
def _(mo):
    phi2 = mo.ui.slider(start=-100, stop=100, step=1, value=0, label="GDD")
    phi2
    return (phi2,)


@app.cell
def _(mo):
    phi3 = mo.ui.slider(start=-300, stop=300, step=1, value=0, label="TOD")
    phi3
    return (phi3,)


@app.cell
def _(Ns, bandwidth, cep, f0, lwe, np, phi2, phi3, plt, showmo, tau):
    dt = 0.25e-15
    Tgrid = 80e-15
    Nt = int(2*Tgrid/dt)
    t = np.linspace(-Tgrid,Tgrid,Nt)
    w = 2*np.pi*np.fft.fftfreq(Nt,dt)
    ws = (w-2*np.pi*f0.value*1e12)
    phi = np.pi*cep.value + (-Tgrid + tau.value*1e-15)*w + 0.5*phi2.value*1e-30*ws**2 + (1.0/6)*phi3.value*1e-45*ws**3
    Ew = np.exp(-ws**Ns.value/(2*np.pi*bandwidth.value*1e12)**Ns.value -1.0j * phi)
    Ew[w<0]=0
    Et = lwe.norma(np.fft.ifft(Ew))
    fwhm_value = 1e15*lwe.fwhm(t,np.abs(Et)**2)
    plt.plot(1e15*t,np.real(Et))
    plt.xlabel("Time (fs)")
    plt.title(f"FWHM duration: {fwhm_value: .2f} fs")
    showmo()
    return Et, Ew, Nt, Tgrid, dt, fwhm_value, phi, t, w, ws


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Numerical model""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        #### Propagation scheme, convergence, stability

        The numerical propagation is performed using the fourth-order Runge-Kutta method. This method is not unconditionally stable, which means that for certain values of the grid parameters, the results will diverge, and lead to NaN values. The simulation will quit if this happens, and tell you. In order to avoid this, you have to set the grid spacing parameters correctly, this means the time step, spatial step and propagation step.

        Although the specifics of how these schemes become unstable are well understood (see the Courant-Friedrichs-Lewy Condition), in practice it's usually best to progressively make your steps smaller until you get a result, then to continue to make them smaller until nothing changes. 

        A step of the RK4 loop of the calculation procedes as follows:
        1. Calculate things that require the time-domain field:
            - Nonlinear polarization
            - Plasma currents
            - Radial Laplacian (if in cylindrical coordinates)
        2. Use results of 1, and prepared linear propagation frequency-domain factors, to provide an estimate of the field after 1/2 a propagation step
        3. Repeat 1 with the new field
        4. Update the 1/2 step estimate with the new polarization
        5. Repeat 1 with the new field
        6. Estimate the final step field with the results
        7. Repeat 1 with the new field
        8. Step to the next propagation step using the steps calculated from 1-7

        You can see that this requires the polarizations to be evaluated 4 times in order to advance one propagation step. This is more effort per step than using the "split step method" you'll hear old textbooks talk about, but ends up being much more accurate in practical calculations.

        When the pulse is propagated in Cartesian coordinates, the field is propagated in a basis of states comprising the plane-wave solutions to Fourier optics as the beam propagates along the z-direction, whose projection into real space (needed to obtain the nonlinear polarization, for example) is given by the standard Fourier optics propagator:

        $$
        \begin{equation}
        E(\mathbf{x} + z\hat{z},\omega) = e^{i k_z z} E(\mathbf{x},\omega),\tag{35}
        \end{equation}
        $$

        where

        $$
        \begin{equation}
        k_z = \sqrt{|k|^2 - k_x^2 - k_y^2}\tag{36}
        \end{equation}
        $$

        The nonlinear wave equation is significantly simplified in this basis, with all of the linear propagation terms, including diffraction, becoming zero. The linear propagation is handled exactly by the evolution of the basis states since we use the Fourier solution, such that only the nonlinear buildup has to be handled numerically.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        #### Cylindrical symmetry tricks

        When the beam can be assumed to be symmetric around the z-axis, we can't use the Fourier basis, but we can still simplify things (and speed up the simulation) a lot, since the grid can be reduced from 3 dimensions, $(x, y, t)$, to 2, $(\rho, t)$. 

        The only thing that makes life more difficult after this coordinate transformation is the transverse Laplacian from eqn. (1). In cylindrical coordinates with radial symmetry, it becomes

        \begin{equation}
        \nabla^2_\perp = \frac{\partial^2}{\partial \rho^2} + \frac{1}{\rho} \frac{\partial}{\partial \rho} \tag{37}
        \end{equation}

        Compared to a Cartesian coordinate, we have an extra term to deal with. And unfortunately, unlike Cartesian coordinates, that term is not a diagonal operator in the Fourier domain.

        One approach is to make a Hankel transform instead of a Fourier transform, in which case the operator will be diagonal. Numerically this can be implemented in a fairly efficient way leveraging FFTs, but only if the data is arranged on a logarithmic grid.

        Another way is to use an evenly spaced grid, but explicitly calculate the $\frac{1}{\rho}$ term in real space. This way one uses FFTs and the same grid as the calculations in Cartesian coordinates, but with the expense of having to do a Fourier transform at each propagation step. Since we have to do that anyway for the inclusion of nonlinearities, it's not a big deal in this case.

        The grid used in this case extends to the negative side of the origin, so that the FFT doesn't see a discontinuity at $\rho = 0$, covering the range $-\rho_{\textrm{max}}$ to $+\rho_{\textrm{max}}$. This would seem to be a waste of memory since the positive and negative sides of the origin are identical by symmetry.

        If we use a grid that is truly symmetric, that is indeed the case: we would be solving everything at $\rho = ... -3, -2, -1, 0, +1, +2, +3...$ etc.

        However, what if we shift the placement of our grid points by $\frac{1}{4}$? Now we have the fields at $\rho = ... -2\frac{3}{4}, -1\frac{3}{4}, -\frac{3}{4}, +\frac{1}{4}, +1\frac{1}{4}, +2\frac{1}{4}, +3\frac{1}{4}$

        This way, the positive and negative sides of the grid have different information: the previously-useless negative-$\rho$ values contain the midpoints between the points on the positive $\rho$ grid. The spacing is still 1 however, meaning that the grid doesn't include high values of the transverse $k$, which often bring numerical instability, in order to get doubled resolution in space. We just have to "interlace" the grid using our knowledge of the symmetry of the system. 

        We make use of this in two places in the cylindrical symmetry propagation mode: to eliminate _aliasing_ in the calculation of the nonlinear polarization, and to improve the accuracy of the finite-difference calculation of the radial derivative included in that $\frac{1}{\rho} \frac{\partial}{\partial \rho}$ term of the Laplacian.

        This results in overall easier convergence of calculations in cylindrical symmetry mode, which is why it's a good place to start when approaching a new system.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Interface

        Here is a view of the user interface, with significant parts labeled in the large colored text

        <img src="https://raw.githubusercontent.com/NickKarpowicz/LightwaveExplorer/refs/heads/master/Documentation/Images/LabeledInterface.svg">

        Let's go through these components one-by-one.

        ### Pulse parameters
        These you find in two columns on the left. Here you can put in all the parameters described in the Pulse Parameters section, plus a few more.
        1. Pulse energy: total energy contained in the pulse; the intensity will be determined by this value
        2. Frequency: $f_0$ as defined above, in THz
        3. Bandwidth: $\Delta_\omega$ as defined above, in THz
        4. SG order: $N_s$ as defined above, giving the order of the supergaussian to use. 2 is a normal gaussian. Only even values are allowed.
        5. CEP/$\pi$: $\phi_0$ as defined above, the carrier-envelope phase. In units of radians/$\pi$, e.g. a value of 1 will invert the field.
        6. Delay: $tau$, in fs. Allows you to move your pulse around the grid - negative values move it earlier in time.
        7. GDD: $\phi_2$, the group delay dispersion, in fs $^2$
        7. GDD: $\phi_3$, the third order dispersion, in fs $^3$
        8. Phase material: You can select from the loaded material in the crystal database, shown in the info box, a material whose dispersion you would like to apply to the pulse. For example, if you leave it as 2, it will apply the phase of passing through BaF2.
        9. Thickness: Thickness in microns of the phase material to apply. If you set it to 1000, the pulse will be stretched (or compressed) as though it passed through 1 mm of barium fluoride. Negative values would allow you to precompensate for the dispersion.
        10. Beamwaist: The Gaussian waist $w_o$. This is not necessarily the size of the beam at the start of the simulation, but what size the focus will have when the beam reaches in via linear propagation.
        11. x offset: amount by which to shift the beam on the spatial grid. Will always be 0 when radial symmetry is applied.
        12. z offset: amount by which to shift the beam relative to the propagation coordinate. 0 means that the beam is focused at the beginning of the simulation. If it is negative, it will be diverging as the simulation starts, and if it's positive, it will be converging.
        13. NC angle: Noncollinear angle of the beam. If it is 0, the beam will be pointing along the z-axis.
        14. Polarization: Angle in degrees of the field direction relative to the x-y plane. 0 means the field is along the x axis, 90 means it is along the y axis.
        15. Circularity: degree of ellipticity of the field. 0 means linearly polarized, 1 means circular.

        #### Grid parameters
        These describe the propagation space of the pulse.

        1. Material index: first column is the number of the material (from the list of loaded database entries displayed in the info box when the program starts). Second column is an alternate material to be applied in certain special modes.
        2. Theta, phi: the crystal propagation angles in degrees. They select the phase matching and nonlinear coefficients of the interaction.
        3. NL absorption: The parameters describing the nonlinear absorption; $\beta$ (in units determined by the order of the interaction), and $\Delta_g$ (in eV).
        4. Drude model parameters: the momentum relaxation time $\gamma$ and the reduced effective mass $\mu$
        5. Grid width, dx: the size and resolution of the spatial grid, respectively. In microns.
        6. Time span, dt: The length and resolution of the temporal grid, in fs.
        7. Length, dz: How far to propagate, and the spatial step size of the simulation, in microns.
        8. Propagation mode: select between full 3D, radial symmetry mode, or 2D cartesian coordinates.
        9. Batch mode: This is a bit complicated, but it allows you to run a sequence of simulations, sweeping one variable.
            - first, select from the drop down the parameter you wish to change. The pulse and material properties you set with floating point values are all adjustable here.
            - in "batch end" put the value to which it should sweep, in the same units as the initial value you gave. For example, if "Delay" on the left-hand-side of the interface is set to "-50", and you put "50" in the "Batch end" box, the series will go from -50 to 50.
            - in Batch steps, put the number of steps in the sweep. If it is 1, only the initial value will be used. If it is zero or less, it will set itself to 1.
        10. Batch mode 2 works the same way, but uses the values in the second column for batch end and batch steps. If you use both, you can do a 2D parameter scan.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Using sequences

        Below the batch mode, you'll find a text box labelel "Crystal sequence". This is where you can program in a more complex series of events to take place rather than simply propagate through one crystal.

        You can program a variety of things to happen, and utilize for loops and variables to simulated complicated systems.

        The syntax is similar to c/c++ - you can write out a series of functions, applying white space as you see fit. Semicolons are not necessary after function calls, but you can add them if you want. 

        Here's an example of a sequence making use of loops in order so simulate a periodically-poled lithium niobate (PPLN) crystal, with two different periods along the propagation direction:

        ```
            set(1,14.6)
            for(0.5*500/v01,0){
                nonlinear(d,d,d,v1,d)
                rotate(180)
            }
            for(17,0){
               nonlinear(d,d,d,9.25,d)
               rotate(180)
            }
        ```

        Let's go through what the different terms mean.

        First we have the set() function. This is defined further below, but in short it allows you to assign a value to a variable. You can access up to 100 variables (named v00 to v99) if you want to - which one you access is the first parameter in the function call. The second parameter of the call is what you set it to. So in this case, we are setting the variable v01 to 14.6.

        Next, we have a for loop. The syntax of a for loop takes two parameters: first, the number of iterations, second, the variable in which to store the counter. This for loop will loop over 0.5*500/v1 iterations, which means that it is reading the variable set in the previous line (v1) to determine how many times to run. This is so the loop can run over a fixed PPLN thickness instead of a fixed number of periods. The variable v0 will contain the index of the for loop (i.e. how many times it has run, starting from 0). This isn't used here, but could for example be used to make a chirped PPLN.

        Next, we have the nonlinear() function. This runs the nonlinear propagation, as given in the function definition below. You might notice that a lot of the parameters just have 'd' instead of a number. This tells the interpreter just to take the default value of the parameter, i.e. to take the value given by the interface. Any value with d will be able to be scanned in a batch mode or changed by the fitting routine, as it will be tied to the currently active value.

        Following that is the rotate() function, which rotates the field relative to the crystal. From the field's point-of-view, it's like the orientation of the crystal was flipped upside down, allowing quasi-phasematching to work.

        The for loop is then closed by curly brackets, and another one starts, doing the same thing, over 17 iterations with half-period of 9.25 microns.

        #### Variables

        There are three types of variables that can be accessed: iXX, vXX and d.

        ##### iXX
        These are interface variables, which correspond to the parameters that are set by the interface. They are called by the letter i followed by a number, where the number is the identifier of the associated batch mode. That is, the number you can see in the pull-down menu when selecting which parameter to scan for a batch. So, for example, i01 will give you the pulse energy in Joules. Setting the length of a crystal to i34^2 will cause its length to be the square of the length (in microns) as entered on the interface. 

        Currently the values of the iXX variables cannot be changed from within the sequence. They might do so if I can convince myself it won't lead to any disasters.

        ##### vXX
        These are user-controlled variables. You can set them to whatever you want with the set() function: set(2,i34) will set the value of v02 to the length of the crystal. You can use them later to determine input parameters of functions. They start with their values initialized to zero, so that's what they'll be if you use them without setting them.

        ##### d
        The keyword d will set the input parameter of a function to the corresponding value of the interface. If you call a nonlinear crystal with nonlinear(d,d,d,d,d), its material index, $\theta$, $\phi$, thickness and propagation step will just be whatever is on the interface, after being potentially modified by a batch scan or fitting command.

        ### Functions
        A number of functions are available for sequences. They are defined here. In order to initialize the electric field grid, a subset of functions must be called at least once before other functions which only modify the field grid are called. See the definition to know which type you're dealing with.

        #### plasma(materialIndex, crystalTheta, crystalPhi, nonlinearAbsorptionStrength, bandGapElectronVolts, drudeGamma, effectiveMass, crystalThickness, propagationStep)
        This function runs a nonlinear propagation through a medium determined by materialIndex.

        A call to plasma() will initialize the electric field grid if it hasn't been done so already.

        #### nonlinear(materialIndex, crystalTheta, crystalPhi, crystalThickness, propagationStep)
        This function will run a propagation through a medium determined by materialIndex, but with no nonlinear absorption or plasma formation.

        A call to nonlinear() will initialize the electric field grid if it hasn't been done so already.

        #### default()
        This function will run a nonlinear propagation through the medium selected on the interface, taking all values from the interface. It is equivalent to calling plasma(d,d,d,d,d,d,d,d,d).

        A call to default() will initialize the electric field grid if it hasn't been done so already.

        #### init()
        This function will create the electric fields in vacuum. This is useful if you want to simulate how the fields refract as they enter the crystal of interest. It's equivalent to calling nonlinear(0,0,0,0,1).

        A call to init() will initialize the electric field grid if it hasn't been done so already. That's probably the only reason you'd want to call it.

        #### addPulse(energy, frequency, bandwidth, sgOrder, cep, delay, gdd, tod, phaseMaterial, phaseThickness, beamwaist, x0, y0, z0, beamAngle, beamAngleY, polarization, circularity, materialIndex, theta, phi)
        This will add a new pulse to the grid, using all of the parameters of the interface. Note that you also have to specify the material in which the pulse is being created. 

        You can repeat this as many times as you want, to create many pulses.

        Note that compared to the interface, you have an extra beam angle and beam position. These are only relevant to 3D calculations and gives the angle and offset in the y-direction. (you can actually specify them on the interface, too, by writing both angles in the text box separated by a semicolon in the box for x0 and the NC angle, e.g. 10;20 will give an x-direction value of 10 and a y-direction value of 20).

        A call to addPulse() will initialize the electric field grid if it hasn't been done so already.

        #### linear(materialIndex, crystalTheta, crystalPhi, crystalThickness, propagationStep)
        This function will propagate *linearly* through the medium determined by materialIndex. The parameter propagationStep is only necessary when in cylindrical symmetry mode. This is because the propagation is carried out analytically in 3D and 2D Cartesian modes since the linear propagator is diagonal in those bases. In cylindrical symmetry mode, it is handled as a numerical propagation with all nonlinearities disabled, which is why specifying a step size is necessary. This is also why this function is significantly faster in 3D and 2D Cartesian coordinates.

        This function will *not* initialize the grid - it is necessary to precede it with one that does.

        #### rotate(angle)
        This will rotate the electric field around its propagation axis by an angle of its input parameter in degrees. For example, rotate(90) will change s-polarized light to p-polarized light.

        This function will *not* initialize the grid - it is necessary to precede it with one that does.

        #### rotateIntoBiaxial(materialIndex, crystalTheta, crystalPhi, frequency)
        This will rotate the field by the angle $\delta$ defined above in the context of the refractive index ellipse for biaxial crystals, aligning the polarization axes of the beam to the refractive index ellipse for the input material, crystal angles, and frequency. The similar function rotateFromBiaxial() will rotate in the opposite direction, and should be applied after propagation.

        This function will *not* initialize the grid - it is necessary to precede it with one that does.

        #### fdtd(timeFactor, dz, frontSpacer, backSpacer, obervationPlane, waitTime, preserveNearField, loadedGridIndex)
        This will use finite difference time domain (FDTD) mode with more fine-grained control over the system. The input parameters are as follows:
         - timeFactor: the difference in time-step between the saved data the one put in as dt on the interface, and the one used in the FDTD simulation. For example, if the value is 5, and dt on the interface is 125 attoseconds, the time step in the simulation will be 25 as. Default value: 5
         - dz: the spatial step along the z-direction to use. Default: value entered on the interface.
         - frontSpacer: if propagating to a crystal, place this distance in vacuum (meters) in front of the surface, so that the light doesn't immediately hit the crystal at the edge of the grid. Default value 1e-6 m
         - backSpacer: place this distance (m) after the crystal before the edge of the grid (e.g. extend the grid in the z-direction by this amount). Default: 1e-6 m.
         - observationPlane: The location (meters) within the numerical grid to record the field as a function of time. The field is recorded in a plane whose surface normal is the z-axis. 0 is the front surface of the grid.
         - waitTime: time (s) to let the simulation progress before recording the field in the observation plane. This is useful to not have a long period of zero field before the field from the start of the grid reaches the observation plane. Default 0.
         - preserveNearField: 1 or 0. 1 will keep the field as it was during the simulation, 0 will filter out spatial frequencies which will not propagate to the far-field.
         - loadedGridIndex: if you have loaded a material map from a file, this should be the index of that file in the loaded optics database. If you set this to d, it will instead construct a flat crystal in the center of the grid, sandwiched between vacuum spacers given by frontSpacer and backSpacer.


        #### fdtdReflection(timeFactor, dz)
        This is a streamlined version of the fdtd function for calculating the reflection from the surface of a crystal. The input parameters are just the first two of the full fdtd() function.


        #### fresnelLoss(materialIndex1, crystalTheta1, crystalPhi1, materialIndex2, crystalTheta2, crystalPhi2)
        This function will apply the Fresnel loss associated with transmission through the interface between the media determined by materialIndex1 and materialIndex2. The first material is the material from which the beam is incident, and the second is the new medium the beam is entering.

        This function is currently not working! It will do nothing if you apply it.

        This function will *not* initialize the grid - it is necessary to precede it with one that does.


        #### sphericalMirror(radiusOfCurvature)
        This will apply a radial phase associated with reflection from a spherical mirror. The radiusOfCurvature parameter sets the radius (and thus focal length x 2) in meters.

        This function will *not* initialize the grid - it is necessary to precede it with one that does.

        #### parabolicMirror(focalLength)
        This will apply a spatial phase associated with reflection (on-axis) from a parabolic mirror. focalLength is the focal length in meters.

        This function will *not* initialize the grid - it is necessary to precede it with one that does.

        #### aperture(diameter, activationParameter)
        This function applies an aperture to the field. The diameter parameter gives the diameter of the opening in meters. The activationParameter determines how "soft" the aperture is. For an activation function $p$, it takes the form:

        $1 - \frac{1}{1 + \exp{\left(-p(r-r_o)\right)}}$

        This function will *not* initialize the grid - it is necessary to precede it with one that does.

        ### farFieldAperture(openingAngle, activationParameter, alphaX, alphaY)
        This applies an aperture in the far field. Conceptually this is somewhat different than the elements as it doesn't correspond to an actual thin optic you can place in your beam. It essentially provides the effect of allowing the field to propagate into the far field, and selecting a conical section of the light determined by the angles openingAngle (e.g. the apex angle of the cone) and the two angles alphaX and alphaY which give its orientation relative to the propagation axis. In other words, you place an aperture far away from the starting position of the field, propagate through it, and then image the field back to the starting point.

        The activation parameter sets the softness of the edges of the cone using the same sigmoid function defined for the aperture() function.

        *NOTE: This function only gives quantitatively correct results in 3D mode at the moment - for the other modes I will fix it, but right now consider its effects to be qualitative-only!*

        This function will *not* initialize the grid - it is necessary to precede it with one that does.

        ### filter(f0, bandwidth, order, inBandAmplitude, outOfBandAmplitude)
        This applies a frequency filter to the field, of the form:

        $F(f) = \mathrm{outOfBandAmplitude} + \mathrm{inBandAmplitude} \times e^{-(f-f_0)^{\mathrm{order}}/(2\times \mathrm{bandwidth}^\mathrm{order})}$

        e.g. with order 4, when outOfBandAmplitude is 0 and inBandAmplitude is 1, it becomes a supergaussian bandpass filter around f0. 

        If outOfBandAmplitude is 1 and inBandAmplitude is -1, it is a bandstop filter.

        If outOfBandAmplitude is 1 and inBandAmplitude is 2, it will apply a gain of 3 in the pass band.

        f0 and bandwidth are in units of THz.

        Filters applied by this function do not modify the phase of the field, and thus do not obey causality. A real filter will distort the spectral phase, so think of any subsequent interactions as taking place in the best-case scenario concerning pulse distortions.

        This function will *not* initialize the grid - it is necessary to precede it with one that does.

        #### set(variableIndex, value)
        This will set the variable determined by variableIndex (e.g. variableIndex=0 will set the value of v0) to have a value of "value".

        This function will *not* initialize the grid, but it doesn't care if the grid is already initialized or not, and will work either way.

        ### energy(variableIndex, type)
        This function will put the energy contained in the field into the variable determined by variableIndex. The type says whether you want just the energy on the x-axis (0), or on the y-axis(1), or all of the energy (2).

        #### for(N, variableIndex)
        This allows you to make a for loop with N iterations. The current iterator value is stored in the variable corresponding to variableIndex

        for(10,1){...} will do whatever ... is 10 times, while counting from 0 to 9 in v1.


        #### <> (Angle brackets)
        These don't act as a function, but they allow you to disable a segment of code while leaving it intact, or add comments. They are similar to /* */ in c++. for example:
        ~~~
            default()
            <The functions below will not run>
            <rotate(90)
            sphericalMirror(1)>
        ~~~
        This code will just run default(). (Note, if you only see "default()" in the above block, there's something funky with how angle brackets are being displayed...) If the braces around the rotation and spherical mirror functions are removed, they will also be applied.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        #### Fitting mode

        The fitting mode allows for automatic adjustment of the light/crystal parameters in order to either maximize the power in a given spectral band, or to match the spectrum of an input reference in that band.

        A fitting command has the form:
        ```
        [beginning of ROI, Hz] [end of ROI, Hz] [Max iterations];
        [Parameter index1] [min1] [max1]; [parameter index2] [min2] [max2]...
        ```
        For example, if you set it to "Maximize x" in the pull-down menu, the command
        ```
        100e12 200e12 1000;
        29 0 90;
        30 0 180;
        ```
        will do the following: do at most 1000 calculations, trying to maximize the intensity in the range from 100 THz to 200 THz, by adjusting the crystal angles $\theta$ and $\phi$, over the full range of values. The units should be in the same units as the rest of the interface.

        Fitting is done with a global optimization, and does not depend on the initial values of the search parameters, only the ranges. Smaller ranges will lead to faster convergence.

        The numbers corresponding to the parameters are the same as those assigned to the different batch modes: you can see them all in the batch pull-down menus.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        #### File formats

        When you run a simulation, with the name simulationName in the path textbox (above the plots), it will automatically save the results when it finishes, creating several files:
        1. simulationName.txt : a text file containing all of the parameters of the simulation. You can load these files using the load button, and all the parameters will be entered into the interface automatically. If you rename it to "DefaultValues.ini" and put it in the same folder as the LightwaveExplorer.exe file, those values will be automatically put in when the program starts.
        2. simulationName_Ext.dat : a binary file (double precision) containing the entire electric field grid.
        3. simulationName_spectrum.dat: a binary file containing the integrated spectra at the end of the simulation
        4. simulationName.py : a python code snippet that will load the binary files and give them the proper shape

        Note that the best way to load these results is not to use the .py or .m files, they are just there as a backup in case you can't open them any other way. The best way is to use the python lightwaveExplorer module:
        ```
        import lightwaveExplorer as lwe
        s = lwe.load("simulationName.txt")
        ```
        will create the object s, whose members contain all of the results and parameters, as well as vectors giving the axes of the grid, spectra, and any batch scans you performed.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        #### Python interface

        There is a python module, LightwaveExplorer, that makes interacting with the results easier. You can install it with
        ```
        pip install LightwaveExplorer
        ```
        It has functions to load the results of a simulation (.load()), as well as different processing methods, and support functions.

        The most important function to know is "load()". You can load a saved result by providing it the name of the associated text file. For example "LWE_result = lwe.load("MyResult.txt")" will give you an object containing the results and input parameters of the simulation.

        The important member values to know for interpreting data are:

         - .Ext_x and .Ext_y: the electric field as a function of space and time, for the x and y polarizations
         - .spectrum_x and .spectrum_y: The polarization-resolved integrated spectra of the results
         - .spectrumTotal: The spectrum of all the resulting light (just the sum of the last two)
         - .timeVector: an array of the time values (e.g. time-axis of Ext_x)
         - .spaceVector: array of the x, y,or r spatial values
         - .frequencyVectorSpectrum: an array giving the frequencies corresponding to the spectra arrays
         - .batchVector: an array of the values of the scanned parameter in the batch mode. That is, if you scanned the z-offset of a beam from -2 to 2 in 5 steps, it will contain [-2,-1,0,1,2].

        If you're using an editor that gives you autocomplete, that's probably the easiest way to find the rest; you also have all of your input parameters, like .pulseEnergy1 and so on.

        The LightwaveExplorer module (which I load as "import LightwavExplorer as lwe" above), also contains some helper functions that I'll describe in more detail soon. These include:
         - sellmeier(): a copy of the sellmeier equations used inside the simulation, so you can give it the same array of parameters and plot the corresponding refractive index (useful when making new entries in the crystal database)
         - getSellmeierFromRII(): you can give this function a link to the YAML file labeled "Full database record" on the site refractiveindex.info, and it will return an array of the associated Sellmeier coefficients in an array that can be used by lwe.sellmeier() or pasted into the CrystalDatabase.txt.
         - getTabulatedDataFromRII(): If the data on refractiveindex.info contains tabulated n and k values, this will retrieve those in a numpy array (first column wavelength, then n, then k)
         - sellmeierFit(): from a tabulated set of n and k vs. wavelength, fit the parameters of a Sellmeier equation to them. I will post an example notebook of how to do this soon!
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r""" """)
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
