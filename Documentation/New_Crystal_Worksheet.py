import marimo

__generated_with = "0.11.22"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Adding a crystal
        This worksheet guides you through all the things you need to make an entry for the crystal database, then produces a text block you can paste inside.

        To get started, just some imports:
        """
    )
    return


@app.cell
def _():
    import numpy as np
    import LightwaveExplorer as lwe
    return lwe, np


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Next, let's write down the info we'll need to make the database entry
        This is a list of what you need to know about your crystal to describe it to LWE.

        This basic form will assume we can get our refractive index directly from refractiveindex.info; I'll also provide a more advanced example where we can fit to a different set of data.
        """
    )
    return


@app.cell
def _(lwe, np):
    #First, the name of the crystal. I'm putting L afterwards to show that the refractive index is going to be stored in complex valued Lorentzians. A G means that the absorption features are Gaussians.
    Name="My new crystal"

    #Next the crystal type, or number of optical axes.
    # 0: isotropic
    # 1: uniaxial
    # 2: biaxial
    CrystalType = 2

    #Next we set the sellmeier equation to use. 
    # 0: Special fitting form (see paper)
    # 1: Set of Lorentzian oscillators
    # 2: Set of Gaussian bands
    SellmeierEquation=0

    # Maybe we can grab the refractive index equations directly from RefractiveIndex.info... if so, it'll look like this:
    if CrystalType == 0:
        scFit = lwe.getSellmeierFromRII("https://refractiveindex.info/database/data/main/ZnTe/nk/Li.yml")
    if CrystalType == 1:
        scFitO = lwe.getSellmeierFromRII("https://refractiveindex.info/database/data/main/BaB2O4/nk/Zhang-o.yml")
        scFitE = lwe.getSellmeierFromRII("https://refractiveindex.info/database/data/main/BaB2O4/nk/Zhang-e.yml")
    if CrystalType == 2:
        scFitX = lwe.getSellmeierFromRII("https://refractiveindex.info/database/data/main/LiB3O5/nk/Chen-alpha.yml")
        scFitY = lwe.getSellmeierFromRII("https://refractiveindex.info/database/data/main/LiB3O5/nk/Chen-beta.yml")
        scFitZ = lwe.getSellmeierFromRII("https://refractiveindex.info/database/data/main/LiB3O5/nk/Chen-gamma.yml")
    #note that you have to fill out the URLs for your crystal type.
    #If you're using SellmeierEquation either 1 or 2, you'll need to do a fitting rather than directly taking the values from refractiveindex.info
    # in that case you should see the LBO Lorentzian example.
    #To locate the coefficients URL, just go to the page for your crystal and copy the link labeled "Full database record"

    #Don't forget to say where your parameters are from
    sellmeierReference = "C Chen et al.. New nonlinear-optical crystal: LiB3O5, J Opt. Soc. Am. 6, 616-621 (1989) via refractiveindex.info"

    #Now we'll add the chi(2) tensor
    #The type is just 0 or 1, 1 means that chi(2) exists for this crystal
    Chi2Type = 1

    #Next we'll add the tensor. It's helpful to look up the crystal in Boyd's book, where there are tables of the symmetry
    # Here's a paper with values:
    dReference = "S. Lin et al., Journal of Applied Physics 67, 634 (1990)"

    # Fill in your values here, and put the tensor in the correct form
    # Example: the tensor for a crystal in the mm2 group would look like this:
    # these values are for LBO
    d15 = -2.405
    d24= 2.605
    d33 = 0.15 
    dTensor = np.array([
        [ 0.0, 0.0, 0.0, 0.0, d15, 0.0],
        [ 0.0, 0.0, 0.0, d24, 0.0, 0.0],
        [ d15, d24, d33, 0.0, 0.0, 0.0]])

    #do you need to apply any rotation operations? Check now!
    # LBO needs them because the optical axes and the crystal axes have a different ordering, so we do this:
    dTensor = lwe.chi2axisSwap(dTensor,1,3,2)
    #MAKE SURE YOU DISABLE THIS LINE IF YOU'RE ADDING A DIFFERENT CRYSTAL!

    #We should also write down the frequencies used to measure these values so that we can use Miller's rule
    #to approximate their dispersion. They used SHG of 1079 nm light, which we put into frequency (Hz) with
    #freq1 + freq2 = freq3
    Chi2Freq1 = 277.8e12
    Chi2Freq2 = Chi2Freq1
    Chi2Freq3 = Chi2Freq1 + Chi2Freq2

    #Next for chi3
    # in principle, we could add a full chi3 tensor, but it's usually very hard to find them in literature.
    # chi3Type 1: we have the full tensor
    # chi3Type 2: we only have an n2 or chi3 value from e.g. a z-scan, under the assumtption that it's symmetric
    chi3Type = 2

    if chi3Type == 1:
    #if the full tensor is available, that's great! fill it in like this:
    #(this is the BaF2 tensor)
        chi3Tensor = np.array([
            [1.59e-22, 0.00e+00, 0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 1.10e-22, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00],
            [0.00e+00, 1.10e-22, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00, 0.00e+00, 1.59e-22, 0.00e+00, 0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 1.10e-22, 0.00e+00],
            [0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 1.10e-22, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00, 0.00e+00, 1.10e-22, 0.00e+00, 0.00e+00, 0.00e+00, 1.59e-22]
        ])
    elif chi3Type == 2:
        #If they give the value as n2, not chi(3), we have to convert it:
        n2 = 2.24e-20
        n0forn2 = 1.566
        chi3 = n2 * (n0forn2**2) / 283
    elif chi3Type == 0:
        chi3 = 0.0
    else:
        assert False, "chi3 type has to be 0, 1, or 2"

    # Give a reference string for your value:
    Chi3Reference = "B. Maingot, et al., Optics Letters 48, 3243 (2023)"

    #Now for the frequencies; since it's a Kerr effect measurement, they're all the same. Should be a frequency in Hz.
    Chi3Freq1 = 289.9e12
    Chi3Freq2 = Chi3Freq1
    Chi3Freq3 = Chi3Freq1
    Chi3Freq4 = Chi3Freq1
    return (
        Chi2Freq1,
        Chi2Freq2,
        Chi2Freq3,
        Chi2Type,
        Chi3Freq1,
        Chi3Freq2,
        Chi3Freq3,
        Chi3Freq4,
        Chi3Reference,
        CrystalType,
        Name,
        SellmeierEquation,
        chi3,
        chi3Tensor,
        chi3Type,
        d15,
        d24,
        d33,
        dReference,
        dTensor,
        n0forn2,
        n2,
        scFit,
        scFitE,
        scFitO,
        scFitX,
        scFitY,
        scFitZ,
        sellmeierReference,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Generate a block of text that can be pasted inside the CrystalDatabase.txt file""")
    return


@app.cell
def _(
    Chi2Freq1,
    Chi2Freq2,
    Chi2Freq3,
    Chi2Type,
    Chi3Freq1,
    Chi3Freq2,
    Chi3Freq3,
    Chi3Freq4,
    Chi3Reference,
    CrystalType,
    Name,
    SellmeierEquation,
    chi3,
    chi3Tensor,
    chi3Type,
    dReference,
    dTensor,
    lwe,
    scFit,
    scFitE,
    scFitO,
    scFitX,
    scFitY,
    scFitZ,
    sellmeierReference,
):
    print("Name:")
    print(Name)
    print("Type:")
    print(CrystalType)
    print("Sellmeier equation:")
    print(SellmeierEquation)
    if CrystalType == 0:
        print("1st axis coefficients:")
        lwe.printSellmeier(scFit)
        print("2nd axis coefficients:")
        lwe.printSellmeier(scFit)
        print("3rd axis coefficients:")
        lwe.printSellmeier(scFit)
    if CrystalType == 1:
        print("1st axis coefficients:")
        lwe.printSellmeier(scFitO)
        print("2nd axis coefficients:")
        lwe.printSellmeier(scFitE)
        print("3rd axis coefficients:")
        lwe.printSellmeier(scFitE)
    if CrystalType == 2:
        print("1st axis coefficients:")
        lwe.printSellmeier(scFitX)
        print("2nd axis coefficients:")
        lwe.printSellmeier(scFitY)
        print("3rd axis coefficients:")
        lwe.printSellmeier(scFitZ)

    print("Sellmeier reference:")
    print(sellmeierReference)
    print("chi2 type:")
    print(Chi2Type)
    print("d:")
    print('\n'.join(' '.join(map(str, row)) for row in dTensor))
    print("d reference:")
    print(dReference)
    print("chi3 type:")
    print(chi3Type)
    print("chi3:")
    if chi3Type == 1:
        print('\n'.join(' '.join(map(str, row)) for row in chi3Tensor))
    else:
        print(chi3)
        print(0)
        print(0)
    print("chi3 reference:")
    print(Chi3Reference)
    print("Spectral file:")
    print("None")
    print("Nonlinear reference frequencies:")
    print(f"{Chi2Freq1:.8g} {Chi2Freq2:.8g} {Chi2Freq3:.8g} {Chi3Freq1:.8g} {Chi3Freq2:.8g} {Chi3Freq3:.8g} {Chi3Freq4:.8g}")
    print("~~~crystal end~~~\n")
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
