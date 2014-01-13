import gammalib

def test_a():
    model = gammalib.GModelSpectralGauss()
    energy = gammalib.GEnergy(1, "TeV")
    time = gammalib.GTime()
    print(model.eval(energy, time))
    
def test_mc():
    gammalib.GModelSpectralGauss.mc
    pass

def test_write():
    gammalib.GModelSpectralGauss.write
    pass

def test_flux():
    gammalib.GModelSpectralGauss.flux
    pass


if __name__ == '__main__':
    test_a()